# PLACER 当前算法流程（与代码一致）

本文档按当前实现整理，代码基线对应：

- `src/main.cpp`
- `src/pipeline/pipeline.cpp`
- `src/pipeline/decision_policy.cpp`
- `src/gate1/gate1_module.cpp`
- `src/component/insert_fragment_module.cpp`
- `src/component/te_quick_classifier.cpp`
- `src/component/te_consensus_module.cpp`
- `src/component/tsd_detector.cpp`
- `include/pipeline.h`

## 1. 入口与运行模式

命令行入口：

```bash
./build/placer <input.bam> <ref.fa> [te.fa]
```

`main()` 做两件事：

1. 解析输入文件并构造 `PipelineConfig`。
2. 读取大量 `PLACER_*` 环境变量覆盖阈值（TE、assembly、TSD、evidence、genotype、bootstrap、并行等）。

随后执行：

- `build_default_pipeline(config)->run()`
- 输出 `scientific.txt`

运行模式：

- `PLACER_PARALLEL=0/未设置`：`run_streaming()`
- `PLACER_PARALLEL=1`：`run_parallel()`（生产者 + worker 队列）

## 2. 全流程总览

每个 bin/component 的核心链路固定为：

1. Gate1 读段预筛
2. 证据点聚窗并构建 component
3. 插入片段提取（CIGAR + split-SA）
4. TE 快速分类（multi-k）
5. Anchor-locked TE 共识（可选）
6. 证据特征计算与硬过滤
7. abPOA 局部组装
8. 组装后 TE 决策（classic/rescue）
9. placeability 打分与 tier
10. 插入接受决策（高/低置信）
11. 分型（0/0,0/1,1/1）
12. TE 置信校准 + TSD 检测 + open-set 状态分类
13. 可选 bootstrap 导出
14. final call 去重与汇总

## 3. Gate1（`SignalFirstGate1Module`）

硬过滤：

- `unmapped` / `secondary` 直接丢弃
- `seq_len < min_seq_len` 丢弃

信号定义：满足任一

- supplementary
- 有 `SA` tag
- 长 soft-clip

若无信号：仅保留 `MAPQ > background_mapq_min`。

若有信号，需通过 3 个 fuse：

1. 连续锚定匹配长度充足（`max_match_block`）
2. 长软剪切两侧锚长度足够
3. 若有 `NM`，要求 `NM / total_match_bases <= max_nm_rate`

## 4. 组件构建（`LinearBinComponentModule`）

不是“每 bin 一个 component”，而是：

1. 从读段抽取证据点（soft-clip/indel/SA hint）并赋权。
2. 对证据点做直方图平滑，找密度峰。
3. 从峰扩展候选窗口并合并相邻窗口。
4. 把 read 按“窗口得分”分配到最佳窗口（含歧义判定）。
5. 每个窗口生成一个 `ComponentCall`，并填充：
   - `read_indices`
   - `soft_clip_read_indices`
   - `split_sa_read_indices`
   - `insertion_read_indices`
   - `breakpoint_candidates`

另外支持跨 bin 上下文（`kCrossBinContextBins=2`），避免长读跨桶丢证据。

## 5. 插入片段提取（Module 2.1）

### 5.1 CIGAR 路径（`CigarInsertionFragmentModule`）

提取三类片段：

- 左/右 soft-clip（`kClipRefLeft/Right`）
- 长插入 `I`（`kCigarInsertion`，每读最多保留最长 2 个）

`short_ins_enable` 时会在 split/indel 主导场景下动态放宽最小插入长度阈值。

### 5.2 Split-SA 路径（`SplitSAFragmentModule`）

两步输出：

1. 稳健路径：
   - 解析 primary + supplementary + `SA:Z`
   - 如有显式 supplementary，对 SA 做一致性校验（query 区间 IOU/边界差）
   - 选 flank（锚长度与 NM/大 indel 惩罚联合排序）
   - 选 mate，定位 query 连接点 `q_junc`
   - 截取 opposite 片段并推断 `ref_junc_pos`/`ref_side`
2. 诊断回退：
   - 若稳健路径失败，保留有限 SA 诊断片段（`max_sa_per_read`）

### 5.3 文件输出

默认会写：

- `ins_fragments.fasta`

模块构造时清空，运行中追加（并发写有互斥锁）。

## 6. TE 快速分类（Module 2.2）

实现：`TEKmerQuickClassifierModule`

### 6.1 索引

- 从 `te.fa` 读取模板
- 支持多 k（`te_kmer_sizes_csv` + `te_kmer_size`）
- 每个 k 建 unique-kmer 索引（歧义 kmer 记为 `-1`）

### 6.2 单片段分类

对每个片段：

1. 低复杂 soft-clip 先过滤（AT 比例/同聚物阈值）。
2. 在多个 k 上计算加权支持，得到 top1/top2 及 `multik_support`。
3. 估计 `aligned_len_est`（最长连续命中 run）。
4. 可触发低 k-mer rescue：
   - 当 top1 支持不足或 top1-top2 间隔过小
   - 对 topN TE 做半全局编辑身份估计
   - 超过 rescue identity 阈值则改判并打 `rescue_used`

输出行可写入：

- `ins_fragment_hits.tsv`

### 6.3 component 投票

`vote_cluster()` 输出：

- `top1/top2`
- `posterior_top1/posterior_top2/posterior_margin`
- `te_name`（满足最小 fragment 数后）
- `vote_fraction`、`median_identity`
- `multik_support`、`rescue_frac`

`passed` 条件：

- `vote_fraction >= te_vote_fraction_min`
- `median_identity >= te_median_identity_min`

## 7. Anchor-locked TE 共识（Module 2.3）

实现：`DeterministicAnchorLockedModule`

核心思想：在目标 TE 模板上对每个片段枚举多个放置，做确定性收敛，输出稳定的 TE 内断点窗口。

步骤：

1. 目标 TE 选择
   - 若 `te_call.passed` 用其 `te_name`
   - 否则按 hits 多数票选
2. 枚举 placement
   - seed k-mer 生成起点候选
   - 对每个起点做半全局仿射对齐（成本 `norm_cost`）
3. CoreCandidate 过滤
   - 需满足 ref side、锚长度、起点不确定性、`delta_tpl`、split-SA 侧向约束等
4. 方向判别 `theta`（FWD/REV）
   - 基于 `phi()` + weighted median + MAD 比较
   - 可能标记 `fail_theta_uncertain`
5. 一轮重分配
   - 目标函数：`norm_cost + λ|phi-c0| + μ*span_start`
6. 生成核心窗口
   - 取 core set
   - 输出 `te_breakpoint_core`、`[win_start, win_end]`
   - 统计 `core_candidate_count/core_set_count/split_sa_core_frac/ref_junc_pos_min/max`

## 8. 证据建模与硬过滤（`collect_evidence`）

汇总维度：

- 覆盖与支持：`global_cov/local_cov/depth/alt/ref`
- 断点离散度：`breakpoint_mad`
- 证据来源计数：`split/clip/indel`
- 低复杂 soft-clip 比例
- TE 一致性（弱门）

关键过滤条件：

- `alt_support_reads >= min_support_required`
- `breakpoint_mad <= evidence_breakpoint_mad_max`
- `low_complex_softclip_frac <= evidence_low_complex_softclip_frac_max`
- `pass_te_consistency`（Stage-A 弱门，允许非强 TE 标签）

任一失败则 `hard_filtered=true`，该 component 终止。

## 9. 组装（`assemble_component`）

- 过滤短片段后，使用 abPOA 组装 consensus
- `identity_est` 用 k-mer overlap 估计
- 依据场景（split/indel dominant、短插入）可放宽 identity 最低阈值

组装 QC 典型标签：

- `NO_CANDIDATE_FRAGMENTS`
- `INSUFFICIENT_POA_READS`
- `CONSENSUS_TOO_SHORT`
- `LOW_IDENTITY`
- `PASS`

## 10. 组装后 TE 决策（`evaluate_post_assembly_te_decision`）

主要分支：

1. `PASS_CLASSIC`：标准 vote + identity 通过
2. `PASS_ASM_RESCUE`：标准失败但 assembly rescue 通过
3. `PASS_SPLIT_INDEL_RESCUE`：无 soft-clip 或 split/indel 主导时的救援通路
4. `PASS_SHORT_INS_RESCUE`：短插入救援
5. pure soft-clip 专门门槛（读段数、片段数、identity）

失败会给出明确 `FAIL_*` QC。

## 11. Placeability + 插入接受 + 分型

### 11.1 Placeability（`score_placeability`）

用 logistic 组合证据特征得到 `evidence_p` 和 `evidence_q`，并分 tier：

- Tier1: `p >= evidence_tier1_prob`
- Tier2: `evidence_tier2_prob <= p < tier1`
- Tier3: 其他

### 11.2 插入接受（`evaluate_insertion_acceptance`）

- Tier1/2 默认接受
- Tier3 仅在 `emit_low_confidence_calls=true` 且满足低置信门槛时接受

输出：`pass` + `low_confidence`。

### 11.3 分型（`genotype_call`）

- 由 alt/ref 深度计算 AF
- 二项似然比较 `0/0, 0/1, 1/1`
- 取最优基因型并给 `GQ`

## 12. TSD、TE open-set 与最终状态

### 12.1 TE 置信度校准

`calibrate_te_confidence_prob()` 结合：

- posterior
- assembly identity
- support
- vote/TE identity
- breakpoint MAD
- rescue 惩罚

输出 `te_conf_prob`。

### 12.2 TSD 检测（`TSDDetector`）

- 从参考序列提取左右断点邻域
- 优先尝试 DUP（左右同序列）
- 再尝试 DEL（断点间缺失序列）
- 计算背景出现率 `bg_p`

结果写入：`tsd_type/tsd_len/tsd_seq/tsd_bg_p`，并对 `te_conf_prob` 做小幅调节。

若开启 `te_fail_on_tsd_inconsistent`，会追加 `FAIL_TSD_INCONSISTENT` 并压低置信度。

### 12.3 Open-set 分类（`classify_te_open_set`）

根据 TE gate、proxy 信号与 `te_conf_prob`，输出：

- `TE_CERTAIN`
- `TE_UNCERTAIN`
- `NON_TE`

必要时提升 `top1` 为 `te_name`，并追加 proxy/QC 标签。

## 13. Bootstrap 导出（可选）

开启 `bootstrap_export_enable` 后，对 `TE_UNCERTAIN`（以及可选 `NON_TE`）导出：

- 共识 FASTA：`pass1_bootstrap_consensus.fasta`
- 元数据 TSV：`pass1_bootstrap_calls.tsv`

用于二次迭代建库。

## 14. FinalCall 去重与输出

`finalize_final_calls()` 会按位置/TE family 去重（默认 50bp 邻域），优先保留支持更强、tier 更高、置信更高、组装更好的 call。

最终输出在 `scientific.txt`：

1. 运行汇总统计
2. 每条 call 的全字段明细（位点、TE 诊断、TSD、置信、tier、GT/AF/GQ、assembly 指标等）

## 15. 当前默认关键参数（节选）

- `bin_size=10000`, `window_size=10000`
- `te_kmer_size=13`, `te_kmer_sizes_csv=9,11,13`
- `te_vote_fraction_min=0.40`, `te_median_identity_min=0.30`
- `te_rescue_vote_fraction_min=0.25`, `te_rescue_median_identity_min=0.20`
- `short_ins_enable=true`, `short_ins_min_len=35`, `short_ins_max_len=300`
- `assembly_min_fragment_len=80`, `assembly_min_consensus_len=80`, `assembly_min_identity_est=0.55`
- `tsd_enable=true`, `tsd_min_len=3`, `tsd_max_len=50`, `tsd_bg_p_max=0.05`
- `emit_low_confidence_calls=false`

如需调参，优先使用 `main.cpp` 支持的 `PLACER_*` 环境变量。
