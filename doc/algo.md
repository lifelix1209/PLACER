# PLACER 项目思路与实现（按当前代码）

本文档基于仓库当前可执行主干实现整理，目标是说明两件事：

1. 项目想解决什么问题（设计思路）。
2. 代码现在具体是怎么做的（实现细节与边界）。

对应核心代码：

- `src/main.cpp`
- `src/stream/bam_io.cpp`
- `src/gate1/gate1_module.cpp`
- `src/pipeline/pipeline.cpp`
- `src/component/insert_fragment_module.cpp`
- `src/component/te_quick_classifier.cpp`
- `src/component/te_consensus_module.cpp`
- `include/pipeline.h`
- `include/bam_io.h`
- `include/gate1_module.h`

---

## 1. 项目思路（Why）

PLACER 的目标是从长读长 BAM（ONT/PacBio）里检测非参考 TE 插入，核心难点来自：

- 重复区域导致定位歧义（同一片段可落多个位置）。
- 对齐器可能把插入“线性化”到 CIGAR，或只留软剪切 / SA 线索。
- 直接做全局重组装成本高，且容易把邻近位点或不同单倍型混在一起。

当前实现采用“分层降本 + 流式处理”的思路：

- 先用 Gate1 做信号优先预筛，尽量减少后续处理量。
- 以 bin 为处理单位串流推进，不把 BAM 全量读入内存。
- 先提取可解释的插入片段，再做 TE 快速分类（k-mer 投票）。
- 在可用时做 anchor-locked 的确定性 TE 断点共识，输出诊断字段。
- local realign / assembly / placeability / genotyping 目前是可运行占位版本，逻辑已内置到 Pipeline 固定流程。

---

## 2. 总体架构（How）

### 2.1 执行入口

`main()`（`src/main.cpp`）流程：

1. 读取命令行参数：`placer <input.bam> <ref.fa> [te.fa]`。
2. 构造 `PipelineConfig`（含 BAM、REF、TE 路径）。
3. 读取环境变量 `PLACER_PARALLEL` 决定串流或并行模式。
4. `build_default_pipeline(config)->run()` 执行主流程。
5. 把结果写到当前目录 `scientific.txt`。

说明：

- `reference_fasta_path` 当前只作为配置输入保留，主流程里尚未实质使用。
- `te_fasta_path` 直接影响 TE 快速分类与 TE 共识模块是否可用。

### 2.2 固定流程模块

当前 `Pipeline` 不再通过 `I*` 抽象接口注入模块，而是固定成员模块 + 固定步骤函数：

- Gate1: `SignalFirstGate1Module`
- Component: `LinearBinComponentModule`
- Fragment 提取: `CigarInsertionFragmentModule`（内部再调用 `SplitSAFragmentModule`）
- TE 快速分类: `TEKmerQuickClassifierModule`
- TE 共识: `DeterministicAnchorLockedModule`
- 证据收集/组装/placeability/分型：`Pipeline` 内部函数
  - `collect_evidence(...)`
  - `assemble_component(...)`
  - `score_placeability(...)`
  - `genotype_call(...)`

`build_default_pipeline(config)` 当前只负责创建 BAM reader 并构造固定流程 `Pipeline`。

---

## 3. 主流程数据通路

### 3.1 BAM 串流与分桶

`Pipeline::run()` 按 `config.enable_parallel` 选择：

- `run_streaming()`：单线程串流。
- `run_parallel()`：生产者（读 BAM）+ 1 个消费者 worker（处理 batch）。

并行模式说明：

- 目前是“单 worker 队列”，不是多 worker 扇出。
- `batch_size` 默认 1000。

### 3.2 `consume_record()` 行为

每条记录依次执行：

1. Gate1 预筛，不通过直接丢弃。
2. 维护 `active_window`（`window_size` 默认 10kb）滑窗状态。
3. 以 `bin_index = pos / bin_size`（默认 10kb）分桶。
4. 染色体或 bin 变化时触发 `flush_current_bin()`。
5. 将记录放入 `current_bin_records`。

说明：

- `active_window` 当前仅维护窗口状态，不直接参与后续打分。

### 3.3 `process_bin_records()` 行为

对每个非空 bin：

1. 组件构建：`component_module_.build(...)`
2. 对每个 component 依次执行：
   - 片段提取：`ins_fragment_module_.extract(...)`
   - TE 分类：`te_classifier_module_.classify(...) + vote_cluster(...)`（启用时）
   - TE 共识：`anchor_locked_module_.resolve(...)`（启用时）
   - 证据收集：`collect_evidence(...)`
   - 组装：`assemble_component(...)`
   - placeability：`score_placeability(...)`
   - 分型：`genotype_call(...)`
3. 组装 `FinalCall` 并追加到 `result.final_calls`。

TE 名称写入规则：

- 若 `te_call.passed`，使用投票结果 `te_call.te_name`。
- 否则若 `anchor_report.has_result` 且有 `anchor_report.te_name`，回退使用共识模块给的 TE 名称。

---

## 4. 各模块实现细节

## 4.1 BAM I/O（`src/stream/bam_io.cpp`）

实现类：`HtslibBamStreamReader`

- 使用 `hts_open`、`sam_hdr_read`、`sam_read1` 流式读取。
- 读阶段先过滤：
  - `BAM_FUNMAP`（未比对）
  - `BAM_FSECONDARY`（次要比对）
- 支持 htslib 解压线程数 `bam_threads`（默认 2）。

`ReadView` 提供：

- 轻量访问 `tid/pos/mapq/flag/cigar/tag`。
- `decode_subsequence(start, len)`：按需局部解码序列。
- `decode_sequence()`：全长解码（仅 assembly 占位模块在用）。

---

## 4.2 Gate1（`SignalFirstGate1Module`）

默认参数（`Gate1SignalConfig`）：

- `min_seq_len = 50`
- `long_soft_clip_min = 100`
- `background_mapq_min = 20`
- `min_anchor_match_bases = 200`
- `min_clip_flank_match_bases = 120`
- `max_nm_rate = 0.20`

判定逻辑：

1. 硬过滤：unmapped / secondary / read 长度不足，直接拒绝。
2. 信号判断（满足其一）：
   - supplementary flag
   - 存在 `SA` tag
   - 最大 soft-clip `>= long_soft_clip_min`
3. 若无信号：仅保留 `MAPQ > background_mapq_min`。
4. 若有信号：通过三道保险丝：
   - Fuse1：最大连续 match 块 `>= min_anchor_match_bases`
   - Fuse2：长 soft-clip 侧邻接锚长度必须足够
   - Fuse3：若有 `NM` tag，要求 `NM / total_match_bases <= max_nm_rate`

---

## 4.3 Component（`LinearBinComponentModule`）

当前实现是占位版本：

- 一个 bin 生成一个 `ComponentCall`。
- `anchor_pos = (bin_start + bin_end) / 2`。
- `read_indices` 包含该 bin 全部读段索引。

同时做三类读段标签：

- `split_sa_read_indices`：有 SA 或 supplementary。
- `soft_clip_read_indices`：最大 soft-clip `>= 20`。
- `insertion_read_indices`：最大 CIGAR 插入 `I >= 50`。

说明：

- `breakpoint_candidates` 结构已预留，当前未填充。

---

## 4.4 插入片段提取（Module 2.1）

主实现：`CigarInsertionFragmentModule::extract()`
补充实现：`SplitSAFragmentModule::extract()`

最终输出：`std::vector<InsertionFragment>`

### 4.4.1 CIGAR 路径（soft-clip + long insertion）

阈值来自 `PipelineConfig`：

- `min_soft_clip_for_seq_extract`（默认 50）
- `min_long_ins_for_seq_extract`（默认 50）

处理规则：

1. 软剪切片段：
   - leading soft-clip -> `source = kClipRefLeft`
   - trailing soft-clip -> `source = kClipRefRight`
2. 长插入片段：
   - 扫 CIGAR 找 `I >= min_long_ins`。
   - 每条 read 按长度降序最多保留 2 个插入片段。
3. 记录 `class_mask`：
   - `kCandidateSoftClip`
   - `kCandidateSplitSaSupplementary`
   - `kCandidateLongInsertion`
4. 片段序列通过 `decode_subsequence()` 延迟提取。

### 4.4.2 Split-SA 路径（稳健提取 + 诊断回退）

核心参数：

- `min_sa_aln_len_for_seq_extract`（默认 50）
- `max_sa_per_read`（默认 3）

内部固定常量（当前写死）：

- `w_ref = max(200, bin_size/2)`
- `w_anchor = 150`
- `mh_max = 20`
- `l_left = l_right = max(min_len, 100)`
- `trim_q = 12`
- `eps_bp = 5`
- `alpha = 1.0`

核心步骤：

1. 先收集同 `read_id` 的显式 supplementary 记录，标准化为 query/ref 区间。
2. 对 primary 记录解析 `SA:Z`，从 SA CIGAR 反推出区间。
3. 若存在显式 supplementary，则做一致性校验：
   - 方向一致，且 query span IOU 高、起止差值小（代码阈值：IOU `>= 0.9`，起止差 `<= 5`）。
4. 在与断点邻域同染色体的对齐里选 `flank`，评分近似为：
   - `score = anchor_len - alpha * (max(0, NM) + 5 * has_large_indel_near_bp)`
5. 选 query 非重叠最大的 `mate`，推断连接点 `q_junc`，围绕连接点截取 opposite 片段。
6. 由 `flank` 方向和 `q_junc` 推断 `ref_junc_pos`，再映射 `ref_side`（Left/Right/Unknown）。
7. 若稳健路径失败，仍可输出 SA 诊断片段（上限 `max_sa_per_read`）。

### 4.4.3 文件输出副作用

若 `ins_fragments_fasta_path` 非空（默认 `ins_fragments.fasta`）：

- 模块构造时清空文件。
- 运行时追加写入片段 FASTA（80 列换行）。
- 通过全局互斥锁保证并行写文件安全。

---

## 4.5 TE 快速分类（Module 2.2）

实现：`TEKmerQuickClassifierModule`

### 4.5.1 索引构建

前提：

- `te_fasta_path` 非空且索引构建成功，否则模块禁用。

规则：

- k-mer 长度 `te_kmer_size`（默认 15）。
- 每条 TE 序列使用“正向 + 反向互补”共同建索引。
- 只接受 `A/C/G/T`。
- k-mer 仅命中单一 TE 时记该 TE id；命中多个 TE 时标为歧义（`-1`）。

### 4.5.2 片段打分

对每个 `InsertionFragment`：

1. 统计 `total_kmers`。
2. 统计各 TE unique k-mer 命中数，取最高 TE 为 `best_te`。
3. 输出：
   - `hit_kmers`
   - `coverage = hit_kmers / total_kmers`
   - `kmer_support = coverage`
   - `aligned_len_est`（best_te 连续命中 run 换算长度）

### 4.5.3 component 投票

规则：

- 仅 `te_name` 非空片段参与投票。
- 若 `fragment_count < te_min_fragments_for_vote`（默认 2）则直接失败。
- 票数最多 TE 为 `te_name`。
- `vote_fraction = best_votes / fragment_count`。
- `median_identity = median(kmer_support of best_te)`。
- 通过条件（默认）：
  - `vote_fraction >= 0.60`
  - `median_identity >= 0.50`

TSV 输出：

- `ins_fragment_hits_tsv_path` 非空时，写 `ins_fragment_hits.tsv`（构造时清空并写表头，运行中追加）。

---

## 4.6 Anchor-locked TE 共识（Module 2.3）

实现：`DeterministicAnchorLockedModule`

目标：在已知 TE 模板下，把“片段落在 TE 内多个可能位置”的不确定性收敛成稳定核心断点窗口，并给出方向诊断 `theta`。

### 4.6.1 启用条件与模板库

- 需 `te_consensus_enable = true` 且 `te_fasta_path` 可成功构建模板库。
- 模板库记录每个 TE 的序列与 seed k-mer 到位置的倒排表（`te_consensus_seed_k`，默认 13）。

### 4.6.2 目标 TE 选择

优先级：

1. 若 TE 快速投票已通过，用 `te_call.te_name`。
2. 否则从 `hits` 里多数票选一个 TE（平票取字典序较小）。

### 4.6.3 片段放置候选 `P_i` 枚举

对每个命中目标 TE 的片段：

1. 基于 seed 命中计数枚举起点（最多 `te_consensus_max_start_candidates`，默认 5）。
2. 对每个起点执行半全局仿射比对 `semiglobal_affine_align_suffix(query, te_suffix)`：
   - mismatch=1, gap_open=2, gap_extend=1。
   - 输出归一化代价 `norm_cost`。
3. 得到若干 `Placement{te_start, te_end, norm_cost}`，按 `norm_cost` 升序排序。

### 4.6.4 CoreCandidate 判定

片段成为 core 候选需同时满足：

- `ref_side != Unknown`
- `anchor_len >= te_consensus_anchor_min`（默认 80）
- `start_span_best <= start_span_bias + start_span_sqrt_scale * sqrt(fragment_len)`
- `delta_tpl >= te_consensus_delta_tpl_min`（默认 0.05）
- 对 split-SA 片段：`|ref_junc_pos - anchor_pos| <= te_consensus_w_side_max`（默认 300）

### 4.6.5 方向判定 `theta0`

定义 `phi(p, side, theta)`：

- `theta = FWD`：
  - `RefLeft -> te_start`
  - `RefRight -> te_end`
- `theta = REV`：
  - `RefLeft -> te_end`
  - `RefRight -> te_start`

对 core 候选按权重：

- `w_i = max(1, anchor_len_i) * exp(-norm_cost_i / te_consensus_temp_t)`

分别计算 FWD/REV 的加权中位中心与 MAD，选择更优方向。

不确定性判定：

- 若 `|mad_fwd - mad_rev| < te_consensus_eps_theta` 且两侧权重和相近（`te_consensus_rel_sumw_eps`），标记 `FAIL_THETA_UNCERTAIN`，输出 `theta0 = Unknown`。
- 代码仍会继续给出核心窗口结果，但 QC 会标记为不确定。

### 4.6.6 一轮确定性重分配 + 核心窗口

1. 先得到 `c0`（core 候选下的加权中位中心）。
2. 每片段在其候选放置里选最优：
   - `obj = norm_cost + lambda*|phi-c0| + mu*|te_start - best_te_start|`
3. 基于重分配结果计算 `c1`。
4. CoreSet 条件：
   - `|phi - c1| <= te_consensus_w_core_gate`（默认 80）。
5. 核心断点与窗口：
   - `te_breakpoint_core = median(phi in CoreSet)`
   - `core_span = 1.4826 * MAD`
   - `w = max(te_consensus_w_min, te_consensus_gamma * core_span)`
   - 输出 `[core-w, core+w]`
6. `has_result` 需要 `core_set_count >= te_consensus_min_core_set`（默认 2）。

输出诊断还包含：

- `mad_fwd/mad_rev`
- `split_sa_core_frac`
- `ref_junc_pos_min/max`
- `core_candidate_count/core_set_count`

---

## 4.7 Local realign（占位）

实现：`SimpleLocalRealignModule`

- 对 component 内每条 read 产出一条 `LocusEvidence`。
- 分数：
  - `normalized_score = max(0, 1 - |pos-center| / 10000)`
- 当前不做真实重比对，仅作为后续模块占位输入。

---

## 4.8 Assembly（占位）

实现：`LazyDecodeAssemblyModule`

- 取 component 第一条 read。
- 解码全长序列并截断到 300bp，作为 `consensus`。
- `pos` 使用 `component.anchor_pos`。

---

## 4.9 Placeability（占位）

实现：`SimplePlaceabilityModule`

- `support_reads = evidence.size()`
- `delta_score = mean(normalized_score) * 100`
- tier 分级：
  - Tier1: `delta_score >= 30`
  - Tier2: `10 <= delta_score < 30`
  - Tier3: 其他

---

## 4.10 Genotyping（占位）

实现：`SimpleGenotypingModule`

- `af = clamp(support_reads / 20.0, 0, 1)`
- `gq = round(delta_score)`
- 分型：
  - `af >= 0.8 -> 1/1`
  - `0.2 <= af < 0.8 -> 0/1`
  - `af < 0.2 -> 0/0`

---

## 5. 输出格式与统计

结果由 `main.cpp` 写入 `scientific.txt`：

1. 汇总计数：

- `total_reads`
- `gate1_passed`
- `processed_bins`
- `components`
- `evidence_rows`
- `assemblies`
- `placeability_calls`
- `genotype_calls`

2. 每条 `FinalCall` 明细字段：

- 位置窗：`chrom tid pos window_start window_end`
- TE 相关：`te te_vote_frac te_median_ident te_fragments`
- TE 共识诊断：`te_theta te_mad_fwd te_mad_rev te_bp_core te_bp_win_start te_bp_win_end te_core_candidates te_core_set split_sa_core_frac te_ref_junc_min te_ref_junc_max te_qc`
- 调用与分型：`tier support_reads gt af gq`

`te_qc` 取值逻辑：

- `DISABLED`：共识模块未启用
- `FAIL_THETA_UNCERTAIN`：方向判定不稳定
- `NO_CORE_RESULT`：未形成有效核心集合
- `PASS`：核心结果可用

---

## 6. 关键默认参数（`PipelineConfig`）

基础调度：

- `bam_threads = 2`
- `progress_interval = 100000`
- `window_size = 10000`
- `bin_size = 10000`
- `enable_parallel = false`
- `batch_size = 1000`

片段提取：

- `ins_fragments_fasta_path = "ins_fragments.fasta"`
- `min_soft_clip_for_seq_extract = 50`
- `min_long_ins_for_seq_extract = 50`
- `min_sa_aln_len_for_seq_extract = 50`
- `max_sa_per_read = 3`

TE 快速分类：

- `ins_fragment_hits_tsv_path = "ins_fragment_hits.tsv"`
- `te_kmer_size = 15`
- `te_vote_fraction_min = 0.60`
- `te_median_identity_min = 0.50`
- `te_min_fragments_for_vote = 2`

TE 共识：

- `te_consensus_enable = true`
- `te_consensus_seed_k = 13`
- `te_consensus_max_start_candidates = 5`
- `te_consensus_delta_tpl_min = 0.05`
- `te_consensus_anchor_min = 80`
- `te_consensus_start_span_sqrt_scale = 2.0`
- `te_consensus_start_span_bias = 6.0`
- `te_consensus_w_side_max = 300`
- `te_consensus_temp_t = 0.25`
- `te_consensus_beta = 0.20`
- `te_consensus_lambda = 0.35`
- `te_consensus_mu = 0.10`
- `te_consensus_w_core_gate = 80.0`
- `te_consensus_w_min = 25.0`
- `te_consensus_gamma = 2.0`
- `te_consensus_eps_theta = 0.05`
- `te_consensus_rel_sumw_eps = 0.10`
- `te_consensus_min_core_set = 2`

---

## 7. 当前实现边界与下一步

当前代码已实现：

- BAM 串流调度、Gate1 预筛、分桶处理、片段提取、TE 快速分类、确定性 anchor-locked 共识、最终文本输出。

当前仍是占位/简化部分：

- Component 仍是“每 bin 一个 component”，未做真实断点聚类。
- local realign / assembly / placeability / genotyping 仍是轻量可运行版本。
- `reference_fasta_path` 尚未进入实质计算链路。
- `third_party/abPOA` 已 vendored，但当前主流程没有直接调用。

因此，PLACER 当前阶段是“可流式运行并输出可诊断结果”的工程骨架版本；高精度科研模型（真实局部重比对、组装、概率分型）可在现有流程骨架上逐步替换升级。
