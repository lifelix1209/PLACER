# PLACER 当前算法逻辑说明（基于现有代码）

> 这份文档只描述“现在代码里真正执行的逻辑”，不写理想化流程。  
> 入口代码：`src/main.cpp`。

---

## 1. 总览：主流程到底在做什么

当前主程序不是纯 streaming 模式，而是：

1. 先把 BAM 读成内存里的 `std::vector<ReadSketch>`。
2. 然后按 phase 串行执行：
   - Phase 2: Gate1 过滤
   - Phase 3: 组件构建（Component）
   - Phase 4: 局部重比对（LocusEvidence）
   - Phase 5: 组装（Contig + 结构折叠）
   - Phase 8: TE reverse index 召回（可选）
   - Phase 6: Placeability
   - Phase 7: Genotyping
   - Enrichment: TE 注释 + TSD + MAPQ
   - TE-specific tier 重打分
   - 写出 `scientific.txt`

主执行顺序在 `PLACERPipeline::run()`，文件：`src/main.cpp`。

---

## 2. 核心数据对象与不变量

### 2.1 `ReadSketch`（`include/bam_reader.h`）
- 是全流程 read 级基础对象。
- 包含：`tid/pos/end_pos`、`cigar_ops`、`sequence`、`SA`、`mapq`、mate 信息（`mtid/mpos/isize/proper_pair`）等。
- `breakpoint` 字段类型 `ReadBreakpoint`，默认构造里 `read_pos/ref_pos` 为 `-1`，`is_valid=false`。

### 2.2 `Component`（`include/component_builder.h`）
- 一个局部断点簇。
- 包含 `read_indices`、`start/end/centroid`、`locus_set` 以及锚点统计。
- 后续 Phase 4/5 都围绕 component 处理。

### 2.3 `LocusEvidence`（`include/local_realign.h`）
- `read_idx × locus_pos` 级证据。
- 含 `normalized_score`、`weighted_identity`、`evidence_bits` 等。

### 2.4 `Contig` 与 `StructuralRepresentative`（`include/assembly.h`）
- `Contig` 是 Phase 5 的直接产物，包含：
  - `sequence`（`up+ins+down` 的 POA 共识）
  - `ins_seq`（仅插入片段集合的共识）
  - `up_flank_seq/down_flank_seq`
- `StructuralRepresentative` 是结构折叠后代表，后续分型、输出都用它。

### 2.5 `LocusCallBundle`（`src/main.cpp`）
- 主流程内部汇总对象：把 representative / evidence / placeability / genotype / TE注释 / TSD /候选统计绑在一起。

---

## 3. 各模块算法细节

## 3.1 BAM 读取（Phase 1）

文件：`src/stream/bam_reader.cpp`

逻辑：
- 逐条读 BAM，填充 `ReadSketch`。
- 解析 CIGAR：
  - 计算 `end_pos`
  - 统计 `total_clip_len`
  - 标记 `has_large_insertion`（`I > 50bp`）
- 抽取 read 序列、SA/MD tag 信息。

注意：
- Phase 1 并不做候选筛选，只做数据标准化和轻量特征预提取。

---

## 3.2 Gate1 过滤（Phase 2）

文件：`src/gate1/gate1.cpp`，调用在 `src/main.cpp` 的 `gate1_filter()`

### 核心方法
- 从 read 抽 probe：
  - 两端 probe（默认 `probe_len=200`）
  - CIGAR 驱动 probe：长 soft-clip、长 insertion 周边
- probe 用 `HashTEIndex` 做 k-mer 命中统计。
- 判定通过条件：
  - 命中满足 `min_hit_count/min_hit_density`，或
  - 总命中超过 `min_total_hits`。

### 当前实现特性
- 有 SA 的 read 直接放行。
- `filtered_reads` 和 `filtered_probes` 强制索引对齐（后面 reverse index 依赖这个对齐）。

---

## 3.3 Component 构建（Phase 3）

文件：`src/component/component_builder.cpp`

### Step A: 多源 Anchor 提取
- 信号来源：
  - clip 锚点（5'/3'）
  - 大插入 `I`
  - 大缺失 `D`
  - SA 锚点（含跨染色体补充锚点）

### Step B: read 内合并
- `anchor_merge_distance` 内合并；
- 优先级 `SA > CLIP > INS > DEL`，位置做均值更新。

### Step C: 密度聚类 + 递归 breaker
- 先按 `cluster_gap` 做线性聚类；
- 簇跨度太大时（`max_cluster_span`）触发递归 breaker：
  - 优先找“密度谷”，回退最大 gap；
  - 子簇比例要达 `min_split_ratio`。

### Step D: 构建 `Component`
- 计算 `centroid`、`density`、`read_count`、证据计数与 `evidence_score`。

---

## 3.4 局部重比对（Phase 4）

文件：`src/local_realign/local_realign.cpp`

主调用：`realign_and_collect(component, reads, genome)`

### Step A: 生成候选 loci
- 从 component 的 read 集合中收集：
  - read `pos`
  - read `end_pos`
  - SA 位置（同染色体）
- 在 `merge_distance=50bp` 内合并，按强度排序，限制 top-k。

### Step B: 证据生成（read × locus）
- 尝试真实序列比对（局部 affine gap）：
  - 从参考取 `locus ± flank_length`
  - 从 read 中取对应片段
  - 算 `normalized_score`、identity、indel
- 比对失败则退化为距离分：
  - `1 - dist/search_window`

### Step C: 聚合成 `PlaceabilityReport`
- 对每个 locus 聚合分数，算 best/second、delta、confidence。
- 同时输出 strand 统计和 tier 初判（1~5）。

当前状态：
- 仍有“无法比对时用距离分回退”的路径，不是纯序列对齐路径。

---

## 3.5 Assembly（Phase 5）

文件：`src/assembly/assembly.cpp`

这是当前候选质量的关键阶段。

### Step A: 断点识别
- `detect_breakpoint_from_cigar()` 从 CIGAR 推断：
  - `INSERTION`（记录 event-time `read_pos/ref_pos/insertion_len`）
  - `SPLIT`（N 操作）
  - `SOFT_CLIP_5P/3P`

### Step B: 按断点抽取 `ReadSegments`
- 统一得到 `up / ins / down`：
  - `SOFT_CLIP_5P`: `ins = [0, bp_pos)`
  - `SOFT_CLIP_3P`: `ins = [bp_pos, read_end)`
  - `INSERTION`: `ins = [bp_pos, bp_pos+insertion_len)`
  - `SPLIT`: 当前直接跳过（不抽 ins）
- 强约束：
  - `ins_len >= ins_min_length`
  - 至少一侧 flank 足够长
- 然后保留“主导断点类型”，避免混合不同语义窗口。

### Step C: 分桶判定（可选）
- 先构造 analysis window：`up_tail + ins + down_head`。
- minimizer Jaccard 统计：
  - split 触发条件：`median<=0.75 && p10<=0.60 && n>=6`
- split 后做自检：
  - `mean_aa >= mean_ab + 0.08`
  - `mean_bb >= mean_ab + 0.08`
- 自检不通过则回退单群组装。

### Step D: 共识构建
- 对每个桶（或全体）做 POA。
- 输出两个层面的序列：
  - `contig.sequence`：`up+ins+down` 的共识
  - `contig.ins_seq`：只用 `ins` 片段做共识（后续 TE-only 更依赖这个）
- `ins_seq` 为空时直接丢弃该 contig（严格路径）。

### Step E: 结构折叠
- 用 fingerprint + RTree + DSU 合并结构近似 contig。
- `StructuralFingerprint` 包含 `tid/breakpoint/ins_length_range/family/orientation`。

---

## 3.6 TE Reverse Index 召回（Phase 8，可选）

文件：`src/te_reverse_index/te_reverse_index.cpp`

启用条件：
- 只有传了 TE fasta，main 才会 `te_reverse_enabled=true`。

逻辑：
1. 对无 SA 的 read（或传入 probes）做 probe→genome k-mer 查询。
2. hit 聚类成 `RescuedLocus`（按 `tid+position` 半径聚类）。
3. 打分 `placeability_score` 并过滤（`min_locus_reads`、`min_hit_density`）。
4. 与已有 component 近邻整合，补充 evidence / representative。

---

## 3.7 Placeability（Phase 6）

文件：`src/placeability/placeability.cpp`

在主流程中：
- 优先使用 Phase 4 已算好的 `PlaceabilityReport`。
- 仅在缺失时 fallback 到 `PlaceabilityScorer::calculate_full_placeability()`。

`PlaceabilityScorer` 的关键点：
- 聚合同 read+locus 证据，避免重复计数。
- 侧翼一致性分两层：
  - `side_consistent`
  - `side_consistency_verified`（是否真的有 read 位置信息支撑）
- `determine_tier()` 对 Tier1 要求 side consistency 已验证。

---

## 3.8 Genotyping（Phase 7）

文件：`src/genotyping/genotyping.cpp`

### 两层模型
1. read-level 三成分混合：ALT / REF / NULL。
2. genotype-level：根据 `E_alt/E_ref` 计算 `0/0,0/1,1/1` 似然与 GQ。

### EM 关键实现
- E-step 用标准形式：
  - `r_ik ∝ pi_k * likelihood_k(x_i) * prior_ik`
- M-step 更新 `pi_alt/pi_ref/pi_null`。
- log-likelihood 与 E-step 同公式，避免目标不一致。
- 输出：
  - `mix_*`（混合权重）
  - `e_*_sum`
  - `AF = E_alt / (E_alt + E_ref)`
  - LRT 显著性、GQ、GT。

---

## 3.9 Enrichment：TE 注释、TSD、读段窗口统计

入口：`annotate_bundles()`，文件：`src/main.cpp`  
TE 注释引擎：`src/annotation/te_annotator.cpp`

### TE 注释
- 先从每个候选 read 抽 query（优先 breakpoint 对应片段）。
- `TEAnnotator::annotate_consensus()`：
  1. 每条 query 做 TE 对齐（正反向 SW）
  2. 过滤低质量命中（`aligned_bases>=50 && identity>=70`）
  3. 选 dominant family + dominant strand
  4. 按长度分布去离群
  5. 仅对 TE 对齐区间做 consensus
- 结果写到 bundle：
  - `family/subfamily/startTE/endTE/strand`
  - `te_match/query_coverage/aligned_bases`
  - `te_consensus`

### TSD 与空位统计
- 先推断断点区间 `bp_left/bp_right`。
- 在 `wiggle=200bp` 内统计：
  - `span_reads_wiggle`
  - `nonref_reads_wiggle`
  - `empty_reads_wiggle`
- TSD 优先级：
  1. representative flanks
  2. 窗口读段提取 flank 共识
  3. 单 read fallback
- `TSDDetector` 检测失败时再用 flank fallback 匹配。

---

## 3.10 TE-specific tier 重打分

入口：`apply_te_specific_tier()`，文件：`src/main.cpp`

机制：
- 先取 base tier。
- 基于 TE 命中质量、coverage、aligned_bases、TSD、双侧 anchor、强信号、candidate_confidence、nonref fraction 等累计 `te_points`。
- 再将 tier 上下调（最多 ±2 级）。
- 同时生成 `te_quality_flags`（例如 `NoTEAnnotation/NoTSD/LowSupport`）。

---

## 3.11 最终输出 TXT 与去重

入口：`write_txt()` + `TXTWriter`，文件：`src/main.cpp`

### 候选硬过滤（写出前）
- 必须有 representative + evidence + placeability。
- 最小 evidence read 数、最小 signal read 数、最小 confidence。
- 要求双侧 anchor 或强 split/indel/discordant 支持。

### 两层去重
1. 事件级去重（坐标、长度带、family、strand 约束）。
2. 序列级去重（同事件内 TE 共识相似度，k-mer Jaccard）。

### TE-only 输出策略（当前实现）
- 若 `has_te_annotation`：
  - `TE_Consensus` 写 TE 共识；
  - `LengthIns` 优先取 TE 共识长度。
- 若无 TE 注释：
  - `TE_Consensus=NA`
  - `INS_Extended=NA`
  - `LengthIns=0`

也就是说，最终主输出优先是 TE-only 语义，不再直接把长窗口序列当最终 TE 共识。

---

## 4. 配置默认值（主流程关注）

定义在 `src/main.cpp::PipelineConfig`：

- Assembly:
  - `assembly_min_reads=3`
  - `assembly_max_reads=50`
  - `assembly_ins_min_length=10`（传给 AssemblyConfig 后用于段长约束）
  - `assembly_soft_clip_core_len=800`
- Placeability:
  - `delta_tier1=30`
  - `delta_tier2=10`
- Genotyping:
  - `max_em_iterations=20`
  - `spatial_lambda=100`

注：某些模块内部还有自己阈值（例如 TEAnnotator 的 `aligned_bases>=50`、`identity>=70`），不完全由 `PipelineConfig` 控制。

---

## 5. 现在代码的已知边界（不是猜测，是实现事实）

1. 主程序当前是“先全量读 BAM 再分 phase”，`stream_processor` 没走主路径。
2. Assembly 中 `SPLIT` 类型目前不产出 `ins` 片段，只作为证据信号。
3. LocalRealigner 仍保留距离回退路径，不是 100% 序列对齐闭环。
4. Genotyping 里的 family 特征目前偏弱：`ReadSketch` 没有 read 级 family，`compute_family_features()` 多数情况下是未知态。
5. TE reverse index 只在给 TE fasta 时启用；不给 TE fasta 时不会做召回。

---

## 6. 调试入口（建议）

- `PLACER_DEBUG_ASM=1`：看 assembly 分段、trim、split 统计。
- `PLACER_DEBUG_ANNOT=1`：看 TE 注释支持数、输出长度、family。
- 重点日志：
  - `[AssemblySplitStats] ... mean_aa/mean_bb/mean_ab ...`
  - `TE annotated loci / TSD-positive loci / te consensus filled`

---

## 7. 代码导航索引（按阶段）

- 主流程编排：`src/main.cpp`
- BAM 读取：`src/stream/bam_reader.cpp`
- Gate1：`src/gate1/gate1.cpp`
- Component：`src/component/component_builder.cpp`
- Local realign：`src/local_realign/local_realign.cpp`
- Assembly：`src/assembly/assembly.cpp`
- Reverse index：`src/te_reverse_index/te_reverse_index.cpp`
- TE 注释：`src/annotation/te_annotator.cpp`
- Placeability：`src/placeability/placeability.cpp`
- Genotyping：`src/genotyping/genotyping.cpp`

