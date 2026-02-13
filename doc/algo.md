# PLACER 项目结构与当前算法逻辑（按现有代码）

本文档只描述当前仓库里已经实现并会执行的逻辑，基于以下主干代码：

- `src/main.cpp`
- `src/pipeline/pipeline.cpp`
- `src/gate1/gate1_module.cpp`
- `src/component/insert_fragment_module.cpp`
- `src/component/te_quick_classifier.cpp`
- `src/stream/bam_io.cpp`
- `include/pipeline.h`
- `include/gate1_module.h`
- `include/bam_io.h`

---

## 1. 项目结构

### 1.1 顶层目录

- `src/`：核心实现
- `include/`：对外头文件与跨模块数据结构
- `doc/`：设计与实现说明
- `test_data/`：本地样例数据
- `third_party/abPOA/`：第三方依赖代码（当前主流程未直接调用）
- `build/`：CMake 构建产物

### 1.2 `src/` 模块划分

- `src/main.cpp`
  - 参数解析
  - 从环境变量读取并行开关
  - 构建并运行 `Pipeline`
  - 输出 `scientific.txt`
- `src/stream/bam_io.cpp`
  - htslib 流式读取 BAM
  - `ReadView` 提供零拷贝字段访问与按需序列解码
- `src/gate1/gate1_module.cpp`
  - Gate1 预筛模块：信号优先 + 三道保险丝
- `src/pipeline/pipeline.cpp`
  - 串流/并行调度
  - 按 `bin` 聚合并驱动各阶段
  - 组装最终 `FinalCall`
- `src/component/insert_fragment_module.cpp`
  - 从组件内 read 提取插入相关片段（clip / I / SA）
  - 可选写入 `ins_fragments.fasta`
- `src/component/te_quick_classifier.cpp`
  - 基于 TE FASTA 构建 unique k-mer 索引
  - 对片段做快速分类并进行组件级投票
  - 可选写入 `ins_fragment_hits.tsv`

### 1.3 `include/` 关键接口

- `include/bam_io.h`
  - `BamStreamReader`、`ReadView`、`BamRecordPtr`
- `include/pipeline.h`
  - 核心数据结构：`ComponentCall`、`InsertionFragment`、`FragmentTEHit`、`FinalCall` 等
  - 模块接口：`IGate1Module` 到 `ITEQuickClassifierModule`
  - 默认模块类声明与 `PipelineBuilder`
- `include/gate1_module.h`
  - `Gate1SignalConfig` 与 `SignalFirstGate1Module`

---

## 2. 端到端执行流程

入口 `main()` 调用 `Pipeline::run()`，主流程如下：

1. 打开 BAM，流式读取每条对齐记录（跳过 unmapped / secondary）。
2. 对每条记录执行 Gate1 预筛。
3. 通过 Gate1 的记录按坐标分桶（`bin_size`，默认 10kb）。
4. 每个 bin 触发一次组件处理（当前实现为“一个 bin 一个 component”）。
5. 对每个 component 顺序执行：
   - 插入片段提取（Module 2.1）
   - TE 快速分类与投票（Module 2.2，可选）
   - 局部证据收集（SimpleLocalRealign）
   - 组装（LazyDecodeAssembly）
   - placeability 打分
   - 基因分型
6. 汇总为 `FinalCall` 列表，最后写出 `scientific.txt`。

---

## 3. 调度与数据流细节

### 3.1 两种运行模式

- 串流模式（默认）：`run_streaming()`
- 并行模式：设置环境变量 `PLACER_PARALLEL=1` 后走 `run_parallel()`
  - 生产者线程：BAM 读取并按 `batch_size` 入队
  - 消费者线程：顺序消费 batch，执行同样的 `consume_record()`

并行模式当前是 1 个 worker，不是多 worker fan-out。

### 3.2 `consume_record()` 行为

对每条记录：

1. `gate1_module_->pass_preliminary()` 判断是否保留。
2. 维护一个 `active_window`（按 `window_size`，默认 10kb，当前主要用于窗口状态维护）。
3. 计算 `bin_index = pos / bin_size`。
4. 染色体或 bin 切换时，先 `flush_current_bin()`。
5. 把记录移动进当前 bin 缓冲。

### 3.3 `process_bin_records()` 行为

1. 生成 `bin_start/bin_end/chrom`。
2. `component_module_->build(...)` 输出 `ComponentCall` 列表。
3. 对每个 component：
   - `ins_fragment_module_->extract(...)`
   - 若 TE 分类启用：
     - `te_classifier_module_->classify(...)`
     - `te_classifier_module_->vote_cluster(...)`
   - `local_realign_module_->collect(...)`
   - `assembly_module_->assemble(...)`
   - `placeability_module_->score(...)`
   - `genotyping_module_->genotype(...)`
4. 组装 `FinalCall`：
   - TE 投票 `passed` 时才写 `te_name`
   - 其余字段总是写入（tier/support/gt/af/gq 等）

---

## 4. 各阶段算法逻辑

## 4.1 BAM I/O（`src/stream/bam_io.cpp`）

- 使用 htslib：`hts_open + sam_hdr_read + sam_read1`
- 过滤规则：
  - 丢弃 `BAM_FUNMAP`
  - 丢弃 `BAM_FSECONDARY`
- `ReadView` 提供：
  - 直接访问 `tid/pos/mapq/flag/cigar/tag`
  - `decode_subsequence(start,len)` 局部解码
  - `decode_sequence()` 全长解码（仅在 assembly 占位模块使用）

这是“零拷贝 + 延迟解码”路线。

## 4.2 Gate1（`SignalFirstGate1Module`）

参数默认值见 `Gate1SignalConfig`：

- `min_seq_len = 50`
- `long_soft_clip_min = 100`
- `background_mapq_min = 20`
- `min_anchor_match_bases = 200`
- `min_clip_flank_match_bases = 120`
- `max_nm_rate = 0.20`

逻辑：

1. 硬过滤：unmapped / secondary / 序列太短直接拒绝。
2. 信号判定：满足以下任一即视为有信号：
   - supplementary flag
   - SA tag
   - 长 soft-clip（`max_soft_clip >= 100`）
3. 无信号时采用背景保留：`MAPQ > 20` 才保留。
4. 有信号时应用三道保险丝：
   - Fuse1：`max_match_block >= 200`
   - Fuse2：长 clip 侧必须有足够邻接匹配锚（`>=120`）
   - Fuse3：若有 NM tag，`NM / total_match_bases <= 0.20`

## 4.3 Component（`LinearBinComponentModule`）

当前实现是 MVP：

- 每个 bin 只产生 1 个 component。
- `anchor_pos = (bin_start + bin_end) / 2`。
- `read_indices` 包含该 bin 全部 read 索引。

同时做三类 read 标记：

- split/SA 类：有 SA tag 或 supplementary。
- soft-clip 类：该 read 的最大 soft-clip `>= 20`。
- long insertion 类：该 read 的最大 CIGAR `I >= 50`。

`breakpoint_candidates` 结构已预留，但当前未填充。

## 4.4 插入片段提取（Module 2.1）

实现：`CigarInsertionFragmentModule + SplitSAFragmentModule`

输入：一个 `ComponentCall` + 对应 `bin_records`  
输出：`std::vector<InsertionFragment>`

片段来源：

1. soft-clip 片段
   - 阈值：`min_soft_clip_for_seq_extract`（默认 50）
   - leading/trailing clip 分别提取
2. 长插入片段（CIGAR `I`）
   - 阈值：`min_long_ins_for_seq_extract`（默认 50）
   - 每条 read 最多保留长度最大的 2 个 `I`
3. SA 片段（MVP）
   - 解析 `SA:Z`
   - 仅主比对记录解析（跳过 supplementary/secondary）
   - 要求 SA CIGAR 能推断 query 区间且长度满足阈值
   - 每条 read 最多 `max_sa_per_read`（默认 3）

输出副作用：

- 若 `ins_fragments_fasta_path` 非空（默认 `ins_fragments.fasta`）：
  - 构造时清空文件
  - 运行时追加写入 FASTA
  - 使用互斥锁保证并行下写文件安全

## 4.5 TE 快速分类（Module 2.2）

实现：`TEKmerQuickClassifierModule`

### 索引构建

- 输入：`te_fasta_path`（为空则模块禁用）
- k-mer 长度：`te_kmer_size`（默认 15）
- 对每条 TE 序列及其反向互补都建索引。
- 仅保留由 `A/C/G/T` 构成的 k-mer。
- 索引语义：
  - k-mer 只出现在一个 TE：记录该 TE id
  - 出现在多个 TE：标记为歧义（`-1`）

### 片段分类

对每个 `InsertionFragment`：

1. 扫片段所有有效 k-mer，统计各 TE 命中次数。
2. 选命中数最多的 TE 作为 best。
3. 计算：
   - `hit_kmers`
   - `total_kmers`
   - `coverage = hit_kmers / total_kmers`
   - `kmer_support = coverage`
   - `aligned_len_est`（best TE 连续命中 run 转换）

### component 投票

- 仅对 `te_name` 非空命中参与投票。
- `fragment_count < te_min_fragments_for_vote`（默认 3）直接失败。
- 票数最高 TE 为 `best_te`。
- `vote_fraction = best_votes / fragment_count`
- `median_identity = median(kmer_support(best_te_hits))`
- 通过条件（默认）：
  - `vote_fraction >= 0.60`
  - `median_identity >= 0.60`

输出副作用：

- 若 `ins_fragment_hits_tsv_path` 非空（默认 `ins_fragment_hits.tsv`）：
  - 构造时清空并写表头
  - 分类时追加结果

## 4.6 局部证据（`SimpleLocalRealignModule`）

当前是占位评分，不做真实局部比对：

- 以 bin 中心为参考点 `center`。
- 每条 read 生成一条 `LocusEvidence`。
- 分数：
  - `normalized_score = max(0, 1 - |read.pos - center| / 10000)`

## 4.7 Assembly（`LazyDecodeAssemblyModule`）

当前是占位组装：

- 取 component 第一条 read。
- 解码整条 read 序列（lazy decode，仅此处解码全长）。
- 截断到 300bp 作为 `consensus`。
- `pos` 使用 component `anchor_pos`。

## 4.8 Placeability（`SimplePlaceabilityModule`）

- `support_reads = evidence.size()`
- `delta_score = mean(normalized_score) * 100`
- tier 划分：
  - `delta_score >= 30` -> Tier 1
  - `10 <= delta_score < 30` -> Tier 2
  - 其他 -> Tier 3

## 4.9 Genotyping（`SimpleGenotypingModule`）

当前是经验阈值分型：

- `af = clamp(support_reads / 20.0, 0, 1)`
- `gq = round(delta_score)`
- `af >= 0.8` -> `1/1`
- `0.2 <= af < 0.8` -> `0/1`
- `af < 0.2` -> `0/0`

---

## 5. 关键配置与当前默认值

来自 `PipelineConfig`（`include/pipeline.h`）：

- `bam_threads = 2`
- `window_size = 10000`
- `bin_size = 10000`
- `enable_parallel = false`
- `batch_size = 1000`
- `ins_fragments_fasta_path = "ins_fragments.fasta"`
- `min_soft_clip_for_seq_extract = 50`
- `min_long_ins_for_seq_extract = 50`
- `min_sa_aln_len_for_seq_extract = 50`
- `max_sa_per_read = 3`
- `ins_fragment_hits_tsv_path = "ins_fragment_hits.tsv"`
- `te_kmer_size = 15`
- `te_vote_fraction_min = 0.60`
- `te_median_identity_min = 0.60`
- `te_min_fragments_for_vote = 3`

---

## 6. 输出与字段

主程序固定写 `scientific.txt`，包含：

- 汇总计数：
  - `total_reads`
  - `gate1_passed`
  - `processed_bins`
  - `components`
  - `evidence_rows`
  - `assemblies`
  - `placeability_calls`
  - `genotype_calls`
- 每条 final call：
  - `chrom/tid/pos/window_start/window_end`
  - `te/te_vote_frac/te_median_ident/te_fragments`
  - `tier/support_reads/gt/af/gq`

---

## 7. 当前实现边界（与长期目标的差异）

1. `reference_fasta_path` 目前在主流程中基本未使用（保留为后续扩展位）。
2. `Component` 仍是“每 bin 一个组件”，未做真正断点聚类。
3. `LocalRealign`/`Assembly`/`Genotyping` 目前均为占位或简化模型。
4. TE 分类是快速 unique k-mer 估计，不是全比对 identity。

以上意味着当前代码重点是“流式架构可跑通 + 模块接口清晰”，并非最终科研精度版本。
