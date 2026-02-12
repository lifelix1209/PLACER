# Component 模块（Gate1 后候选聚合 + 2.1/2.2）

本文档描述 **当前代码已实现** 的 Component 阶段逻辑与接口，重点是：

- Gate1 之后如何把 reads 聚合成 component（当前实现为按 bin 聚合的占位版本）
- Module 2.1：从 component 内 reads 提取插入片段序列（soft-clip / long insertion）
- Module 2.2：对插入片段做快速 TE 分类（unique k-mer）并对 component 投票

对应代码位置：

- 类型与接口：`include/pipeline.h`
- Component（当前实现）：`src/pipeline/pipeline.cpp` (`LinearBinComponentModule`)
- 2.1 插入片段提取：`src/component/insert_fragment_module.cpp` (`CigarInsertionFragmentModule`)
- 2.2 TE quick 分类：`src/component/te_quick_classifier.cpp` (`TEKmerQuickClassifierModule`)

---

## 1. 输入与输出

### 1.1 输入

- Pipeline 以流式方式读取 BAM，并在 `Pipeline::consume_record(...)` 中按 `bin_size` 聚合出一个 bin 的 `bin_records`（`std::vector<BamRecordPtr>`）。
- `BamRecordPtr` 是 `bam1_t` 的所有权封装；使用 `ReadView` 对字段做零拷贝访问。

### 1.2 输出

Component 阶段输出的是 `std::vector<ComponentCall>`，其中每个 `ComponentCall` 至少包含：

- `chrom`, `tid`
- `bin_start`, `bin_end`
- `anchor_pos`
- `read_indices`
- 三类 read 索引（用于后续按证据类型处理）：
  - `soft_clip_read_indices`
  - `split_sa_read_indices`
  - `insertion_read_indices`
- `breakpoint_candidates`

说明：

- `BreakpointCandidate`/`breakpoint_candidates` 当前主要作为后续“断点聚类”版本的接口预留；在当前实现里不会被填充。

---

## 2. 当前 Component 构建实现：`LinearBinComponentModule`

### 2.1 聚合策略（当前占位）

当前 `LinearBinComponentModule::build(...)` 的策略是：

- **每个 bin 仅生成 1 个 `ComponentCall`**（也就是“bin == component”的占位实现）
- `anchor_pos = (bin_start + bin_end) / 2`
- `read_indices = [0, 1, ..., bin_records.size()-1]`

这保证了后续 pipeline 的调用链可用，但尚未实现“断点定位 + 窗口聚类”。

### 2.2 三类证据的粗分类（当前实现）

对 bin 内每条 read（索引为 `i`）：

Split / SA / Supplementary：

- 若 `ReadView::has_sa_tag()` 为真，或 BAM flag 包含 `BAM_FSUPPLEMENTARY`，则将 `i` 加入 `split_sa_read_indices`。

Soft-clip：

- 扫描 CIGAR，计算该 read 的 `max_soft_clip`。
- 若 `max_soft_clip >= 20`，则将 `i` 加入 `soft_clip_read_indices`。

Long insertion：

- 扫描 CIGAR，计算该 read 的 `max_insertion`。
- 若 `max_insertion >= 50`，则将 `i` 加入 `insertion_read_indices`。

注意：

- 这些阈值（20/50）当前在 `src/pipeline/pipeline.cpp` 中写死，后续应迁移到配置（例如 `ComponentConfig` 或 `PipelineConfig`）。
- 一条 read 可以同时落入多个类别。

### 2.3 复杂度

- 主要成本是对每条 read 扫一次 CIGAR：`O(n_reads * n_cigar)`。

---

## 3. Module 2.1：插入片段提取（`CigarInsertionFragmentModule` + `SplitSAFragmentModule`）

该模块在 `Pipeline::process_bin_records(...)` 中，对每个 component 执行片段提取：

- 输出结构：`std::vector<InsertionFragment>`
- 可选输出文件：`PipelineConfig.ins_fragments_fasta_path`（默认 `ins_fragments.fasta`）
- 默认会合并两类来源：CIGAR（soft-clip / long insertion）和 SA split 片段

### 3.1 启动时的文件行为

若 `ins_fragments_fasta_path` 非空：

- 构造函数会以 `std::ios::trunc` 清空该文件（便于每次运行生成干净输出）
- 每次 `extract(...)` 会先缓冲本批 FASTA 文本，再在互斥锁保护下以 `std::ios::app` 追加写入

将 `ins_fragments_fasta_path` 置空可禁用文件写入，但仍会返回 `InsertionFragment` 向量。

### 3.2 Soft-clip 片段提取

参数：

- `PipelineConfig.min_soft_clip_for_seq_extract`（默认 50）

实现细节：

- 仅检查 leading/trailing soft-clip（跳过 hard-clip，取首个/最后一个非 hard-clip CIGAR op 是否为 `S`）。
- leading clip：`start=0, length=leading_soft_clip`
- trailing clip：`start=seq_len - trailing_soft_clip, length=trailing_soft_clip`
- 片段序列通过 `ReadView::decode_subsequence(start, length)` 延迟解码获得。
- `InsertionFragment.source` 采用**参考侧语义**：
  - leading（首个非 hard-clip 端）-> `clipRefLeft`
  - trailing（最后一个非 hard-clip 端）-> `clipRefRight`
  - 不随 `BAM_FREVERSE` 变换（CIGAR 端点本身已是参考侧定义）
- `InsertionFragment.is_reverse` 记录 BAM `BAM_FREVERSE`。

### 3.3 Long insertion（CIGAR `I`）片段提取

参数：

- `PipelineConfig.min_long_ins_for_seq_extract`（默认 50）

实现细节：

- 扫描 CIGAR 时维护 query 游标 `qpos`（凡是消耗 query 的 op 都会推进 `qpos`）。
- 遇到 `I` 且 `len >= min_long_ins_for_seq_extract`，记录 `{start=qpos, len}`。
- 若一条 read 有多个长插入事件，当前实现仅保留 **长度最大的前 2 个**（避免片段爆炸）。
- 对每个事件使用 `ReadView::decode_subsequence(start, len)` 提取插入序列。

### 3.4 Split-SA 片段提取（MVP）

参数：

- `PipelineConfig.min_sa_aln_len_for_seq_extract`（默认 50）
- `PipelineConfig.max_sa_per_read`（默认 3）

实现细节：

- 仅在 primary record 上解析 `SA:Z`（跳过 supplementary/secondary，避免重复片段）
- 解析每条 SA 的 `CIGAR`，取 query 区间：
  - `qstart = leading soft-clip(S) length`
  - `qlen = sum(M,=,X,I)`（不把 `S` 放进片段长度）
  - 若 SA CIGAR 两端都没有 `S`，MVP 直接跳过（无法仅凭 SA 稳定定位 query 起点）
- 通过 `ReadView::decode_subsequence(qstart, qlen)` 切出片段，`source=kSplitSa`
- 每条 read 最多输出 `max_sa_per_read` 条 SA 片段，且要求 `qlen >= min_sa_aln_len_for_seq_extract`

### 3.5 Fragment ID 与 FASTA 格式

- `fragment_id` 带必要的溯源信息，例如：
  - `chrom:anchor_pos|read=<qname>|idx=<idx>|strand=+|src=clipRefL|len=...`
  - `chrom:anchor_pos|read=<qname>|idx=<idx>|strand=-|src=ins|start=...|len=...`
  - `chrom:anchor_pos|read=<qname>|idx=<idx>|strand=+|src=splitSA|sa=chr1:12345+|qlen=...`
- `fragment_id` 保留 `|` 分隔结构，便于后续解析字段。
- FASTA 输出时只对 header 做安全化（空白/换行替换），不破坏 `fragment_id` 的字段分隔。
- FASTA 序列按 80 列折行写出。

### 3.6 当前边界

- 对 SA CIGAR 为纯 `M`（无软剪） 的情况，MVP 仍会跳过；后续可结合 supplementary 记录补强。

---

## 4. Module 2.2：快速 TE 分类（unique k-mer）与 component 投票

该模块用于对 2.1 产出的插入片段做“快速但粗糙”的 TE 分类，占位用于打通 pipeline。

### 4.1 索引构建（TE FASTA -> unique k-mer）

参数：

- `PipelineConfig.te_fasta_path`（为空则禁用该模块）
- `PipelineConfig.te_kmer_size`（默认 15）

索引构建规则：

- 从 TE FASTA 逐条读取：
  - TE 名称取 header 的第一个 token（按空白分割）
  - 序列转为大写
- 对每条 TE 序列，分别对正向序列与反向互补序列提取 k-mer（都映射到同一个 TE id）：
  - 仅接受完全由 `A/C/G/T` 组成的 k-mer（含 N 等字符会跳过）
  - 若 k-mer 第一次出现：记录到该 TE id
  - 若 k-mer 已存在但来自不同 TE：标记为“歧义”（值置为 `-1`）

查询时：

- 只有映射到 **唯一 TE** 的 k-mer（值 `>=0`）才会参与计数
- 歧义 k-mer（`-1`）会被忽略

### 4.2 片段分类（per-fragment）

对每条 `InsertionFragment`：

- 统计该片段内所有“有效 k-mer”的数量 `total_kmers`
- 统计每个 TE 的 unique k-mer 命中次数，选择命中最多的 TE 作为 `best`：
  - `hit_kmers = best_hits`
  - `coverage = best_hits / total_kmers`（快速估计）
  - `kmer_support = coverage`（当前实现直接复用 coverage，不代表序列 identity）
- `aligned_len_est`：对 `best` TE 做二次扫描，取“连续命中 best 的 k-mer 的最长 run”：
  - `aligned_len_est = k + max_run - 1`（若 `max_run==0` 则为 0）

可选输出文件：

- `PipelineConfig.ins_fragment_hits_tsv_path`（默认 `ins_fragment_hits.tsv`）
- 构造函数会清空并写入表头；每次 `classify(...)` 会追加写入记录行：
  - `fragment_id, te, frag_len, hit_kmers, total_kmers, coverage, aligned_len_est, kmer_support`

### 4.3 component 投票（per-component）

投票输入是该 component 的所有 `FragmentTEHit`：

- 先过滤掉 `te_name` 为空的 hit（表示该片段没有命中任何 TE）
- 要求 `fragment_count >= PipelineConfig.te_min_fragments_for_vote`（默认 3），否则直接不通过
- 统计每个 `te_name` 的票数，选择票数最多的 `best_te`
- `vote_fraction = best_votes / fragment_count`
- `median_identity`：对 `best_te` 对应的 hits 的 `kmer_support` 取中位数
- `passed` 条件（默认）：
  - `vote_fraction >= PipelineConfig.te_vote_fraction_min`（0.60）
  - 且 `median_identity >= PipelineConfig.te_median_identity_min`（0.60）

---

## 5. 与 pipeline 的集成位置

集成入口在 `Pipeline::process_bin_records(...)`（`src/pipeline/pipeline.cpp`）：

1. `component_module_->build(...)` 生成 `ComponentCall`
2. 对每个 component：
   - 2.1 `ins_fragment_module_->extract(...)` 得到 `InsertionFragment` 列表
   - 若 2.2 启用（存在 TE k-mer index）：
     - `te_classifier_module_->classify(...)`
     - `te_classifier_module_->vote_cluster(...)`
   - TE 投票结果写入最终的 `FinalCall.te_*` 字段
3. 后续继续走 `LocalRealign -> Assembly -> Placeability -> Genotyping`

---

## 6. TODO（后续计划，不是当前实现）

以下内容是你提出的“断点定位 + 窗口聚类”的目标方向，但 **目前代码尚未实现**：

- 从三类证据（soft-clip / SA/supp / long I）提取参考侧断点 `BreakpointCandidate`
- 按 `cluster_gap` 对断点做窗口聚类，将一个 bin 内的 reads 收敛成多个 `ComponentCall`
- 输出约束（`min_support`, `MAD/IQR` 等）用于“入口宽，输出严”
- 从 split/SA/supp 里更精确地提取“段间未对齐区段”作为插入片段
