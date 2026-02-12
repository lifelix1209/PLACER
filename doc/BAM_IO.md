# PLACER BAM_IO 模块文档（零拷贝版本）

本文档描述当前代码中的 BAM I/O 实现与调用方式，基于：
- `include/bam_io.h`
- `src/stream/bam_io.cpp`
- `include/pipeline.h`
- `src/pipeline/pipeline.cpp`

---

## 1. 设计目标

当前 BAM_IO 目标是：
1. **单遍流式读取** BAM。
2. **零拷贝传递记录**（传递 `bam1_t` 所有权，而不是解码到 `std::string`/`std::vector`）。
3. 通过 `ReadView` 提供按需访问，避免早期深拷贝。
4. 让上层 pipeline 在真正需要时（如插入片段提取、assembly）再做序列解码，并尽量只做局部解码。

---

## 2. 核心类型

### 2.1 `BamRecordPtr`
定义：
- `using BamRecordPtr = std::unique_ptr<bam1_t, BamRecordDeleter>`

行为：
- 独占所有权，RAII 自动释放。
- 删除器调用 `bam_destroy1`。

### 2.2 `ReadView`
`ReadView` 不拥有数据，只包装 `const bam1_t*`：
- 坐标与质量：`tid()/pos()/mapq()/flag()`
- read 名称：`qname()`（返回 `std::string_view`，生命周期受 `bam1_t` 约束）
- 原始序列访问：`seq_encoded()/seq_len()`
- 原始 CIGAR 访问：`cigar()/n_cigar()`
- tag 检查：`has_tag()/has_sa_tag()/has_md_tag()`
- tag 读取：`get_int_tag()`（例如读取 `NM`）
- 需要时才解码：`decode_sequence()`、`decode_subsequence(start, length)`

这意味着：
- 绝大多数路径可以在压缩/原始表示上决策。
- 只有少量路径需要调用 `decode_sequence()` / `decode_subsequence()` 触发字符串分配。

### 2.3 回调接口
- `RecordHandler = std::function<void(BamRecordPtr&&)>`
- `ProgressHandler = std::function<bool(int64_t processed, int32_t current_tid)>`

---

## 3. Reader 抽象与实现

### 3.1 抽象
`BamStreamReader` 提供：
- `is_valid()`
- `bam_path()`
- `chromosome_count()`
- `chromosome_name(tid)`
- `stream(record_handler, progress_handler, progress_interval)`

工厂函数：
- `make_bam_reader(bam_path, decompression_threads)`

### 3.2 `HtslibBamStreamReader`（当前实现）
构造阶段：
1. `hts_open`
2. 可选 `hts_set_threads`
3. `sam_hdr_read`

`stream()` 主循环：
1. 每轮 `bam_init1()` 新建记录对象（用于安全移动所有权）。
2. `sam_read1` 读入。
3. 过滤：
   - 跳过 `BAM_FSECONDARY`
   - 跳过 `BAM_FUNMAP`
4. 将 `BamRecordPtr` move 给上层回调。
5. 按 `progress_interval` 回调进度。
6. 结束后再触发一次最终进度回调。

返回值：
- 成功：返回处理记录数（过滤后二级计数）。
- reader 无效或 handler 为空：返回 `-1`。

---

## 4. `ReadView` 的零拷贝访问细节

1. `seq_encoded()` 直接返回 `bam_get_seq(record_)` 指针。
2. `cigar()` 直接返回 `bam_get_cigar(record_)` 指针。
3. `decode_sequence()` / `decode_subsequence()` 才进行 4-bit 到字符解码，字符表为 `"=ACMGRSVTWYHKDBN"`。
4. `qname()` 返回 `bam_get_qname(record_)` 的 `std::string_view`，上层如需跨 record 生命周期持有请显式拷贝为 `std::string`。

当前上层实现中：
- `decode_subsequence()` 主要用于 Module 2.1 的插入片段提取（soft-clip / long insertion）。
- `decode_sequence()` 只在 assembly 默认模块中调用（占位实现 + lazy decode）。

---

## 5. 与 Pipeline 的集成模式

当前 pipeline 已从“全量加载”改为“滑动窗口流式处理”：

1. `Pipeline::run()` 按配置选择：
   - `run_streaming()`
   - `run_parallel()`（生产者-消费者）

2. 记录消费入口：`consume_record(BamRecordPtr&&)`
   - Gate1 预过滤：基于 `ReadView`，不解码序列。
   - 维护 `active_window`（`deque<WindowCoord>`）。
   - 线性 bin 分桶（`bin_index = pos / bin_size`）。
   - 发生 bin 切换时立即结算上一 bin。

3. bin 结算入口：`process_bin_records(...)`
   - `component -> insertion_fragment(2.1) -> te_quick_classifier(2.2, optional) -> local_realign -> assembly -> placeability -> genotyping`
   - 只保留轻量结果到 `PipelineResult.final_calls`
   - 不缓存全量 read

---

## 6. 线性 bin（替代 `std::map`）

当前实现不再使用 `std::map` 做分桶聚合。

算法：
1. 维护 `current_tid` 与 `current_bin_index`。
2. 新 read 到达后计算新 `bin_index`。
3. 若与当前 bin 不同，则 flush 当前 bin 并开始新 bin。
4. 同 bin 继续 append。

复杂度：
- O(N) 单次线性扫描（依赖 BAM 的 coordinate-sorted 特性）。

---

## 7. 并行模式（生产者-消费者）

`PipelineConfig.enable_parallel=true` 时启用。

模型：
1. 生产者（reader 线程）读取 BAM，按 `batch_size` 组 batch。
2. batch 推入线程安全 `SafeQueue<std::vector<BamRecordPtr>>`。
3. 消费者 worker 线程拉取 batch，执行与 streaming 模式相同的 `consume_record` 和 flush。

特点：
- 解耦 I/O 与业务计算。
- 保持单 worker 顺序消费，避免跨 batch 重排序问题。

---

## 8. 默认 Gate1 预过滤（当前实现）

当前默认模块：`SignalFirstGate1Module`（`src/gate1/gate1_module.cpp`）。

规则为“信号优先，质量兜底”：
1. 硬过滤：
   - `unmapped` 丢弃
   - `secondary` 丢弃
   - `seq_len < 50` 丢弃
2. 信号判定（任一满足即视为有信号）：
   - `supplementary` flag
   - `SA` tag
   - 长 soft-clip（默认 `S >= 100`）
3. 若有信号：忽略 MapQ，但加三条保险丝：
   - 保险丝 1：至少一侧参考 anchor 像样（`max_match_block >= 200`）
   - 保险丝 2：若是长 soft-clip，clip 邻近 flank 需满足最小匹配锚长度（默认 `>=120`）
   - 保险丝 3：若有 `NM` tag，`NM / total_match_bases` 过高（默认 `>0.20`）则丢弃
4. 若无信号：仅保留 `MAPQ > 20` 的背景 reads。

该逻辑全部基于 `ReadView` 的零拷贝访问，不会触发序列解码。
说明：当前 pipeline 的 Module 2.1 会对少量候选片段调用 `decode_subsequence()` 做局部解码；assembly 默认模块仍会对极少量 reads 调用 `decode_sequence()`（占位实现）。

---

## 9. 资源与内存行为

1. BAM_IO 层不构造 `ReadSketch` 大对象。
2. 主路径持有的是 `BamRecordPtr`（原始记录）。
3. `active_window` 仅存 `(tid,pos)` 轻量坐标。
4. bin flush 后对应 `BamRecordPtr` 批次释放。
5. 避免了“全量通过 reads 存入 vector”内存线性增长问题。

---

## 10. 当前边界

1. 并行模式当前为单消费者 worker（1 producer + 1 consumer）。
2. SA 解析未在 BAM_IO 层展开为完整结构，仅通过 `has_sa_tag` 做快速判断。
3. 默认模块仍是占位实现，重点是数据流与接口架构，不是最终算法质量。

---

## 11. 运行方式

串行流式：
- `./build/placer <input.bam> <ref.fa> [te.fa]`

并行模式：
- `PLACER_PARALLEL=1 ./build/placer <input.bam> <ref.fa> [te.fa]`

输出：
- `scientific.txt`（汇总 + 每个 call 的轻量结果）
