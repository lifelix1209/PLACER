# PLACER BAM I/O 技术文档（按当前实现）

本文档对应当前代码实现：

- `include/bam_io.h`
- `src/stream/bam_io.cpp`
- （调用侧）`src/pipeline/pipeline.cpp`

目标是准确描述 BAM I/O 层当前“已经实现并实际生效”的行为。

---

## 1. 模块职责

BAM I/O 模块负责：

1. 打开 BAM 并读取 header。
2. 以流式方式逐条输出 `bam1_t` 记录。
3. 在 I/O 层提前过滤 `unmapped` 和 `secondary`。
4. 提供 `ReadView` 作为轻量访问包装，支持按需序列解码。

该模块不负责：

- 区域检索（region fetch）。
- 断点/TE 业务判定。
- 全量缓存所有 reads。

---

## 2. 核心数据类型

## 2.1 `BamRecordPtr`

定义：

- `using BamRecordPtr = std::unique_ptr<bam1_t, BamRecordDeleter>`

语义：

- 独占 `bam1_t` 所有权。
- 释放时调用 `bam_destroy1`（RAII）。
- 在 pipeline 中通过 move 传递，避免深拷贝记录。

## 2.2 `ReadView`

`ReadView` 只包装 `const bam1_t*`，不拥有底层内存。

主要接口：

- 坐标/质量：`tid()`、`pos()`、`mapq()`、`flag()`
- 序列/CIGAR 原始访问：`seq_len()`、`seq_encoded()`、`n_cigar()`、`cigar()`
- 名称：`qname()`（`std::string_view`）
- tag 相关：`has_tag()`、`has_sa_tag()`、`has_md_tag()`、`get_int_tag()`、`get_string_tag()`
- 解码：`decode_subsequence(start, length)`、`decode_sequence()`

生命周期注意：

- `qname()` 返回的是对 `bam1_t` 内部缓冲区的 view。
- 一旦对应 `bam1_t` 被释放或 move 后失效，view 不能继续使用。

---

## 3. Reader 抽象接口

`BamStreamReader` 抽象定义：

- `is_valid()`
- `bam_path()`
- `chromosome_count()`
- `chromosome_name(tid)`
- `stream(record_handler, progress_handler, progress_interval)`

工厂函数：

- `make_bam_reader(bam_path, decompression_threads)`

回调类型：

- `RecordHandler = std::function<void(BamRecordPtr&&)>`
- `ProgressHandler = std::function<bool(int64_t processed, int32_t current_tid)>`

---

## 4. 当前实现：`HtslibBamStreamReader`

## 4.1 构造和销毁

构造函数行为：

1. `hts_open(bam_path, "r")`
2. 若 `decompression_threads > 1`，调用 `hts_set_threads`
3. `sam_hdr_read`
4. 成功后 `valid_ = true`

失败路径：

- 打开 BAM 失败或读取 header 失败时，`valid_ = false`，并输出错误日志。

析构函数：

- 若 `header_ != nullptr`，调用 `bam_hdr_destroy(header_)`
- 若 `file_ != nullptr`，调用 `hts_close(file_)`

## 4.2 染色体信息接口

- `chromosome_count()` 返回 `header_->n_targets`（无 header 返回 0）。
- `chromosome_name(tid)` 在 `tid` 越界或 header 不可用时返回空字符串 `""`。

## 4.3 `stream()` 实际语义

前置检查：

- 若 reader 无效或 `record_handler` 为空，直接返回 `-1`。

主循环：

1. 每轮新建一个 `bam1_t`（`bam_init1`）。
2. 调用 `sam_read1(file_, header_, record.get())`。
3. 若返回值 `< 0`，结束循环（实现上不区分 EOF 与读错误）。
4. 过滤：
   - `BAM_FSECONDARY`：跳过
   - `BAM_FUNMAP`：跳过
5. 通过 `record_handler(std::move(record))` 交给上层。
6. `processed++`（计数的是“交给上层后的记录数”）。
7. 满足进度条件时触发 `progress_handler`：
   - 需要 `progress_handler != nullptr`
   - 需要 `progress_interval > 0`
   - 需要 `processed - last_progress >= progress_interval`
   - 若回调返回 `false`，提前停止读取

循环结束后：

- 如果提供了 `progress_handler`，会再调用一次最终回调。

返回值：

- 正常返回 `processed`。
- 记录对象分配失败时返回 `-1`。

---

## 5. `ReadView` 解码与边界行为

## 5.1 序列字符映射

解码表使用：

- `=ACMGRSVTWYHKDBN`

这是 htslib 4-bit 编码到字符的直接映射，保留了非 ACGT 符号。

## 5.2 `decode_subsequence(start, length)`

行为：

1. 若 `record_` 为空、`length<=0`、read 长度 `<=0`，返回空串。
2. `start < 0` 会被夹到 0。
3. `start >= read_len` 返回空串。
4. 终点 `end = min(start + length, read_len)`。
5. 仅解码 `[start, end)`。

## 5.3 `decode_sequence()`

- 解码整条 read 序列并返回 `std::string`。
- 若记录为空或 `l_qseq<=0` 返回空串。

## 5.4 tag 读取

- `has_tag` 基于 `bam_aux_get` 判断存在性。
- `get_int_tag` 使用 `bam_aux2i` 读取整数型值到 `int64_t`。
- `get_string_tag` 使用 `bam_aux2Z` 读取字符串；不存在或类型不匹配时返回 `false`。

---

## 6. 与 Pipeline 的实际集成

`Pipeline` 通过 `make_bam_reader` 获取 reader，并在两种模式下使用同一 `stream()`：

- `run_streaming()`：读取线程内直接处理每条记录。
- `run_parallel()`：读取线程产出 batch，worker 线程消费 batch。

集成要点：

1. BAM I/O 层已过滤 `unmapped/secondary`。
2. Gate1 层仍会再次检查这两类 flag（属于上层冗余保护，不影响正确性）。
3. progress 回调在 pipeline 中用于打印：
   - `processed`
   - `current_tid`

---

## 7. 性能与内存特征

1. 单遍流式读取，不做全量 reads 缓存。
2. 记录以 `BamRecordPtr` move 传递，避免对象复制。
3. 业务层多数逻辑可通过 `ReadView` 访问原始 CIGAR/tag，不必立即解码序列。
4. 仅在需要序列时调用 `decode_subsequence`/`decode_sequence` 分配字符串。

---

## 8. 当前边界与注意事项

1. `stream()` 对 `sam_read1 < 0` 不区分 EOF 与 I/O 错误。
2. 当前只实现顺序流式读取；未提供按区域随机访问接口。
3. `HtslibBamStreamReader` 不是并发读同一实例的线程安全接口（当前调用方式也未这样使用）。
4. 模块未在 I/O 层解析 SA/MD 语义，只提供 tag 原始访问能力。

---

## 9. 运行示例

串行模式（默认）：

```bash
./build/placer <input.bam> <ref.fa> [te.fa]
```

并行模式（生产者 + 单 worker）：

```bash
PLACER_PARALLEL=1 ./build/placer <input.bam> <ref.fa> [te.fa]
```
