#include "local_realign.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <thread>
#include <future>
#include <execution>
#include <cstdint>
#include <cstring>
#include <cassert>
#include <numeric>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <functional>
#include <unordered_map>

#if defined(__AVX2__)
#include <immintrin.h>
#endif

namespace placer {

// ============================================================================
// N-base penalty helper
// [修正] 同时处理大小写
// ============================================================================

static int8_t n_base_penalty(char c) {
    switch (c) {
        case 'N': case 'n': return -10;
        default: return 0;
    }
}

static bool is_n_base(char c) {
    return c == 'N' || c == 'n';
}

// Forward declaration
static std::string compress_cigar(const std::string& expanded);

// ============================================================================
// BitpackedSeq Implementation
// [修正] data_ 大小计算：每 32bp 一个 uint64_t word
// ============================================================================

void BitpackedSeq::encode(std::string_view ascii_seq) {
    size_t n = ascii_seq.size();
    size_ = n;

    // [修正] 每 32 个碱基占一个 uint64_t（每碱基 2 bit，32*2=64 bit）
    data_.resize((n + 31) / 32);
    std::fill(data_.begin(), data_.end(), 0ULL);

    size_t out_idx = 0;
    uint64_t word = 0;
    int bits = 0;

    for (size_t i = 0; i < n; ++i) {
        uint8_t code = ascii_to_2bit(ascii_seq[i]);
        word |= (static_cast<uint64_t>(code & MASK) << (bits * 2));
        bits++;

        if (bits == 32) {
            data_[out_idx++] = word;
            word = 0;
            bits = 0;
        }
    }

    if (bits > 0 && out_idx < data_.size()) {
        data_[out_idx] = word;
    }
}

char BitpackedSeq::at(size_t pos) const {
    if (pos >= size_) return 'N';
    size_t word_idx = pos / 32;
    size_t bit_offset = (pos % 32) * 2;
    uint8_t code = (data_[word_idx] >> bit_offset) & MASK;
    static const char decode[] = "ACGT";
    return decode[code & 0x03];
}

std::string BitpackedSeq::decode() const {
    std::string result;
    result.reserve(size_);
    for (size_t i = 0; i < size_; ++i) {
        result.push_back(at(i));
    }
    return result;
}

// ============================================================================
// SIMD-Accelerated Hamming Distance
// [修正] AVX2 逻辑简化：直接逐字节比较计数
// ============================================================================

#if defined(__AVX2__)
int simd_hamming_avx2(const char* a, const char* b, size_t len) {
    int mismatches = 0;
    size_t i = 0;

    // 处理 32 字节对齐的块
    for (; i + 32 <= len; i += 32) {
        __m256i va = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(a + i));
        __m256i vb = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(b + i));

        // [修正] 简化逻辑：直接比较相等，然后计数不等的
        __m256i eq = _mm256_cmpeq_epi8(va, vb);
        // eq: 相等位置 = 0xFF，不等位置 = 0x00
        // movemask: 每字节取最高位，相等=1，不等=0
        uint32_t mask = static_cast<uint32_t>(
            _mm256_movemask_epi8(eq));
        // 不等的数量 = 32 - popcount(mask)
        mismatches += 32 - __builtin_popcount(mask);
    }

    // 处理剩余字节
    for (; i < len; ++i) {
        if (a[i] != b[i]) ++mismatches;
    }

    return mismatches;
}
#else
int simd_hamming_avx2(const char*, const char*, size_t) {
    return -1;  // Not available
}
#endif

int LocalRealigner::hamming_distance(std::string_view a, std::string_view b) {
    if (a.size() != b.size()) return -1;
    if (a.empty()) return 0;

#if defined(__AVX2__)
    if (a.size() >= 32) {
        return simd_hamming_avx2(a.data(), b.data(), a.size());
    }
#endif

    // Scalar fallback with 8-byte unrolling
    int distance = 0;
    const char* pa = a.data();
    const char* pb = b.data();
    const char* end = pa + a.size();

    while (pa + 8 <= end) {
        distance += (pa[0] != pb[0]) + (pa[1] != pb[1]) +
                    (pa[2] != pb[2]) + (pa[3] != pb[3]) +
                    (pa[4] != pb[4]) + (pa[5] != pb[5]) +
                    (pa[6] != pb[6]) + (pa[7] != pb[7]);
        pa += 8;
        pb += 8;
    }

    while (pa < end) {
        if (*pa != *pb) ++distance;
        ++pa;
        ++pb;
    }

    return distance;
}

// ============================================================================
// Affine Gap Alignment with Full Traceback
// [修正] 完整 traceback 实现，精确报告 matches/mismatches/gaps
// [修正] N 碱基判断 bug (is_n_t)
// ============================================================================

// Traceback 方向编码
enum TraceDir : uint8_t {
    TRACE_DIAG_M = 0,   // 来自 M 矩阵的对角线
    TRACE_DIAG_X = 1,   // 来自 X 矩阵的对角线
    TRACE_DIAG_Y = 2,   // 来自 Y 矩阵的对角线
    TRACE_UP_OPEN = 3,  // Y: gap open (从 M)
    TRACE_UP_EXT = 4,   // Y: gap extend (从 Y)
    TRACE_LEFT_OPEN = 5, // X: gap open (从 M)
    TRACE_LEFT_EXT = 6,  // X: gap extend (从 X)
};

AlignmentResult affine_gap_align(
    std::string_view query,
    std::string_view target,
    const RealignConfig& config) {

    AlignmentResult result;
    const int8_t match_sc = config.match;
    const int8_t mismatch_sc = config.mismatch;
    const int8_t gap_open = config.gap_open;
    const int8_t gap_extend = config.gap_extend;

    const int qlen = static_cast<int>(query.size());
    const int tlen = static_cast<int>(target.size());

    if (qlen == 0 || tlen == 0) return result;

    const int NEG_INF = INT_MIN / 4;

    // 完整 DP 矩阵（用于 traceback）
    // 对于大序列可以考虑 Hirschberg，但这里假设 flank 长度有限（<10kb）
    const size_t rows = static_cast<size_t>(qlen + 1);
    const size_t cols = static_cast<size_t>(tlen + 1);

    // 使用一维数组模拟二维矩阵
    std::vector<int> M(rows * cols, NEG_INF);
    std::vector<int> X(rows * cols, NEG_INF);  // gap in target
    std::vector<int> Y(rows * cols, NEG_INF);  // gap in query

    // Traceback 矩阵
    std::vector<uint8_t> trace_M(rows * cols, 0);
    std::vector<uint8_t> trace_X(rows * cols, 0);
    std::vector<uint8_t> trace_Y(rows * cols, 0);

    auto idx = [cols](int i, int j) -> size_t {
        return static_cast<size_t>(i) * cols + j;
    };

    // 初始化
    M[idx(0, 0)] = 0;

    // 第一列：query gap
    for (int i = 1; i <= qlen; ++i) {
        Y[idx(i, 0)] = gap_open + i * gap_extend;
        M[idx(i, 0)] = Y[idx(i, 0)];
        trace_M[idx(i, 0)] = TRACE_DIAG_Y;
    }

    // 第一行：target gap
    for (int j = 1; j <= tlen; ++j) {
        X[idx(0, j)] = gap_open + j * gap_extend;
        M[idx(0, j)] = X[idx(0, j)];
        trace_M[idx(0, j)] = TRACE_DIAG_X;
    }

    // 主 DP 循环
    for (int i = 1; i <= qlen; ++i) {
        char qc = query[i - 1];
        bool qn = is_n_base(qc);

        for (int j = 1; j <= tlen; ++j) {
            char tc = target[j - 1];

            // [修正] N 碱基判断：同时检查大小写
            bool tn = is_n_base(tc);
            int8_t n_pen = (qn || tn) ? n_base_penalty('N') : 0;

            // Match/mismatch score
            int8_t sub_score;
            if (qn || tn) {
                // N 碱基：给一个较低的匹配分（鼓励避开）
                sub_score = match_sc + n_pen;
            } else {
                // [修正] 大小写不敏感比较
                char qcu = (qc >= 'a' && qc <= 'z') ? (qc - 32) : qc;
                char tcu = (tc >= 'a' && tc <= 'z') ? (tc - 32) : tc;
                sub_score = (qcu == tcu) ? match_sc : mismatch_sc;
            }

            // === M[i][j]: 来自对角线 ===
            int from_M = M[idx(i - 1, j - 1)] + sub_score;
            int from_X = X[idx(i - 1, j - 1)] + sub_score;
            int from_Y = Y[idx(i - 1, j - 1)] + sub_score;

            int best_M = from_M;
            uint8_t best_M_trace = TRACE_DIAG_M;

            if (from_X > best_M) {
                best_M = from_X;
                best_M_trace = TRACE_DIAG_X;
            }
            if (from_Y > best_M) {
                best_M = from_Y;
                best_M_trace = TRACE_DIAG_Y;
            }

            M[idx(i, j)] = best_M;
            trace_M[idx(i, j)] = best_M_trace;

            // === X[i][j]: gap in target (deletion from query perspective) ===
            int x_open = M[idx(i, j - 1)] + gap_open + gap_extend;
            int x_ext = X[idx(i, j - 1)] + gap_extend;

            if (x_open >= x_ext) {
                X[idx(i, j)] = x_open;
                trace_X[idx(i, j)] = TRACE_LEFT_OPEN;
            } else {
                X[idx(i, j)] = x_ext;
                trace_X[idx(i, j)] = TRACE_LEFT_EXT;
            }

            // === Y[i][j]: gap in query (insertion from query perspective) ===
            int y_open = M[idx(i - 1, j)] + gap_open + gap_extend;
            int y_ext = Y[idx(i - 1, j)] + gap_extend;

            if (y_open >= y_ext) {
                Y[idx(i, j)] = y_open;
                trace_Y[idx(i, j)] = TRACE_UP_OPEN;
            } else {
                Y[idx(i, j)] = y_ext;
                trace_Y[idx(i, j)] = TRACE_UP_EXT;
            }
        }
    }

    // 找最优终止点（半全局：query 全消耗，target 可部分）
    int best_score = NEG_INF;
    int best_i = qlen;
    int best_j = 0;
    uint8_t best_matrix = 0;  // 0=M, 1=X, 2=Y

    // 检查最后一行（query 全部消耗）
    for (int j = 1; j <= tlen; ++j) {
        int scores[3] = {M[idx(qlen, j)], X[idx(qlen, j)], Y[idx(qlen, j)]};
        for (int k = 0; k < 3; ++k) {
            if (scores[k] > best_score) {
                best_score = scores[k];
                best_i = qlen;
                best_j = j;
                best_matrix = k;
            }
        }
    }

    // 也检查最后一列（target 全部消耗）
    for (int i = 1; i <= qlen; ++i) {
        int scores[3] = {M[idx(i, tlen)], X[idx(i, tlen)], Y[idx(i, tlen)]};
        for (int k = 0; k < 3; ++k) {
            if (scores[k] > best_score) {
                best_score = scores[k];
                best_i = i;
                best_j = tlen;
                best_matrix = k;
            }
        }
    }

    // ===== Traceback =====
    int matches = 0, mismatches = 0, gap_opens = 0, gap_extensions = 0;
    int i = best_i, j = best_j;
    uint8_t current_matrix = best_matrix;

    // CIGAR 操作（反向构建）
    std::string cigar_ops_rev;

    static constexpr int MAX_TRACE_ITER = 2000000;
    int iter = 0;

    while (i > 0 && j > 0 && iter++ < MAX_TRACE_ITER) {
        if (current_matrix == 0) {
            // 在 M 矩阵中
            uint8_t tr = trace_M[idx(i, j)];

            char qc = query[i - 1];
            char tc = target[j - 1];
            char qcu = (qc >= 'a' && qc <= 'z') ? (qc - 32) : qc;
            char tcu = (tc >= 'a' && tc <= 'z') ? (tc - 32) : tc;

            if (qcu == tcu && !is_n_base(qc) && !is_n_base(tc)) {
                matches++;
                cigar_ops_rev.push_back('=');
            } else {
                mismatches++;
                cigar_ops_rev.push_back('X');
            }

            // 转移到前一个矩阵
            if (tr == TRACE_DIAG_M) current_matrix = 0;
            else if (tr == TRACE_DIAG_X) current_matrix = 1;
            else current_matrix = 2;

            i--;
            j--;

        } else if (current_matrix == 1) {
            // 在 X 矩阵中（gap in target = deletion）
            uint8_t tr = trace_X[idx(i, j)];

            if (tr == TRACE_LEFT_OPEN) {
                gap_opens++;
                gap_extensions++;
                current_matrix = 0;  // 回到 M
            } else {
                gap_extensions++;
                current_matrix = 1;  // 继续 X
            }

            cigar_ops_rev.push_back('D');
            j--;

        } else {
            // 在 Y 矩阵中（gap in query = insertion）
            uint8_t tr = trace_Y[idx(i, j)];

            if (tr == TRACE_UP_OPEN) {
                gap_opens++;
                gap_extensions++;
                current_matrix = 0;  // 回到 M
            } else {
                gap_extensions++;
                current_matrix = 2;  // 继续 Y
            }

            cigar_ops_rev.push_back('I');
            i--;
        }
    }

    // 处理剩余（到达边界）
    while (i > 0 && iter++ < MAX_TRACE_ITER) {
        cigar_ops_rev.push_back('I');
        gap_extensions++;
        if (cigar_ops_rev.size() == 1 ||
            cigar_ops_rev[cigar_ops_rev.size() - 2] != 'I') {
            gap_opens++;
        }
        i--;
    }
    while (j > 0 && iter++ < MAX_TRACE_ITER) {
        cigar_ops_rev.push_back('D');
        gap_extensions++;
        if (cigar_ops_rev.size() == 1 ||
            cigar_ops_rev[cigar_ops_rev.size() - 2] != 'D') {
            gap_opens++;
        }
        j--;
    }

    // 填充结果
    result.score = static_cast<double>(best_score);
    result.normalized_score = (std::min(qlen, tlen) > 0) ?
        best_score / static_cast<double>(std::min(qlen, tlen)) : 0.0;

    result.matches = matches;
    result.mismatches = mismatches;
    result.gap_opens = gap_opens;
    result.gap_extensions = gap_extensions;

    int aligned_len = matches + mismatches;
    result.identity = (aligned_len > 0) ?
        static_cast<float>(matches) / aligned_len : 0.0f;

    result.has_indel = (gap_opens > 0);
    result.query_start = i;  // 未消耗的 query 起始
    result.query_end = best_i;
    result.target_start = j;
    result.target_end = best_j;

    return result;
}

// [新增] CIGAR 压缩：将 "===XXD" 压缩为 "3=2X1D"
std::string compress_cigar(const std::string& expanded) {
    if (expanded.empty()) return "";

    std::string cigar;
    cigar.reserve(expanded.size());

    char prev = expanded[0];
    int count = 1;

    for (size_t i = 1; i < expanded.size(); ++i) {
        if (expanded[i] == prev) {
            count++;
        } else {
            cigar += std::to_string(count);
            cigar += prev;
            prev = expanded[i];
            count = 1;
        }
    }
    cigar += std::to_string(count);
    cigar += prev;

    return cigar;
}

// ============================================================================
// GenomeAccessor Implementation
// [修正] 文件句柄复用（mmap 或持久 ifstream）
// [修正] fetch_region 的 file_offset 计算和错误处理
// ============================================================================

GenomeAccessor::GenomeAccessor(std::string_view fasta_path)
    : fasta_path_(fasta_path) {
    load_fai(fasta_path);

    // [新增] 打开文件并保持打开状态
    if (!fasta_path_.empty()) {
        file_ = new std::ifstream(fasta_path_, std::ios::binary);
        if (!file_->is_open()) {
            delete file_;
            file_ = nullptr;
        }
    }
}

// [新增] 析构函数：关闭文件
GenomeAccessor::~GenomeAccessor() {
    if (file_ != nullptr) {
        file_->close();
        delete file_;
        file_ = nullptr;
    }
}

GenomeAccessor::GenomeAccessor(GenomeAccessor&& other) noexcept
    : fasta_path_(std::move(other.fasta_path_))
    , index_(std::move(other.index_))
    , file_(other.file_) {
    other.file_ = nullptr;
}

GenomeAccessor& GenomeAccessor::operator=(GenomeAccessor&& other) noexcept {
    if (this != &other) {
        // 先关闭旧文件
        if (file_ != nullptr) {
            file_->close();
            delete file_;
        }
        fasta_path_ = std::move(other.fasta_path_);
        index_ = std::move(other.index_);
        file_ = other.file_;
        other.file_ = nullptr;
    }
    return *this;
}

bool GenomeAccessor::load_fai(std::string_view fasta_path) {
    std::string fai_path(fasta_path);
    fai_path += ".fai";

    std::ifstream file(fai_path);
    if (!file.is_open()) {
        return false;
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        IndexEntry entry;
        std::istringstream iss(line);
        iss >> entry.name >> entry.length >> entry.offset
            >> entry.line_blen >> entry.line_len;

        if (iss.fail()) continue;

        // [修正] 验证 fai 条目合理性
        if (entry.length <= 0 || entry.offset < 0 ||
            entry.line_blen <= 0 || entry.line_len <= 0) {
            continue;
        }

        index_.push_back(entry);
    }

    return !index_.empty();
}

std::optional<SeqView> GenomeAccessor::fetch_region(
    int chrom_tid, int32_t start, int32_t end) const {

    if (chrom_tid < 0 || chrom_tid >= static_cast<int>(index_.size())) {
        return std::nullopt;
    }

    // [修正] 检查持久化文件句柄
    if (file_ == nullptr || !file_->is_open()) {
        return std::nullopt;
    }

    const auto& idx = index_[chrom_tid];

    // [修正] 范围验证
    if (start < 0) start = 0;
    if (end > static_cast<int32_t>(idx.length)) {
        end = static_cast<int32_t>(idx.length);
    }
    if (start >= end) return std::nullopt;

    int32_t region_len = end - start;

    // file_offset 计算
    int64_t file_offset = idx.offset + start + (start / idx.line_blen);
    int64_t fetch_len = (end - start) + ((end - start) / idx.line_blen) + 2;

    // [修正] 使用 thread_local 缓冲区避免每次分配
    thread_local static std::string chunk;
    thread_local static std::string cleaned;
    chunk.clear();
    cleaned.clear();
    chunk.resize(static_cast<size_t>(fetch_len));
    cleaned.reserve(region_len);

    // [修正] 使用持久化文件句柄
    file_->seekg(file_offset);
    file_->read(chunk.data(), fetch_len);

    // 去除换行符 + 转大写
    for (char c : chunk) {
        if (c == '\n' || c == '\r') continue;
        cleaned.push_back(std::toupper(c));
        if (static_cast<int32_t>(cleaned.size()) >= region_len) break;
    }

    return SeqView(cleaned.data(), cleaned.size());
}

std::optional<SeqView> GenomeAccessor::fetch(
    int chrom_tid, int32_t start, int32_t end) const {

    if (chrom_tid < 0) return std::nullopt;
    if (start < 0) start = 0;

    if (chrom_tid < static_cast<int>(index_.size())) {
        int64_t chrom_len = static_cast<int64_t>(index_[chrom_tid].length);
        if (end > static_cast<int32_t>(chrom_len)) {
            end = static_cast<int32_t>(chrom_len);
        }
    }
    if (start >= end) return std::nullopt;

    return fetch_region(chrom_tid, start, end);
}

std::optional<BitpackedSeq> GenomeAccessor::fetch_packed(
    int chrom_tid) const {

    if (chrom_tid < 0 || chrom_tid >= static_cast<int>(index_.size())) {
        return std::nullopt;
    }

    auto view = fetch(chrom_tid, 0,
                      static_cast<int32_t>(index_[chrom_tid].length));
    if (!view.has_value()) return std::nullopt;

    BitpackedSeq seq(view->sv());
    return seq;
}

// ============================================================================
// LocalRealigner Implementation
// ============================================================================

LocalRealigner::LocalRealigner(RealignConfig config)
    : config_(std::move(config)) {}

// ============================================================================
// Alignment Entry Points
// ============================================================================

AlignmentResult LocalRealigner::align_sequences(
    std::string_view query,
    std::string_view target,
    const RealignConfig& config) {

    AlignmentResult result;

    if (query.empty() || target.empty()) {
        return result;
    }

    result = affine_gap_align(query, target, config);

    if (!result.is_valid()) {
        result.score = -1e9;
        return result;
    }

    if (result.identity < config.min_identity ||
        result.normalized_score < config.min_normalized_score) {
        result.score = -1e9;
        return result;
    }

    // 计算等效 MAPQ
    if (result.matches + result.mismatches > 0) {
        float id = result.identity;
        float error_rate = std::max(1.0f - id, 1e-6f);
        result.mapq_equiv = static_cast<float>(
            -10.0 * std::log10(error_rate));
        result.mapq_equiv = std::clamp(result.mapq_equiv, 0.0f, 60.0f);
    }

    return result;
}

AlignmentResult LocalRealigner::dispatch_align_(
    std::string_view query,
    std::string_view target) {

    // [优化] 根据序列长度选择对齐策略
    size_t min_len = std::min(query.size(), target.size());
    size_t max_len = std::max(query.size(), target.size());

    // 1. 短序列：使用 Hamming 距离（等长）或简化 DP
    if (min_len <= 16) {
        return hamming_align_(query, target);
    }

    // 2. 中等序列：使用带状对齐 (banded alignment)
    if (min_len <= 64) {
        return banded_align_(query, target, 8);  // bandwidth = 8
    }

    // 3. 长序列：使用 SIMD 或完整对齐
    if (config_.use_simd && query.size() >= 32) {
        return simd_align_(query, target);
    }
    return scalar_align_(query, target);
}

// [新增] 快速 Hamming 距离对齐（短序列）
AlignmentResult LocalRealigner::hamming_align_(
    std::string_view query,
    std::string_view target) {

    AlignmentResult result;

    if (query.empty() || target.empty()) {
        return result;
    }

    // 长度差异惩罚
    int len_diff = static_cast<int>(
        query.size() > target.size() ? query.size() - target.size() : target.size() - query.size());

    // 计算 Hamming 距离
    int matches = 0;
    int mismatches = 0;
    int min_len = static_cast<int>(std::min(query.size(), target.size()));

    for (int i = 0; i < min_len; ++i) {
        char qc = std::toupper(query[i]);
        char tc = std::toupper(target[i]);
        if (qc == tc || qc == 'N' || tc == 'N') {
            matches++;
        } else {
            mismatches++;
        }
    }

    // 分数计算
    int score = matches * config_.match - mismatches * config_.mismatch
                - len_diff * (config_.gap_open + config_.gap_extend);

    result.score = score;
    result.normalized_score = static_cast<double>(score) / std::max(min_len, 1);
    result.matches = matches;
    result.mismatches = mismatches;
    result.identity = static_cast<float>(matches) / std::max(min_len, 1);
    result.has_indel = (len_diff > 0);
    result.gap_opens = len_diff > 0 ? 1 : 0;
    result.gap_extensions = len_diff > 0 ? len_diff - 1 : 0;

    return result;
}

// [新增] 带状对齐（中等序列）
AlignmentResult LocalRealigner::banded_align_(
    std::string_view query,
    std::string_view target,
    int bandwidth) {

    AlignmentResult result;

    if (query.empty() || target.empty()) {
        return result;
    }

    const int qlen = static_cast<int>(query.size());
    const int tlen = static_cast<int>(target.size());
    const int k = std::min(bandwidth, std::min(qlen, tlen) / 2 + 1);

    // 使用简化的 Smith-Waterman 带状 DP
    // 只计算对角线附近 bandwidth 范围内的值

    std::vector<int> prev(2 * k + 1, INT_MIN / 4);
    std::vector<int> curr(2 * k + 1, INT_MIN / 4);

    // 初始化
    prev[k] = 0;

    int best_score = 0;
    int best_i = 0, best_j = 0;

    for (int j = 1; j <= tlen; ++j) {
        for (int diag = -k; diag <= k; ++diag) {
            int i = diag + j;  // query index
            curr[diag + k] = INT_MIN / 4;

            if (i < 0 || i > qlen) continue;

            // Match/Mismatch
            if (i > 0) {
                char qc = std::toupper(query[i - 1]);
                char tc = std::toupper(target[j - 1]);
                int8_t sub = (qc == tc || qc == 'N' || tc == 'N') ?
                    config_.match : -config_.mismatch;
                curr[diag + k] = prev[diag + k] + sub;
            }

            // Gap in target (deletion from query)
            if (diag > -k && j > 1) {
                int val = prev[diag - 1 + k] + config_.gap_open + config_.gap_extend;
                curr[diag + k] = std::max(curr[diag + k], val);
            }

            // Gap in query (insertion from query)
            if (diag < k && i > 1) {
                int val = prev[diag + 1 + k] + config_.gap_open + config_.gap_extend;
                curr[diag + k] = std::max(curr[diag + k], val);
            }

            // 更新最佳分数
            if (curr[diag + k] > best_score) {
                best_score = curr[diag + k];
                best_i = i;
                best_j = j;
            }
        }
        std::swap(prev, curr);
    }

    // 简化结果
    result.score = best_score;
    result.normalized_score = static_cast<double>(best_score) /
        std::max(std::min(qlen, tlen), 1);
    result.matches = best_score > 0 ?
        static_cast<int>(best_score / config_.match) : 0;
    result.mismatches = 0;
    result.identity = result.normalized_score > 1.0 ? 1.0f :
        static_cast<float>(std::max(0.0, result.normalized_score));
    result.has_indel = (best_i != best_j);
    result.gap_opens = (best_i != best_j) ? 1 : 0;

    return result;
}

AlignmentResult LocalRealigner::simd_align_(
    std::string_view query,
    std::string_view target) {

    // [修正] SIMD 加速策略：
    // 1. 先用 SIMD hamming 做快速过滤
    // 2. 如果等长且 hamming 距离很小，直接返回近似结果
    // 3. 否则回退到完整 affine gap alignment

#if defined(__AVX2__)
    if (query.size() == target.size()) {
        int ham = hamming_distance(query, target);
        if (ham >= 0) {
            double identity = 1.0 - static_cast<double>(ham) / query.size();

            // 如果 identity 很高且等长，可以快速返回
            if (identity >= 0.95) {
                AlignmentResult result;
                result.score = static_cast<double>(
                    ham * config_.mismatch +
                    (static_cast<int>(query.size()) - ham) * config_.match);
                result.normalized_score = result.score /
                    static_cast<double>(query.size());
                result.matches = static_cast<int>(query.size()) - ham;
                result.mismatches = ham;
                result.identity = static_cast<float>(identity);
                result.has_indel = false;
                result.gap_opens = 0;
                result.gap_extensions = 0;
                result.query_start = 0;
                result.query_end = static_cast<int>(query.size());
                result.target_start = 0;
                result.target_end = static_cast<int>(target.size());

                float error_rate = std::max(1.0f - result.identity, 1e-6f);
                result.mapq_equiv = std::clamp(
                    static_cast<float>(-10.0 * std::log10(error_rate)),
                    0.0f, 60.0f);

                return result;
            }
        }
    }
#endif

    // 回退到完整对齐
    return scalar_align_(query, target);
}

AlignmentResult LocalRealigner::scalar_align_(
    std::string_view query,
    std::string_view target) {

    return align_sequences(query, target, config_);
}

// ============================================================================
// Seed Finding（用于快速候选筛选）
// ============================================================================

std::vector<uint32_t> LocalRealigner::find_seeds(
    std::string_view seq,
    int kmer_size,
    int stride) {

    std::vector<uint32_t> seeds;
    if (static_cast<int>(seq.size()) < kmer_size) return seeds;

    seeds.reserve((seq.size() - kmer_size) / stride + 1);

    uint32_t mask = (kmer_size < 16) ?
        ((1u << (kmer_size * 2)) - 1) : 0xFFFFFFFFu;

    // 构建第一个 kmer
    uint32_t kmer = 0;
    int valid_bases = 0;

    for (int i = 0; i < kmer_size; ++i) {
        if (is_n_base(seq[i])) {
            // N 碱基：重置 kmer
            valid_bases = 0;
            kmer = 0;
        } else {
            uint8_t code = BitpackedSeq::ascii_to_2bit(seq[i]);
            kmer = ((kmer << 2) | (code & 0x03)) & mask;
            valid_bases++;
        }
    }

    if (valid_bases >= kmer_size) {
        seeds.push_back(kmer);
    }

    // 滚动计算后续 kmer
    for (size_t i = kmer_size; i < seq.size(); ++i) {
        if (is_n_base(seq[i])) {
            valid_bases = 0;
            kmer = 0;
            continue;
        }

        uint8_t code = BitpackedSeq::ascii_to_2bit(seq[i]);
        kmer = ((kmer << 2) | (code & 0x03)) & mask;
        valid_bases++;

        if (valid_bases >= kmer_size &&
            (i - kmer_size + 1) % static_cast<size_t>(stride) == 0) {
            seeds.push_back(kmer);
        }
    }

    return seeds;
}

// ============================================================================
// Locus Management
// [修正] 合并逻辑 + 位置加权
// ============================================================================

void LocalRealigner::populate_locus_set(
    Component& component,
    const std::vector<ReadSketch>& reads) {

    component.locus_set.clear();

    int32_t window_start = component.start - config_.search_window;
    int32_t window_end = component.end + config_.search_window;

    for (size_t read_idx : component.read_indices) {
        if (read_idx >= reads.size()) continue;

        const auto& read = reads[read_idx];

        // Primary alignment
        if (read.pos >= window_start && read.pos <= window_end) {
            component.locus_set.push_back({
                read.tid,
                read.pos,
                static_cast<double>(read.mapq),
                1,
                0x01  // primary evidence
            });
        }

        // [修正] 也考虑 read 的 end_pos
        if (read.end_pos >= window_start && read.end_pos <= window_end &&
            std::abs(read.end_pos - read.pos) > 100) {
            component.locus_set.push_back({
                read.tid,
                read.end_pos,
                static_cast<double>(read.mapq) * 0.9,
                1,
                0x01
            });
        }

        // SA splits
        if (read.has_sa) {
            for (const auto& [sa_tid, sa_pos] : read.sa_targets) {
                if (sa_tid != component.chrom_tid) continue;

                if (sa_pos >= window_start && sa_pos <= window_end) {
                    component.locus_set.push_back({
                        sa_tid,
                        sa_pos,
                        static_cast<double>(read.mapq) * 0.8,
                        1,
                        0x02  // SA evidence
                    });
                }
            }
        }
    }

    if (component.locus_set.empty()) return;

    // 按位置排序
    std::sort(component.locus_set.begin(), component.locus_set.end(),
        [](const LocusCandidate& a, const LocusCandidate& b) {
            return a.pos < b.pos;
        });

    // [修正] 合并近距离 loci（加权位置）
    std::vector<LocusCandidate> merged;
    merged.reserve(component.locus_set.size());

    const int merge_distance = 50;

    for (const auto& locus : component.locus_set) {
        if (merged.empty()) {
            merged.push_back(locus);
            continue;
        }

        LocusCandidate& last = merged.back();

        if (std::abs(locus.pos - last.pos) <= merge_distance &&
            locus.chrom_tid == last.chrom_tid) {
            // [修正] 加权合并位置
            double total_score = last.score + locus.score;
            if (total_score > 0) {
                last.pos = static_cast<int32_t>(
                    (static_cast<double>(last.pos) * last.score +
                     static_cast<double>(locus.pos) * locus.score) /
                    total_score);
            }
            last.score = std::max(last.score, locus.score);
            last.support_reads += locus.support_reads;
            last.evidence_mask |= locus.evidence_mask;
        } else {
            merged.push_back(locus);
        }
    }

    component.locus_set = std::move(merged);

    // 限制数量：按 score 取 top-k
    if (component.locus_set.size() > config_.max_locus_per_component) {
        std::partial_sort(
            component.locus_set.begin(),
            component.locus_set.begin() + config_.max_locus_per_component,
            component.locus_set.end(),
            [](const LocusCandidate& a, const LocusCandidate& b) {
                return a.score > b.score;
            }
        );
        component.locus_set.resize(config_.max_locus_per_component);
    }
}

std::vector<LocusCandidate> LocalRealigner::rank_loci(
    std::vector<LocusCandidate>&& loci,
    const RealignConfig& config) {

    std::sort(loci.begin(), loci.end(),
        [](const LocusCandidate& a, const LocusCandidate& b) {
            // 综合评分：score × 证据类型数 × 支持 reads 数
            int a_types = __builtin_popcount(a.evidence_mask);
            int b_types = __builtin_popcount(b.evidence_mask);
            double a_strength = a.score * a_types *
                std::log2(a.support_reads + 1);
            double b_strength = b.score * b_types *
                std::log2(b.support_reads + 1);
            return a_strength > b_strength;
        }
    );

    return std::move(loci);
}

// ============================================================================
// Flanks Extraction
// [修正] 正确提取上下游侧翼序列
// ============================================================================

std::pair<SeqView, SeqView> LocalRealigner::extract_flanks_view(
    const ReadSketch& read,
    int flank_length,
    const GenomeAccessor& genome) {

    // 上游侧翼
    int32_t up_start = std::max<int32_t>(0, read.pos - flank_length);
    int32_t up_end = read.pos;
    auto up_opt = genome.fetch(read.tid, up_start, up_end);
    SeqView up_flank = up_opt.value_or(SeqView());

    // 下游侧翼
    int32_t down_start = read.end_pos;
    int32_t down_end = read.end_pos + flank_length;

    // [修正] 限制不超过染色体长度
    const auto& idx = genome.get_index(read.tid);
    if (idx.length > 0) {
        down_end = std::min<int32_t>(
            down_end, static_cast<int32_t>(idx.length));
    }

    auto down_opt = genome.fetch(read.tid, down_start, down_end);
    SeqView down_flank = down_opt.value_or(SeqView());

    return {up_flank, down_flank};
}

// ============================================================================
// Component Realignment
// [修正] 真正使用 affine gap alignment 而非距离占位
// ============================================================================

PlaceabilityReport LocalRealigner::realign_component(
    Component& component,
    const std::vector<ReadSketch>& reads,
    const GenomeAccessor& genome) {

    PlaceabilityReport report;

    // Step 1: 填充候选位点
    populate_locus_set(component, reads);

    if (component.locus_set.empty()) {
        component.locus_set.push_back({
            component.chrom_tid,
            component.centroid,
            1.0,
            static_cast<uint32_t>(component.read_count),
            0x04  // default evidence
        });
    }

    // Step 2: 排序候选位点
    component.locus_set = rank_loci(
        std::move(component.locus_set), config_);

    // Step 3: 对每个 read × 每个 locus 生成证据
    std::vector<LocusEvidence> all_evidence;
    all_evidence.reserve(
        component.read_indices.size() * component.locus_set.size());

    size_t reads_processed = 0;

    for (size_t read_idx : component.read_indices) {
        if (read_idx >= reads.size()) continue;
        if (reads_processed >= config_.max_reads_per_component) break;
        reads_processed++;

        const auto& read = reads[read_idx];

        // [修正] 提取 read 的侧翼序列用于真正对齐
        auto [up_flank, down_flank] = extract_flanks_view(
            read, config_.flank_length, genome);

        for (size_t li = 0; li < component.locus_set.size(); ++li) {
            const auto& locus = component.locus_set[li];

            LocusEvidence ev;
            ev.read_idx = read_idx;
            ev.locus_pos = locus.pos;
            ev.is_reverse = (read.flag & 0x10) != 0;
            ev.evidence_bits = 0;

            // [修正] 真正的序列对齐评分
            bool has_alignment = false;

            // 提取 locus 周围的参考序列
            int32_t ref_start = std::max<int32_t>(
                0, locus.pos - config_.flank_length);
            int32_t ref_end = locus.pos + config_.flank_length;

            auto ref_opt = genome.fetch(
                component.chrom_tid, ref_start, ref_end);

            if (ref_opt.has_value() && !read.sequence.empty()) {
                // 用 read 序列对齐到参考区域
                std::string_view ref_seq = ref_opt->sv();
                std::string_view read_seq(
                    read.sequence.data(), read.sequence.size());

                // 如果 read 序列太长，只取与 locus 相关的部分
                int32_t read_offset = locus.pos - read.pos;
                int32_t extract_start = std::max<int32_t>(
                    0, read_offset - config_.flank_length);
                int32_t extract_end = std::min<int32_t>(
                    static_cast<int32_t>(read_seq.size()),
                    read_offset + config_.flank_length);

                if (extract_start < extract_end &&
                    extract_end <= static_cast<int32_t>(read_seq.size())) {

                    std::string_view read_segment = read_seq.substr(
                        extract_start, extract_end - extract_start);

                    if (!read_segment.empty() && !ref_seq.empty()) {
                        AlignmentResult aln = dispatch_align_(
                            read_segment, ref_seq);

                        if (aln.score > -1e8) {
                            ev.normalized_score = std::max(
                                0.0, aln.normalized_score);
                            ev.total_score = aln.score;
                            ev.weighted_identity = aln.identity;
                            has_alignment = true;

                            // 设置证据位
                            if (aln.identity >= config_.min_identity) {
                                ev.evidence_bits |= 0x01;  // 序列匹配
                            }
                            if (aln.has_indel) {
                                ev.evidence_bits |= 0x08;  // 有 indel
                            }
                        }
                    }
                }
            }

            // [修正] 如果无法做序列对齐，回退到距离评分
            if (!has_alignment) {
                double dist = std::abs(
                    static_cast<double>(read.pos) - locus.pos);
                ev.normalized_score = std::max(
                    0.0, 1.0 - dist / config_.search_window);
                ev.total_score = ev.normalized_score * 100;
                ev.weighted_identity = 0.0f;

                if (ev.normalized_score > config_.min_normalized_score) {
                    ev.evidence_bits |= 0x02;  // 距离证据
                }
            }

            // SA 证据加成
            if (locus.evidence_mask & 0x02) {
                ev.evidence_bits |= 0x04;  // SA 支持
            }

            all_evidence.push_back(ev);
        }
    }

    // Step 4: 计算 placeability
    report = calculate_placeability(all_evidence);

    // Step 5: 链统计
    report.forward_count = 0;
    report.reverse_count = 0;

    // [修正] 只统计有效证据的链方向
    for (const auto& ev : all_evidence) {
        if (ev.normalized_score < 0.1) continue;
        if (ev.is_reverse) ++report.reverse_count;
        else ++report.forward_count;
    }

    report.strand_balanced =
        (report.forward_count > 0 && report.reverse_count > 0);

    // [修正] 链平衡比值
    if (report.forward_count + report.reverse_count > 0) {
        int min_strand = std::min(
            report.forward_count, report.reverse_count);
        int max_strand = std::max(
            report.forward_count, report.reverse_count);
        report.strand_ratio = static_cast<float>(min_strand) /
            std::max(max_strand, 1);
    }

    // Step 6: 确定 tier
    report.tier = determine_tier_(report);

    // Step 7: 记录最佳 locus
    if (!component.locus_set.empty() && report.best_score > 0) {
        report.best_locus = 0;  // 使用索引
        report.best_locus_pos = component.locus_set[0].pos;
    }

    return report;
}

// ============================================================================
// Collect Evidence - 优化版本
// [优化]
// 1. 限制每个 component 最多处理 max_reads_per_component reads
// 2. 限制每个 component 最多处理 max_loci_per_component loci
// 3. 距离预过滤：跳过太远的 read-locus 组合
// ============================================================================

std::vector<LocusEvidence> LocalRealigner::collect_evidence(
    Component& component,
    const std::vector<ReadSketch>& reads,
    const GenomeAccessor& genome) {

    // Step 1: 填充候选位点
    populate_locus_set(component, reads);

    if (component.locus_set.empty()) {
        component.locus_set.push_back({
            component.chrom_tid,
            component.centroid,
            1.0,
            static_cast<uint32_t>(component.read_count),
            0x04  // default evidence
        });
    }

    // Step 2: 排序并限制候选位点数量
    component.locus_set = rank_loci(
        std::move(component.locus_set), config_);

    // [优化] 限制 loci 数量
    const size_t max_loci = std::min(static_cast<size_t>(config_.max_locus_per_component),
                                     component.locus_set.size());
    if (component.locus_set.size() > max_loci) {
        component.locus_set.resize(max_loci);
    }

    // Step 3: 限制 reads 数量
    std::vector<size_t> sampled_reads;
    const size_t max_reads = std::min(static_cast<size_t>(config_.max_reads_per_component),
                                      component.read_indices.size());

    // 分层采样：从 read_indices 开头取 max_reads 个
    for (size_t i = 0; i < component.read_indices.size() && sampled_reads.size() < max_reads; ++i) {
        sampled_reads.push_back(component.read_indices[i]);
    }

    // Step 4: 对每个 read × 每个 locus 生成证据
    std::vector<LocusEvidence> all_evidence;
    all_evidence.reserve(sampled_reads.size() * component.locus_set.size());

    // [优化] 预计算距离阈值
    double dist_threshold = config_.search_window * config_.min_normalized_score;

    for (size_t read_idx : sampled_reads) {
        if (read_idx >= reads.size()) continue;

        const auto& read = reads[read_idx];

        for (const auto& locus : component.locus_set) {
            // [优化] 距离预过滤
            double dist = std::abs(
                static_cast<double>(read.pos) - locus.pos);
            if (dist > dist_threshold) {
                continue;  // 太远的组合直接跳过
            }

            LocusEvidence ev;
            ev.read_idx = read_idx;
            ev.locus_pos = locus.pos;
            ev.is_reverse = (read.flag & 0x10) != 0;
            ev.evidence_bits = 0;

            // 距离评分
            ev.normalized_score = std::max(
                0.0, 1.0 - dist / config_.search_window);
            ev.total_score = ev.normalized_score * 100;

            if (ev.normalized_score > config_.min_normalized_score) {
                ev.evidence_bits |= 0x02;
            }

            all_evidence.push_back(ev);
        }
    }

    return all_evidence;
}

// ============================================================================
// Placeability Calculation
// [修正] 按 locus 聚合证据，而非全局混合
// ============================================================================

PlaceabilityReport LocalRealigner::calculate_placeability(
    const std::vector<LocusEvidence>& evidence) {

    PlaceabilityReport report;

    if (evidence.empty()) {
        return report;
    }

    // [修正] 按 locus 聚合分数
    std::unordered_map<int32_t, std::vector<double>> locus_scores;

    for (const auto& ev : evidence) {
        if (ev.normalized_score >= 0.1) {
            locus_scores[ev.locus_pos].push_back(ev.normalized_score);
        }
    }

    if (locus_scores.empty()) {
        return report;
    }

    // 计算每个 locus 的综合分数
    struct LocusScore {
        int32_t pos;
        double avg_score;
        double max_score;
        int support_count;
        double combined_score;
    };

    std::vector<LocusScore> scored_loci;
    scored_loci.reserve(locus_scores.size());

    for (const auto& [pos, scores] : locus_scores) {
        LocusScore ls;
        ls.pos = pos;
        ls.support_count = static_cast<int>(scores.size());

        // 平均分
        double sum = 0.0;
        double max_s = 0.0;
        for (double s : scores) {
            sum += s;
            max_s = std::max(max_s, s);
        }
        ls.avg_score = sum / scores.size();
        ls.max_score = max_s;

        // 综合分数：加权平均 + 支持数量 bonus
        ls.combined_score = ls.avg_score * 0.6 + ls.max_score * 0.4;
        ls.combined_score *= std::log2(ls.support_count + 1);

        scored_loci.push_back(ls);
    }

    // 按综合分数排序
    std::sort(scored_loci.begin(), scored_loci.end(),
        [](const LocusScore& a, const LocusScore& b) {
            return a.combined_score > b.combined_score;
        });

    // 填充报告
    report.best_score = scored_loci[0].combined_score;
    report.best_normalized = scored_loci[0].avg_score;
    report.best_locus_pos = scored_loci[0].pos;
    report.best_support_count = scored_loci[0].support_count;

    if (scored_loci.size() > 1) {
        report.second_score = scored_loci[1].combined_score;
        report.second_normalized = scored_loci[1].avg_score;
        report.delta_score = report.best_score - report.second_score;
    } else {
        report.second_score = 0.0;
        report.second_normalized = 0.0;
        report.delta_score = report.best_score;
    }

    // [修正] 置信度计算：综合 delta、best score、支持数量
    double delta_factor = std::min(1.0, report.delta_score /
        std::max(report.best_score * 0.5, 1.0));

    double score_factor = std::min(1.0, report.best_normalized);

    double support_factor = std::min(1.0,
        std::log2(scored_loci[0].support_count + 1) / 4.0);

    report.confidence = static_cast<float>(
        delta_factor * 0.4 + score_factor * 0.4 + support_factor * 0.2);

    report.confidence = std::clamp(report.confidence, 0.0f, 1.0f);

    // 唯一性：如果第二名也很强，降低置信度
    if (scored_loci.size() > 1) {
        double ratio = report.second_score /
            std::max(report.best_score, 1e-6);
        if (ratio > 0.8) {
            // 第二名接近第一名，不唯一
            report.confidence *= static_cast<float>(1.0 - ratio);
            report.is_ambiguous = true;
        }
    }

    return report;
}

// ============================================================================
// Tier Determination
// [修正] 完整的分层逻辑
// ============================================================================

int LocalRealigner::determine_tier_(const PlaceabilityReport& report) const {
    // Tier 1: 高置信度、唯一定位、链平衡
    if (report.confidence >= 0.8 &&
        report.best_normalized >= 0.7 &&
        report.delta_score >= report.best_score * 0.3 &&
        report.strand_balanced &&
        report.strand_ratio >= 0.2 &&
        !report.is_ambiguous) {
        return 1;
    }

    // Tier 2: 中等置信度，或链不平衡但其他指标好
    if (report.confidence >= 0.5 &&
        report.best_normalized >= 0.5 &&
        report.delta_score >= report.best_score * 0.15) {
        return 2;
    }

    // Tier 3: 低置信度但有一些证据
    if (report.confidence >= 0.2 &&
        report.best_normalized >= 0.3 &&
        report.best_support_count >= 2) {
        return 3;
    }

    // Tier 4: 非常弱的证据
    if (report.best_score > 0) {
        return 4;
    }

    // Tier 5: 无有效证据
    return 5;
}

// ============================================================================
// Thread Pool Implementation
// ============================================================================

namespace {

class ThreadPool {
public:
    explicit ThreadPool(size_t num_threads) : stop_(false) {
        size_t n = num_threads > 0 ?
            num_threads : std::thread::hardware_concurrency();
        n = std::max(n, size_t(1));

        workers_.reserve(n);
        for (size_t i = 0; i < n; ++i) {
            workers_.emplace_back([this] {
                while (true) {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(queue_mutex_);
                        condition_.wait(lock, [this] {
                            return stop_ || !tasks_.empty();
                        });
                        if (stop_ && tasks_.empty()) return;
                        task = std::move(tasks_.front());
                        tasks_.pop();
                    }
                    task();
                }
            });
        }
    }

    ~ThreadPool() {
        {
            std::unique_lock<std::mutex> lock(queue_mutex_);
            stop_ = true;
        }
        condition_.notify_all();
        for (auto& t : workers_) {
            if (t.joinable()) t.join();
        }
    }

    // 禁止拷贝
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    template<typename F>
    auto enqueue(F&& f) -> std::future<decltype(f())> {
        using R = decltype(f());
        auto task = std::make_shared<std::packaged_task<R()>>(
            std::forward<F>(f));
        std::future<R> future = task->get_future();
        {
            std::unique_lock<std::mutex> lock(queue_mutex_);
            if (stop_) {
                throw std::runtime_error(
                    "enqueue on stopped ThreadPool");
            }
            tasks_.emplace([task]() { (*task)(); });
        }
        condition_.notify_one();
        return future;
    }

    size_t size() const { return workers_.size(); }

private:
    std::vector<std::thread> workers_;
    std::queue<std::function<void()>> tasks_;
    std::mutex queue_mutex_;
    std::condition_variable condition_;
    bool stop_;
};

}  // anonymous namespace

// ============================================================================
// Batch Processing
// [修正] 安全的引用捕获 + 异常处理
// ============================================================================

std::vector<PlaceabilityReport> LocalRealigner::realign_batch(
    std::vector<Component>& components,
    const std::vector<ReadSketch>& reads,
    const GenomeAccessor& genome) {

    std::vector<PlaceabilityReport> reports;
    reports.resize(components.size());

    if (components.empty()) return reports;

    // 小批量直接串行处理
    if (components.size() <= 4 || config_.num_threads <= 1) {
        for (size_t i = 0; i < components.size(); ++i) {
            reports[i] = realign_component(
                components[i], reads, genome);
        }
        return reports;
    }

    // [修正] 使用索引而非引用捕获，避免悬空引用风险
    ThreadPool pool(config_.num_threads);
    std::vector<std::future<PlaceabilityReport>> futures;
    futures.reserve(components.size());

    for (size_t i = 0; i < components.size(); ++i) {
        futures.push_back(pool.enqueue(
            [this, i, &components, &reads, &genome]()
                -> PlaceabilityReport {
                try {
                    return realign_component(
                        components[i], reads, genome);
                } catch (const std::exception& e) {
                    // 异常安全：返回空报告
                    PlaceabilityReport empty_report;
                    empty_report.tier = 5;
                    return empty_report;
                }
            }));
    }

    // 收集结果（保持顺序）
    for (size_t i = 0; i < futures.size(); ++i) {
        if (futures[i].valid()) {
            reports[i] = futures[i].get();
        } else {
            reports[i].tier = 5;
        }
    }

    return reports;
}

// ============================================================================
// Process Component（内部包装）
// ============================================================================

PlaceabilityReport LocalRealigner::process_component_(
    Component& component,
    const std::vector<ReadSketch>& reads,
    const GenomeAccessor& genome) {

    return realign_component(component, reads, genome);
}

// ============================================================================
// realign_and_collect: 真正的比对 + evidence 收集
// [新增] 使用真实序列比对而非距离评分
// ============================================================================

RealignResult LocalRealigner::realign_and_collect(
    Component& component,
    const std::vector<ReadSketch>& reads,
    const GenomeAccessor& genome) {

    RealignResult result;

    // Step 1: 填充候选位点
    populate_locus_set(component, reads);

    if (component.locus_set.empty()) {
        component.locus_set.push_back({
            component.chrom_tid,
            component.centroid,
            1.0,
            static_cast<uint32_t>(component.read_count),
            0x04
        });
    }

    // Step 2: 排序并限制候选位点数量
    component.locus_set = rank_loci(
        std::move(component.locus_set), config_);

    const size_t max_loci = std::min(
        static_cast<size_t>(config_.max_locus_per_component),
        component.locus_set.size());
    if (component.locus_set.size() > max_loci) {
        component.locus_set.resize(max_loci);
    }

    // Step 3: 限制 reads 数量
    std::vector<size_t> sampled_reads;
    const size_t max_reads = std::min(
        static_cast<size_t>(config_.max_reads_per_component),
        component.read_indices.size());

    for (size_t i = 0;
         i < component.read_indices.size() && sampled_reads.size() < max_reads;
         ++i) {
        sampled_reads.push_back(component.read_indices[i]);
    }

    // Step 4: 对每个 read × 每个 locus 生成证据（使用真正比对）
    result.evidence.reserve(sampled_reads.size() * component.locus_set.size());

    for (size_t read_idx : sampled_reads) {
        if (read_idx >= reads.size()) continue;

        const auto& read = reads[read_idx];

        // 提取 read 的侧翼序列
        auto [up_flank, down_flank] = extract_flanks_view(
            read, config_.flank_length, genome);

        for (const auto& locus : component.locus_set) {
            LocusEvidence ev;
            ev.read_idx = read_idx;
            ev.locus_pos = locus.pos;
            ev.is_reverse = (read.flag & 0x10) != 0;
            ev.evidence_bits = 0;

            // 真正的序列对齐评分
            bool has_alignment = false;

            // 提取 locus 周围的参考序列
            int32_t ref_start = std::max<int32_t>(
                0, locus.pos - config_.flank_length);
            int32_t ref_end = locus.pos + config_.flank_length;

            auto ref_opt = genome.fetch(
                component.chrom_tid, ref_start, ref_end);

            if (ref_opt.has_value() && !read.sequence.empty()) {
                std::string_view ref_seq = ref_opt->sv();
                std::string_view read_seq(
                    read.sequence.data(), read.sequence.size());

                int32_t read_offset = locus.pos - read.pos;
                int32_t extract_start = std::max<int32_t>(
                    0, read_offset - config_.flank_length);
                int32_t extract_end = std::min<int32_t>(
                    static_cast<int32_t>(read_seq.size()),
                    read_offset + config_.flank_length);

                if (extract_start < extract_end &&
                    extract_end <= static_cast<int32_t>(read_seq.size())) {

                    std::string_view read_segment = read_seq.substr(
                        extract_start, extract_end - extract_start);

                    if (!read_segment.empty() && !ref_seq.empty()) {
                        AlignmentResult aln = dispatch_align_(
                            read_segment, ref_seq);

                        if (aln.score > -1e8) {
                            ev.normalized_score = std::max(0.0, aln.normalized_score);
                            ev.total_score = aln.score;
                            ev.weighted_identity = aln.identity;
                            has_alignment = true;

                            if (aln.identity >= config_.min_identity) {
                                ev.evidence_bits |= 0x01;
                            }
                            if (aln.has_indel) {
                                ev.evidence_bits |= 0x08;
                            }
                        }
                    }
                }
            }

            // 如果无法做序列对齐，回退到距离评分
            if (!has_alignment) {
                double dist = std::abs(
                    static_cast<double>(read.pos) - locus.pos);
                ev.normalized_score = std::max(
                    0.0, 1.0 - dist / config_.search_window);
                ev.total_score = ev.normalized_score * 100;
                ev.weighted_identity = 0.0f;

                if (ev.normalized_score > config_.min_normalized_score) {
                    ev.evidence_bits |= 0x02;
                }
            }

            // SA 证据加成
            if (locus.evidence_mask & 0x02) {
                ev.evidence_bits |= 0x04;
            }

            result.evidence.push_back(ev);
        }
    }

    // Step 5: 计算 placeability
    result.report = calculate_placeability(result.evidence);

    // Step 6: 链统计
    result.report.forward_count = 0;
    result.report.reverse_count = 0;

    for (const auto& ev : result.evidence) {
        if (ev.normalized_score < 0.1) continue;
        if (ev.is_reverse) ++result.report.reverse_count;
        else ++result.report.forward_count;
    }

    result.report.strand_balanced =
        (result.report.forward_count > 0 && result.report.reverse_count > 0);

    if (result.report.forward_count + result.report.reverse_count > 0) {
        int min_strand = std::min(
            result.report.forward_count, result.report.reverse_count);
        int max_strand = std::max(
            result.report.forward_count, result.report.reverse_count);
        result.report.strand_ratio = static_cast<float>(min_strand) /
            std::max(max_strand, 1);
    }

    // Step 7: 确定 tier
    result.report.tier = determine_tier_(result.report);

    // Step 8: 记录最佳 locus
    if (!component.locus_set.empty() && result.report.best_score > 0) {
        result.report.best_locus = 0;
        result.report.best_locus_pos = component.locus_set[0].pos;
    }

    return result;
}

}  // namespace placer


