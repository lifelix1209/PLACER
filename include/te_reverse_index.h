#ifndef PLACER_TE_REVERSE_INDEX_H
#define PLACER_TE_REVERSE_INDEX_H

#include "gate1.h"
#include "component_builder.h"
#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <array>
#include <memory>
#include <limits>

namespace placer {

// ============================================================================
// TE Reverse Index Configuration
// ============================================================================

struct TEReverseIndexConfig {
    // k-mer 参数（与 Gate1 一致）
    uint8_t kmer_size = 15;           // k-mer 长度
    uint8_t minimizer_window = 10;      // minimizer 窗口大小

    // 索引参数
    int max_genome_hits_per_kmer = 20; // 每个 k-mer 在基因组中最大返回位置数
    int min_genome_hits = 3;           // 最小基因组命中数（用于过滤）

    // 聚类参数
    int32_t locus_cluster_radius = 50; // 位点聚类半径 (bp)
    int min_locus_reads = 2;          // 最小支持 reads 数

    // 召回参数
    int min_probe_hits = 3;           // 最小探针命中数
    double min_hit_density = 0.1;       // 最小命中密度

    // 性能参数
    int num_threads = 4;
    size_t max_index_memory_gb = 8;    // 最大内存使用 (GB)

    // Gate1 配置（用于探针提取）
    Gate1Config gate1_config;
};

// ============================================================================
// Genome K-mer Index（基因组 k-mer 索引，用于反向查找）
// ============================================================================

/**
 * GenomeKmerIndex: 基因组 k-mer 索引
 *
 * 目的：从 TE 探针序列快速查找基因组中的同源位置
 * 实现：对参考基因组构建 minimizer 索引
 *       查询时对探针提取 minimizers，返回基因组位置列表
 */
class GenomeKmerIndex {
public:
    explicit GenomeKmerIndex(const TEReverseIndexConfig& config);
    ~GenomeKmerIndex();

    // 从 FASTA 构建索引
    bool build_from_fasta(const std::string& fasta_path);

    // 从序列构建索引
    void build_from_sequences(const std::vector<std::string>& sequences);

    // 查询：返回基因组位置列表
    // 返回 (chrom_tid, position) 对的向量
    std::vector<std::pair<int32_t, int32_t>> query(
        const std::string& sequence) const;

    // 查询：返回命中计数
    int query_hit_count(const std::string& sequence) const;

    // 获取基因组大小
    int64_t genome_size() const { return total_genome_bases_; }

    // 获取染色体数
    int num_chromosomes() const { return static_cast<int>(chrom_starts_.size()); }

    // 获取染色体起始位置
    int64_t chrom_start(int32_t tid) const {
        return (tid >= 0 && tid < static_cast<int32_t>(chrom_starts_.size()))
            ? chrom_starts_[tid] : -1;
    }

    // 检查是否已构建
    bool empty() const { return kmer_to_positions_.empty(); }

private:
    TEReverseIndexConfig config_;

    // 2-bit 编码
    static inline uint8_t char_to_2bit(char c) {
        switch (c) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': return 3;
            default: return 4;  // N 或其他
        }
    }

    static inline uint64_t roll_kmer(uint64_t prev, uint8_t outgoing, uint8_t incoming, uint64_t mask) {
        return ((prev << 2) | incoming) & mask;
    }

    // 提取 minimizer
    uint64_t extract_minimizer(const std::string& seq, size_t& pos) const;

    // 添加序列到索引
    void add_sequence(const std::string& seq, int32_t chrom_tid, int64_t offset);

    // 索引结构：k-mer -> 位置列表
    // 使用 vector 存储多个位置，避免 unordered_set 的内存开销
    std::unordered_map<uint64_t, std::vector<int64_t>> kmer_to_positions_;

    // Reservoir sampling 计数器：记录每个 k-mer 看到的总位置数
    std::unordered_map<uint64_t, size_t> kmer_total_counts_;

    // 位置编码：(tid << 40) | position，用于压缩存储
    // 使用 40 bits 表示位置（足够 1TB 基因组）

    // 染色体信息
    std::vector<int64_t> chrom_starts_;  // 每个染色体的起始偏移
    std::vector<std::string> chrom_names_;
    int64_t total_genome_bases_ = 0;

    // k-mer 掩码
    uint64_t kmer_mask_ = 0;
};

// ============================================================================
// Rescued Locus Candidate
// ============================================================================

/**
 * RescuedLocus: 通过 TE 反向索引挽救的候选位点
 */
struct RescuedLocus {
    int32_t chrom_tid = -1;
    int32_t position = -1;           // 基因组位置（取中位数）
    int32_t cluster_id = -1;          // 聚类 ID

    // 支持信息
    int support_reads = 0;           // 支持的 reads 数
    int total_hits = 0;             // 总命中数
    double hit_density = 0.0;        // 命中密度

    // 探针来源
    std::vector<size_t> probe_read_indices;  // 支持的 read 索引
    std::vector<int32_t> probe_positions;     // 探针在 read 上的位置

    // [修正] genome 坐标列表（用于 merge 时的 position 和 density 计算）
    std::vector<int32_t> genome_hit_positions;

    // cluster span（用于 merge）
    int32_t cluster_start = 0;
    int32_t cluster_end = 0;

    // 评估信息
    double placeability_score = 0.0;
    bool passed_filter = false;
};

// ============================================================================
// Read Rescue Results
// ============================================================================

/**
 * ReadRescueResult: 单个 Read 的挽救结果
 */
struct ReadRescueResult {
    size_t read_idx = 0;

    // 探针信息
    std::string probe_sequence;
    int probe_start = -1;
    int probe_end = -1;

    // 基因组命中
    std::vector<std::pair<int32_t, int32_t>> genome_hits;  // (tid, position)

    // 是否成功挽救
    bool rescued = false;
    int num_valid_hits = 0;
};

// ============================================================================
// Hit 结构体：统一的命中记录（内部使用）
// ============================================================================

struct Hit {
    int32_t tid;
    int32_t pos;          // 估计的 read 在 genome 上的起点
    size_t  read_idx;
    int32_t probe_offset; // probe 在 read 中的偏移
};

// ============================================================================
// TE Reverse Index（主类）
// ============================================================================

/**
 * TEReverseIndex: TE 反向索引
 *
 * 目的：
 * 1. 对没有 Secondary/Supplementary Alignment 的 Reads 查找候选位点
 * 2. 使用 TE 探针序列查询基因组 k-mer 索引
 * 3. 聚类命中位置，生成候选 locus 集合
 *
 * 数据流：
 * Read (无 SA) → 提取探针 → GenomeKmerIndex 查询 → 聚类 → RescuedLocus
 */
class TEReverseIndex {
public:
    explicit TEReverseIndex(const TEReverseIndexConfig& config);
    ~TEReverseIndex();

    // 初始化：从 TE FASTA 构建基因组索引
    bool initialize(const std::string& genome_fasta);

    // 批量处理：挽救无 SA 的 reads
    // 返回挽救的候选位点列表
    std::vector<RescuedLocus> rescue_reads(
        const std::vector<ReadSketch>& reads,
        const std::vector<Component>& existing_components,
        const std::vector<bool>& has_sa_flags);

    // 批量处理：带有证据的 reads
    std::vector<RescuedLocus> rescue_with_evidence(
        const std::vector<ReadSketch>& reads,
        const std::vector<std::vector<ProbeFragment>>& all_probes);

    // 集成到现有组件：补充 locus 集合
    // 返回需要更新的组件 ID 列表
    std::vector<int32_t> integrate_with_components(
        std::vector<Component>& components,
        const std::vector<ReadSketch>& reads);

    // 添加外部 locus 到索引（用于跨组件共享）
    void add_external_loci(const std::vector<std::pair<int32_t, int32_t>>& loci);

    // 统计信息
    struct Stats {
        int64_t total_reads_processed = 0;
        int64_t reads_with_sa = 0;
        int64_t reads_rescued = 0;
        int64_t total_rescued_loci = 0;
        int64_t loci_merged = 0;
        int64_t loci_discarded = 0;
        double avg_hits_per_rescued_read = 0.0;
    };
    const Stats& stats() const { return stats_; }

    // 重置统计
    void reset_stats() { stats_ = Stats(); }

private:
    TEReverseIndexConfig config_;
    std::unique_ptr<GenomeKmerIndex> genome_index_;
    Gate1 gate1_;  // 用于提取探针
    Stats stats_;

    // 从探针收集 hits
    std::vector<Hit> collect_hits_from_probes(
        const std::vector<ProbeFragment>& probes,
        size_t read_idx) const;

    // 聚类命中位置（新版：使用 Hit 结构体）
    std::vector<RescuedLocus> cluster_hits(
        const std::vector<Hit>& hits);

    // 合并两个相邻的 loci
    void merge_loci(
        RescuedLocus& target,
        const RescuedLocus& source) const;

    // 评估挽救的位点
    void evaluate_locus(RescuedLocus& locus) const;

    // 过滤低质量位点
    std::vector<RescuedLocus> filter_loci(
        std::vector<RescuedLocus>& loci);

    // 提取 Read 的探针（如果有）
    bool extract_probes_from_read(
        const ReadSketch& read,
        std::vector<ProbeFragment>& probes) const;

    // 并行处理辅助函数
    std::vector<ReadRescueResult> rescue_batch(
        const std::vector<ReadSketch>& reads,
        const std::vector<std::vector<ProbeFragment>>& all_probes,
        size_t start, size_t end);
};

// ============================================================================
// Utility Functions
// ============================================================================

// 计算两个位置的距离（考虑染色体）
inline int32_t position_distance(
    int32_t chrom1, int32_t pos1,
    int32_t chrom2, int32_t pos2) {
    if (chrom1 != chrom2) return std::numeric_limits<int32_t>::max();
    return std::abs(pos1 - pos2);
}

// 聚类位置列表
std::vector<std::vector<std::pair<int32_t, int32_t>>> cluster_positions(
    const std::vector<std::pair<int32_t, int32_t>>& positions,
    int32_t radius);

}  // namespace placer

#endif  // PLACER_TE_REVERSE_INDEX_H
