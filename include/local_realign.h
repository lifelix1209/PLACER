#ifndef PLACER_LOCALREALIGN_H
#define PLACER_LOCALREALIGN_H

#include "component_builder.h"
#include "bam_reader.h"
#include <vector>
#include <string>
#include <cstdint>
#include <string_view>
#include <optional>
#include <array>
#include <concepts>
#include <memory>
#include <htslib/faidx.h>

namespace placer {

/**
 * 2-bit encoded sequence for memory efficiency
 * A=0, C=1, G=2, T=3, N=4(invalid)
 */
class BitpackedSeq {
public:
    static constexpr uint8_t INVALID = 4;
    static constexpr uint8_t MASK = 0x03;

    BitpackedSeq() = default;

    explicit BitpackedSeq(std::string_view ascii_seq) {
        encode(ascii_seq);
    }

    void encode(std::string_view ascii_seq);  // Declaration - implementation in .cpp

    size_t size() const { return size_; }
    const std::vector<uint64_t>& data() const { return data_; }

    static uint8_t ascii_to_2bit(char c) {
        switch (c) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': return 3;
            default: return INVALID;
        }
    }

    static char bit2_to_ascii(uint8_t code) {
        static constexpr char map[5] = {'A', 'C', 'G', 'T', 'N'};
        return map[code & 0x07];
    }

    // Get 2-bit code at position
    uint8_t operator[](size_t idx) const {
        size_t word_idx = idx / 32;
        int bit_idx = (idx % 32) * 2;
        return static_cast<uint8_t>((data_[word_idx] >> bit_idx) & MASK);
    }

    // Compare two sequences (Hamming distance)
    size_t hamming_distance(const BitpackedSeq& other) const {
        size_t count = 0;
        size_t min_len = std::min(size_, other.size_);
        for (size_t i = 0; i < min_len; ++i) {
            if ((*this)[i] != other[i]) ++count;
        }
        return count;
    }

    // Convert to string (for debugging)
    std::string to_string() const {
        std::string result;
        result.reserve(size_);
        for (size_t i = 0; i < size_; ++i) {
            result.push_back(bit2_to_ascii((*this)[i]));
        }
        return result;
    }

private:
    std::vector<uint64_t> data_;
    size_t size_ = 0;
};

/**
 * Zero-copy sequence view (no allocation)
 */
struct SeqView {
    const char* data = nullptr;
    size_t length = 0;

    SeqView() = default;
    SeqView(const char* d, size_t l) : data(d), length(l) {}

    explicit SeqView(std::string_view sv) : data(sv.data()), length(sv.size()) {}

    bool empty() const { return length == 0; }
    char operator[](size_t i) const { return data[i]; }

    std::string_view sv() const { return std::string_view(data, length); }
};

/**
 * Reference genome accessor with FAI index support
 * Uses manual file access (portable, no htslib dependency)
 */
class GenomeAccessor {
public:
    struct IndexEntry {
        std::string name;
        uint64_t offset;
        uint64_t length;
        int line_len;
        int line_blen;
    };

    explicit GenomeAccessor(std::string_view fasta_path);
    ~GenomeAccessor() = default;

    // Move semantics
    GenomeAccessor(GenomeAccessor&&) noexcept;
    GenomeAccessor& operator=(GenomeAccessor&&) noexcept;

    // Delete copy
    GenomeAccessor(const GenomeAccessor&) = delete;
    GenomeAccessor& operator=(const GenomeAccessor&) = delete;

    // Get reference sequence as view
    std::optional<SeqView> fetch(int chrom_tid, int32_t start, int32_t end) const;

    // Get full chromosome as BitpackedSeq
    std::optional<BitpackedSeq> fetch_packed(int chrom_tid) const;

    const IndexEntry& get_index(int chrom_tid) const { return index_[chrom_tid]; }
    size_t num_chroms() const { return index_.size(); }

private:
    std::string fasta_path_;
    std::vector<IndexEntry> index_;

    bool load_fai(std::string_view fasta_path);
    std::optional<SeqView> fetch_region(int chrom_tid, int32_t start, int32_t end) const;
};

/**
 * Industrial-grade alignment result
 */
struct AlignmentResult {
    double score = 0.0;           // Alignment score
    double normalized_score = 0.0; // Score per base
    float identity = 0.0f;        // Similarity (0-1)
    int matches = 0;
    int mismatches = 0;
    int gap_opens = 0;
    int gap_extensions = 0;
    int32_t query_start = 0;
    int32_t query_end = 0;
    int32_t target_start = 0;
    int32_t target_end = 0;
    bool has_indel = false;

    // Quality metrics
    float mapq_equiv = 0.0f;  // Equivalent MAPQ

    bool is_valid() const { return score > -1e9; }
};

/**
 * Locus evidence with full support metrics
 */
struct LocusEvidence {
    size_t read_idx = 0;
    int32_t locus_pos = 0;

    // Alignment metrics
    double up_score = 0.0;
    double down_score = 0.0;
    double total_score = 0.0;
    double normalized_score = 0.0;

    // Identity breakdown
    float up_identity = 0.0f;
    float down_identity = 0.0f;
    float weighted_identity = 0.0f;

    // Strand and direction
    bool is_reverse = false;
    bool consistent_strand = true;

    // Match quality
    int up_match_bases = 0;
    int down_match_bases = 0;
    int total_match_bases = 0;

    // Evidence strength
    uint32_t evidence_bits = 0;  // bit0=up, bit1=down, bit2=indel, bit3=high_id
};

/**
 * Placeability assessment for a component
 */
struct PlaceabilityReport {
    double delta_score = 0.0;      // Best - Second Best
    double best_normalized = 0.0;
    double second_normalized = 0.0;
    float confidence = 0.0f;        // 0-1 confidence score

    int best_locus = -1;
    double best_score = 0.0;
    double second_score = 0.0;

    // Strand consistency
    int forward_count = 0;
    int reverse_count = 0;
    bool strand_balanced = false;

    // Tier determination
    int tier = 3;  // 1, 2, or 3
};

/**
 * Industrial-grade realignment configuration
 */
struct RealignConfig {
    // Flank parameters
    int flank_length = 1000;        // Target flank length

    // Search space
    int32_t search_window = 10000; // ±N bp around candidate

    // Alignment parameters (affine gap)
    int8_t match = 2;
    int8_t mismatch = -3;
    int8_t gap_open = -5;
    int8_t gap_extend = -1;

    // Filtering
    float min_identity = 0.85f;
    double min_normalized_score = 0.3;
    float min_mapq_equiv = 10.0f;

    // Limits
    uint32_t max_locus_per_component = 20;
    uint32_t max_reads_per_component = 100;

    // Performance
    bool use_simd = true;          // Enable SIMD if available
    int num_threads = 4;            // Parallel threads
};

/**
 * Industrial-grade Local Realigner
 *
 * Design principles:
 * - Zero-copy sequence access via views
 * - SIMD-accelerated alignment (KSW2/Edlib interface)
 * - Batch processing for parallel execution
 * - Proper N-base and masking handling
 */
class LocalRealigner {
public:
    explicit LocalRealigner(RealignConfig config = RealignConfig());

    // =========================================================================
    // Batch Processing (Industrial-grade)
    // =========================================================================

    /**
     * Process multiple components in parallel
     */
    std::vector<PlaceabilityReport> realign_batch(
        std::vector<Component>& components,
        const std::vector<ReadSketch>& reads,
        const GenomeAccessor& genome);

    /**
     * Process single component (thread-safe)
     */
    PlaceabilityReport realign_component(
        Component& component,
        const std::vector<ReadSketch>& reads,
        const GenomeAccessor& genome);

    // =========================================================================
    // Core Alignment Interface
    // =========================================================================

    /**
     * SIMD-accelerated alignment (interface to KSW2/Edlib)
     * Returns optimal alignment score and metrics
     */
    static AlignmentResult align_sequences(
        std::string_view query,
        std::string_view target,
        const RealignConfig& config);

    /**
     * Fast Hamming distance for equal-length sequences
     * Uses SIMD when available
     */
    static int hamming_distance(std::string_view a, std::string_view b);

    /**
     * Seed finding with minimizer support
     */
    static std::vector<uint32_t> find_seeds(
        std::string_view seq,
        int kmer_size = 10,
        int stride = 1);

    // =========================================================================
    // Locus Management
    // =========================================================================

    /**
     * Populate locus set from multiple evidence sources:
     * - Primary alignment positions
     * - SA tag split points
     * - Discordant pairs (if available)
     * - TE k-mer reverse index hits
     */
    void populate_locus_set(
        Component& component,
        const std::vector<ReadSketch>& reads);

    /**
     * Filter and rank loci by evidence strength
     */
    static std::vector<LocusCandidate> rank_loci(
        std::vector<LocusCandidate>&& loci,
        const RealignConfig& config);

    // =========================================================================
    // Utility
    // =========================================================================

    /**
     * Extract flanks as views (zero-copy)
     */
    static std::pair<SeqView, SeqView> extract_flanks_view(
        const ReadSketch& read,
        int flank_length,
        const GenomeAccessor& genome);

    /**
     * Calculate placeability score (Δ = best - second best)
     */
    static PlaceabilityReport calculate_placeability(
        const std::vector<LocusEvidence>& evidence);

    const RealignConfig& config() const { return config_; }

private:
    RealignConfig config_;

    // Internal alignment dispatcher (selects optimal algorithm)
    AlignmentResult dispatch_align_(
        std::string_view query,
        std::string_view target);

    // SIMD-accelerated path
    AlignmentResult simd_align_(
        std::string_view query,
        std::string_view target);

    // Fallback scalar path
    AlignmentResult scalar_align_(
        std::string_view query,
        std::string_view target);

    // Process evidence for a component
    PlaceabilityReport process_component_(
        Component& component,
        const std::vector<ReadSketch>& reads,
        const GenomeAccessor& genome);

    // Tier determination
    int determine_tier_(const PlaceabilityReport& report) const;
};

// =========================================================================
// Inline implementations for performance
// =========================================================================

inline int LocalRealigner::determine_tier_(const PlaceabilityReport& report) const {
    // Tier 1: High confidence unique placement
    if (report.delta_score > 30.0 && report.confidence > 0.9f &&
        report.best_normalized > 0.5 && report.strand_balanced) {
        return 1;
    }

    // Tier 2: Multiple placements but consistent structure
    if (report.delta_score > 10.0 && report.confidence > 0.6f) {
        return 2;
    }

    // Tier 3: Low confidence or inconsistent
    return 3;
}

}  // namespace placer

#endif  // PLACER_LOCALREALIGN_H
