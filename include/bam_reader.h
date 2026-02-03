#ifndef PLACER_BAM_READER_H
#define PLACER_BAM_READER_H

#include <string>
#include <vector>
#include <cstdint>
#include <functional>
#include <htslib/sam.h>

namespace placer {

/**
 * ReadSketch v5: Lightweight read representation for streaming
 * Contains minimum info needed for downstream processing
 */
struct ReadSketch {
    std::string qname;                    // Read name (for primary/supplementary linking)
    int32_t tid;                           // Target ID (chromosome index)
    int32_t pos;                           // 0-based leftmost position
    int32_t end_pos;                       // 0-based rightmost position (calculated from CIGAR)
    uint16_t flag;                         // BAM flags
    uint8_t mapq;                          // Mapping quality

    // Gate 1 requirements
    std::string sequence;                  // Full read sequence for TE-proxy

    // Pre-extracted signals (avoid re-parsing CIGAR downstream)
    int32_t total_clip_len = 0;            // Total soft-clip length
    bool has_large_insertion = false;      // Any insertion > 50bp in CIGAR

    // CIGAR operations in compressed form
    std::vector<std::pair<char, int>> cigar_ops;

    // SA tag presence and targets (TID + pos pairs)
    bool has_sa = false;
    std::vector<std::pair<int32_t, int32_t>> sa_targets;  // (tid, pos) for each split

    // MD tag presence
    bool has_md = false;
};

using ReadCallback = std::function<void(const ReadSketch&)>;

/**
 * BamReader: Single-pass BAM streaming reader
 * Reads BAM sequentially without full genome scan
 */
class BamReader {
public:
    explicit BamReader(const std::string& bam_path);
    ~BamReader();

    // Disable copy, enable move
    BamReader(const BamReader&) = delete;
    BamReader& operator=(const BamReader&) = delete;
    BamReader(BamReader&& other) noexcept;
    BamReader& operator=(BamReader&& other) noexcept;

    /**
     * Stream through BAM, calling callback for each alignment
     * @param callback Function to call for each ReadSketch
     * @return Number of records processed, -1 if invalid
     */
    int64_t stream(const ReadCallback& callback);

    /**
     * Get header information
     */
    const std::string& get_bam_path() const { return bam_path_; }
    int64_t get_total_records() const { return total_records_; }

    /**
     * Check if file is valid
     */
    bool is_valid() const { return valid_; }

private:
    std::string bam_path_;
    htsFile* hts_file_ = nullptr;
    bam_hdr_t* header_ = nullptr;
    bam1_t* aln_ = nullptr;
    int64_t total_records_ = 0;
    bool valid_ = false;

    /**
     * Fast extraction of ReadSketch from BAM record
     */
    void extract_readsketch_fast(const bam1_t* aln, ReadSketch& sketch);

    /**
     * Extract CIGAR operations in compressed form
     */
    std::vector<std::pair<char, int>> extract_cigar(const bam1_t* aln) const;

    /**
     * Extract SA tag targets (tid, pos pairs)
     */
    void parse_sa_targets(const bam1_t* aln, ReadSketch& sketch);
};

}  // namespace placer

#endif  // PLACER_BAM_READER_H
