#ifndef PLACER_BAM_READER_H
#define PLACER_BAM_READER_H

#include <string>
#include <vector>
#include <cstdint>
#include <functional>
#include <htslib/sam.h>

namespace placer {

/**
 * Read-level breakpoint information
 * Extracted from CIGAR soft-clips or split alignments
 */
struct ReadBreakpoint {
    // Breakpoint position in read coordinate (0-based)
    int32_t read_pos = -1;
    // Type of breakpoint signal
    enum Type {
        NONE = 0,
        SOFT_CLIP_5P = 1,   // 5' soft-clip boundary
        SOFT_CLIP_3P = 2,   // 3' soft-clip boundary
        SPLIT = 3,          // Split alignment breakpoint
        INSERTION = 4       // Large insertion boundary
    } type = NONE;
    // Reference coordinate at this breakpoint (if available)
    int32_t ref_pos = -1;
    // For INSERTION signal: inserted length from the triggering I-op
    int32_t insertion_len = 0;
    // Is this a valid evidence breakpoint?
    bool is_valid = false;
};

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

    // Read-level breakpoint (for INS extraction)
    ReadBreakpoint breakpoint;

    // MD tag presence
    bool has_md = false;
};

/**
 * ReadCallback: Function called for each read during streaming
 */
using ReadCallback = std::function<void(const ReadSketch&)>;

/**
 * StreamProgressCallback: Called periodically with progress info
 * @param processed Number of records processed so far
 * @param current_chrom Current chromosome being processed (-1 if unknown)
 * @return True to continue streaming, False to abort
 */
using StreamProgressCallback = std::function<bool(int64_t processed, int32_t current_chrom)>;

/**
 * BamReader: Single-pass BAM streaming reader
 * Reads BAM sequentially without full genome scan
 *
 * Supports two streaming modes:
 * 1. Callback-only: Just calls callback for each read
 * 2. Streaming mode: Integrates with WindowBuffer for immediate window processing
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
     * Stream with progress reporting
     * @param callback Function to call for each ReadSketch
     * @param progress_callback Called every N records (default: 100000)
     * @param progress_interval Number of records between progress updates
     * @return Number of records processed, -1 if invalid
     */
    int64_t stream_with_progress(
        const ReadCallback& callback,
        StreamProgressCallback progress_callback = nullptr,
        int64_t progress_interval = 100000);

    /**
     * Get header information
     */
    const std::string& get_bam_path() const { return bam_path_; }
    int64_t get_total_records() const { return total_records_; }

    /**
     * Get chromosome name from TID
     * @param tid Chromosome ID from BAM header
     * @return Chromosome name, or empty string if invalid
     */
    std::string get_chrom_name(int32_t tid) const;

    /**
     * Get number of chromosomes in BAM header
     */
    int32_t get_num_chromosomes() const { return header_ ? header_->n_targets : 0; }

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
