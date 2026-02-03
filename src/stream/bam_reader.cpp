#include "bam_reader.h"
#include <iostream>
#include <cstring>
#include <algorithm>

namespace placer {

// BAM 4-bit encoded sequence character lookup
static const char BAM_SEQ_CHARS[] = "=ACMGRSVTWYHKDBN";

BamReader::BamReader(const std::string& bam_path) : bam_path_(bam_path) {
    hts_file_ = hts_open(bam_path.c_str(), "r");
    if (!hts_file_) {
        std::cerr << "[Error] Failed to open BAM: " << bam_path << std::endl;
        valid_ = false;
        return;
    }

    // Enable multi-threaded decompression if available
    hts_set_threads(hts_file_, 4);

    header_ = sam_hdr_read(hts_file_);
    if (!header_) {
        std::cerr << "[Error] Failed to read BAM header" << std::endl;
        valid_ = false;
        return;
    }

    aln_ = bam_init1();
    valid_ = true;
}

BamReader::~BamReader() {
    if (aln_) bam_destroy1(aln_);
    if (header_) bam_hdr_destroy(header_);
    if (hts_file_) hts_close(hts_file_);
}

BamReader::BamReader(BamReader&& other) noexcept
    : bam_path_(std::move(other.bam_path_)),
      hts_file_(other.hts_file_),
      header_(other.header_),
      aln_(other.aln_),
      total_records_(other.total_records_),
      valid_(other.valid_) {
    other.hts_file_ = nullptr;
    other.header_ = nullptr;
    other.aln_ = nullptr;
    other.valid_ = false;
}

BamReader& BamReader::operator=(BamReader&& other) noexcept {
    if (this != &other) {
        // Release existing resources
        if (aln_) bam_destroy1(aln_);
        if (header_) bam_hdr_destroy(header_);
        if (hts_file_) hts_close(hts_file_);

        // Move resources from other
        bam_path_ = std::move(other.bam_path_);
        hts_file_ = other.hts_file_;
        header_ = other.header_;
        aln_ = other.aln_;
        total_records_ = other.total_records_;
        valid_ = other.valid_;

        // Invalidate other
        other.hts_file_ = nullptr;
        other.header_ = nullptr;
        other.aln_ = nullptr;
        other.valid_ = false;
    }
    return *this;
}

int64_t BamReader::stream(const ReadCallback& callback) {
    if (!valid_) {
        std::cerr << "[Error] BamReader is not valid, cannot stream" << std::endl;
        return -1;
    }

    total_records_ = 0;
    int ret;

    // Reuse ReadSketch to reduce memory allocation overhead
    ReadSketch sketch;
    sketch.sequence.reserve(20000);  // Typical long-read length
    sketch.cigar_ops.reserve(100);

    while ((ret = sam_read1(hts_file_, header_, aln_)) >= 0) {
        // Filter: only process primary alignments
        // Secondary alignments (0x100) are alternative mappings, skip them
        // Supplementary alignments (0x800) are split reads, keep them (critical for TE detection)
        if (aln_->core.flag & BAM_FSECONDARY) {
            continue;
        }

        // Optional: filter unmapped reads (can be disabled for some workflows)
        if (aln_->core.flag & BAM_FUNMAP) {
            continue;
        }

        // Optional: filter very low quality reads
        // if (aln_->core.qual < 5) continue;

        // Clear previous data, keeping allocated capacity
        sketch.cigar_ops.clear();
        sketch.sa_targets.clear();
        // sequence will be resized in extract_readsketch_fast

        extract_readsketch_fast(aln_, sketch);
        callback(sketch);
        total_records_++;
    }

    if (ret < -1) {
        std::cerr << "[Error] BAM read error at record " << total_records_ << std::endl;
    }

    return total_records_;
}

void BamReader::extract_readsketch_fast(const bam1_t* aln, ReadSketch& sketch) {
    // Basic alignment info
    sketch.qname = bam_get_qname(aln);
    sketch.tid = aln->core.tid;
    sketch.pos = aln->core.pos;
    sketch.flag = aln->core.flag;
    sketch.mapq = aln->core.qual;

    // CIGAR parsing with end_pos calculation
    uint32_t* cigar = bam_get_cigar(aln);
    int n_cigar = aln->core.n_cigar;

    sketch.cigar_ops.reserve(n_cigar);
    int64_t ref_consumed = 0;
    sketch.total_clip_len = 0;
    sketch.has_large_insertion = false;

    for (int i = 0; i < n_cigar; ++i) {
        int op = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);
        char op_char = bam_cigar_opchr(cigar[i]);

        sketch.cigar_ops.emplace_back(op_char, len);

        // Track soft-clip length
        if (op == BAM_CSOFT_CLIP) {
            sketch.total_clip_len += len;
        }
        // Track large insertions (>50bp threshold from dev_plan)
        else if (op == BAM_CINS && len > 50) {
            sketch.has_large_insertion = true;
        }

        // Calculate reference consumption
        if (bam_cigar_type(op) & 2) {
            ref_consumed += len;
        }
    }
    sketch.end_pos = sketch.pos + static_cast<int32_t>(ref_consumed);

    // Sequence extraction (required for Gate 1 TE-proxy)
    int32_t l_qseq = aln->core.l_qseq;
    if (l_qseq > 0) {
        uint8_t* seq_ptr = bam_get_seq(aln);
        sketch.sequence.resize(l_qseq);
        for (int i = 0; i < l_qseq; ++i) {
            sketch.sequence[i] = BAM_SEQ_CHARS[bam_seqi(seq_ptr, i)];
        }
    }

    // SA tag parsing
    uint8_t* sa_ptr = bam_aux_get(aln, "SA");
    if (sa_ptr) {
        sketch.has_sa = true;
        parse_sa_targets(aln, sketch);
    } else {
        sketch.has_sa = false;
    }

    // MD tag presence check
    sketch.has_md = (bam_aux_get(aln, "MD") != nullptr);
}

std::vector<std::pair<char, int>> BamReader::extract_cigar(const bam1_t* aln) const {
    std::vector<std::pair<char, int>> cigar_ops;

    int n_cigar = aln->core.n_cigar;
    uint32_t* cigar = bam_get_cigar(aln);

    for (int i = 0; i < n_cigar; ++i) {
        char op = bam_cigar_opchr(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);
        cigar_ops.emplace_back(op, len);
    }

    return cigar_ops;
}

void BamReader::parse_sa_targets(const bam1_t* aln, ReadSketch& sketch) {
    uint8_t* sa_ptr = bam_aux_get(aln, "SA");
    if (!sa_ptr) return;

    const char* sa_str = bam_aux2Z(sa_ptr);
    if (!sa_str) return;

    sketch.sa_targets.clear();
    const char* p = sa_str;

    // Reusable buffer for rname to avoid repeated allocations
    static std::string rname_buffer;

    while (*p) {
        // 1. Parse rname (Chromosome Name)
        const char* comma1 = strchr(p, ',');
        if (!comma1) break;

        size_t rname_len = comma1 - p;
        if (rname_len == 0) break;

        // Extract rname string for TID lookup
        rname_buffer.assign(p, rname_len);

        // Lookup TID via header (required for spatial prior in Phase 7)
        // sam_hdr_name2tid is a hash lookup, fast enough for sparse SA tags
        int32_t target_tid = sam_hdr_name2tid(header_, rname_buffer.c_str());

        // 2. Parse Position (SA is 1-based, convert to 0-based)
        p = comma1 + 1;
        char* end_ptr;
        long pos = strtol(p, &end_ptr, 10);

        // 3. Skip Strand, CIGAR, MapQ, NM (4 commas total)
        // Format: rname,pos,strand,CIGAR,mapQ,NM;
        p = end_ptr;
        for (int i = 0; i < 4; ++i) {
            const char* next = strchr(p, ',');
            if (!next) { p = nullptr; break; }
            p = next + 1;
        }
        if (!p) break;

        // 4. Skip NM to semicolon
        const char* semicolon = strchr(p, ';');

        // Only add valid TID entries (required for spatial prior calculation)
        if (target_tid >= 0) {
            sketch.sa_targets.emplace_back(target_tid, static_cast<int32_t>(pos - 1));
        }

        if (semicolon) {
            p = semicolon + 1;
        } else {
            break;
        }
    }
}

}  // namespace placer
