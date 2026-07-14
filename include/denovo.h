#ifndef PLACER_DENOVO_H
#define PLACER_DENOVO_H

#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include <htslib/sam.h>

namespace placer {

struct DenovoConfig {
    std::string child_scientific_path;
    std::string parent_bam_list_path;
    std::vector<std::string> parent_bam_paths;
    std::string reference_fasta_path;
    std::string te_fasta_path;
    std::string out_prefix = "trio_denovo";

    int32_t bam_threads = 1;
    int32_t child_min_support_reads = 2;

    int32_t fetch_window = 500;
    int32_t default_match_window = 30;
    int32_t max_match_window = 100;
    int32_t min_softclip_len = 50;
    int32_t min_insertion_len = 50;
    int32_t parent_mapq_min = 0;
    bool family_match_veto = true;
    bool emit_review_status = true;
    bool dry_run = false;
};

struct DenovoChildCandidate {
    std::string chrom;
    int32_t pos = -1;
    std::string te_name;
    int32_t support_reads = 0;
    int32_t event_start = -1;
    int32_t event_end = -1;
};

struct ParentReadEvidence {
    std::string parent_bam_path;
    std::string read_name;
    std::string evidence_type;
    int32_t breakpoint_pos = -1;
    int32_t fragment_len = 0;
    std::string matched_te = "NA";
    std::string matched_family = "NA";
    std::string veto_reason = "NA";
};

struct ParentEvidenceSummary {
    int32_t total_support_reads = 0;
    int32_t exact_te_reads = 0;
    int32_t family_te_reads = 0;
    int32_t ambiguous_signal_reads = 0;
    std::vector<std::string> support_bams;
    std::vector<ParentReadEvidence> veto_reads;
    std::string scanner_status = "NOT_RUN";
};

struct DenovoCall {
    DenovoChildCandidate child;
    ParentEvidenceSummary parent_summary;
    std::string status = "PENDING_PARENT_SCAN";
    std::string de_novo = "NA";
};

struct DenovoResult {
    std::string implementation_status = "SKELETON_PARENT_SCAN_NOT_IMPLEMENTED";
    int64_t child_rows_total = 0;
    int64_t child_candidates_considered = 0;
    int64_t parent_bams = 0;
    int64_t calls_written = 0;
    int64_t inherited_prefilter_hits = 0;
    int64_t parent_veto_calls = 0;
    int64_t review_calls = 0;
    int64_t denovo_pass_calls = 0;
    std::vector<DenovoCall> calls;
};

class IndexedBamReader {
public:
    explicit IndexedBamReader(std::string bam_path, int32_t decompression_threads = 1);
    ~IndexedBamReader();

    IndexedBamReader(const IndexedBamReader&) = delete;
    IndexedBamReader& operator=(const IndexedBamReader&) = delete;

    bool is_valid() const { return valid_; }
    bool has_index() const { return index_ != nullptr; }
    const std::string& bam_path() const { return bam_path_; }
    int32_t chromosome_count() const;
    std::string chromosome_name(int32_t tid) const;

    bool fetch(
        const std::string& chrom,
        int32_t start,
        int32_t end,
        const std::function<bool(const bam1_t*)>& handler) const;

private:
    std::string bam_path_;
    htsFile* file_ = nullptr;
    bam_hdr_t* header_ = nullptr;
    hts_idx_t* index_ = nullptr;
    bool valid_ = false;
};

class ParentPoolScanner {
public:
    explicit ParentPoolScanner(DenovoConfig config);
    ~ParentPoolScanner();

    ParentPoolScanner(const ParentPoolScanner&) = delete;
    ParentPoolScanner& operator=(const ParentPoolScanner&) = delete;

    const std::string& implementation_status() const { return implementation_status_; }

    ParentEvidenceSummary scan_candidate(const DenovoChildCandidate& candidate) const;

private:
    struct Impl;
    DenovoConfig config_;
    std::vector<std::unique_ptr<IndexedBamReader>> readers_;
    std::unique_ptr<Impl> impl_;
    std::string implementation_status_ = "SKELETON_PARENT_SCAN_NOT_IMPLEMENTED";
};

bool parse_denovo_cli_args(
    int argc,
    char** argv,
    DenovoConfig& config,
    std::string& error_message);

std::vector<DenovoChildCandidate> load_denovo_child_candidates(
    const DenovoConfig& config,
    int64_t* total_rows = nullptr);

DenovoResult run_denovo(const DenovoConfig& config);

void write_denovo_outputs(const DenovoConfig& config, const DenovoResult& result);

int run_denovo_cli(int argc, char** argv);

}  // namespace placer

#endif  // PLACER_DENOVO_H
