#include "pipeline.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <iostream>
#include <limits>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <htslib/sam.h>
#include <abpoa.h>

namespace placer {
namespace {

constexpr int32_t kSoftClipSignalMin = 20;
constexpr int32_t kLongInsertionSignalMin = 50;
constexpr int32_t kLargeIndelEvidenceMin = 40;
constexpr double kSaHintWeight = 0.35;
constexpr int32_t kHistBinSize = 20;
constexpr int32_t kHistPadding = 120;
constexpr double kPeakMinWeight = 1.2;
constexpr int32_t kPeakMergeDistance = 120;
constexpr int32_t kWindowMergeGap = 80;
constexpr int32_t kWindowMaxExpandBins = 30;
constexpr double kWindowDropRatio = 0.35;
constexpr int32_t kMinWindowSpan = 120;
constexpr double kReadAssignMinScore = 0.8;
constexpr double kReadAmbiguousRatio = 0.75;
constexpr double kReadAssignDecayBp = 150.0;

enum class EvidenceKind : uint8_t {
    kSoftClip = 0,
    kIndel = 1,
    kSAHint = 2
};

struct EvidencePoint {
    size_t read_index = 0;
    int32_t pos = -1;
    double weight = 0.0;
    EvidenceKind kind = EvidenceKind::kSoftClip;
    int32_t signal_len = 0;
    uint8_t class_mask = 0;
};

struct ReadSignalSummary {
    std::string read_id;
    bool is_reverse = false;
    bool has_sa_or_supp = false;
    int32_t max_soft_clip = 0;
    int32_t max_ins = 0;
    uint8_t class_mask = 0;
};

struct EvidenceBundle {
    std::vector<EvidencePoint> points;
    std::vector<ReadSignalSummary> read_summaries;
};

struct DensityPeak {
    int32_t pos = -1;
    double weight = 0.0;
};

struct CandidateWindow {
    int32_t start = -1;
    int32_t end = -1;
    int32_t center = -1;
    double peak_weight = 0.0;
};

const char* theta_tag(InsertionTheta theta) {
    switch (theta) {
        case InsertionTheta::kFwd: return "FWD";
        case InsertionTheta::kRev: return "REV";
        default: return "NA";
    }
}

double sigmoid(double x) {
    if (x >= 0.0) {
        const double z = std::exp(-x);
        return 1.0 / (1.0 + z);
    }
    const double z = std::exp(x);
    return z / (1.0 + z);
}

double log_binomial_likelihood(int32_t alt_count, int32_t ref_count, double alt_rate) {
    const double p = std::clamp(alt_rate, 1e-6, 1.0 - 1e-6);
    return (static_cast<double>(alt_count) * std::log(p)) +
           (static_cast<double>(ref_count) * std::log(1.0 - p));
}

int32_t median_i32(std::vector<int32_t> values) {
    if (values.empty()) {
        return -1;
    }
    std::sort(values.begin(), values.end());
    const size_t mid = values.size() / 2;
    if (values.size() % 2 == 1) {
        return values[mid];
    }
    return static_cast<int32_t>((static_cast<int64_t>(values[mid - 1]) + values[mid]) / 2);
}

double mad_from_positions(const std::vector<int32_t>& positions) {
    if (positions.empty()) {
        return 0.0;
    }
    const int32_t center = median_i32(positions);
    std::vector<int32_t> abs_dev;
    abs_dev.reserve(positions.size());
    for (int32_t pos : positions) {
        abs_dev.push_back(std::abs(pos - center));
    }
    return static_cast<double>(median_i32(std::move(abs_dev)));
}

bool is_softclip_source(InsertionFragmentSource source) {
    return source == InsertionFragmentSource::kClipRefLeft ||
           source == InsertionFragmentSource::kClipRefRight;
}

int32_t max_homopolymer_run(const std::string& seq) {
    if (seq.empty()) {
        return 0;
    }
    int32_t best = 1;
    int32_t run = 1;
    for (size_t i = 1; i < seq.size(); ++i) {
        if (seq[i] == seq[i - 1]) {
            ++run;
            best = std::max(best, run);
        } else {
            run = 1;
        }
    }
    return best;
}

double at_fraction(const std::string& seq) {
    int32_t at_count = 0;
    int32_t total = 0;
    for (char c : seq) {
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
            continue;
        }
        ++total;
        if (c == 'A' || c == 'T') {
            ++at_count;
        }
    }
    if (total <= 0) {
        return 0.0;
    }
    return static_cast<double>(at_count) / static_cast<double>(total);
}

bool is_low_complexity_softclip(
    const InsertionFragment& fragment,
    const PipelineConfig& config) {
    if (!is_softclip_source(fragment.source)) {
        return false;
    }
    if (fragment.sequence.empty()) {
        return true;
    }
    const double at_frac = at_fraction(fragment.sequence);
    const int32_t homopolymer = max_homopolymer_run(fragment.sequence);
    return at_frac >= config.te_softclip_low_complexity_at_frac_min ||
           homopolymer >= config.te_softclip_low_complexity_homopolymer_min;
}

uint8_t char_to_2bit(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return 4;
    }
}

bool build_kmer(const std::string& s, int32_t start, int32_t k, uint64_t& out) {
    out = 0;
    for (int32_t i = 0; i < k; ++i) {
        const uint8_t code = char_to_2bit(s[static_cast<size_t>(start + i)]);
        if (code > 3) {
            return false;
        }
        out = (out << 2) | code;
    }
    return true;
}

std::unordered_set<uint64_t> kmer_set(const std::string& seq, int32_t k) {
    std::unordered_set<uint64_t> out;
    if (k <= 0 || static_cast<int32_t>(seq.size()) < k) {
        return out;
    }
    out.reserve(static_cast<size_t>(seq.size()));
    for (int32_t i = 0; i + k <= static_cast<int32_t>(seq.size()); ++i) {
        uint64_t key = 0;
        if (!build_kmer(seq, i, k, key)) {
            continue;
        }
        out.insert(key);
    }
    return out;
}

double kmer_overlap_fraction(
    const std::unordered_set<uint64_t>& a,
    const std::unordered_set<uint64_t>& b) {
    if (a.empty() || b.empty()) {
        return 0.0;
    }
    const auto* smaller = &a;
    const auto* larger = &b;
    if (a.size() > b.size()) {
        smaller = &b;
        larger = &a;
    }
    int32_t inter = 0;
    for (uint64_t key : *smaller) {
        if (larger->find(key) != larger->end()) {
            ++inter;
        }
    }
    const int32_t denom = std::max<int32_t>(1, std::min<int32_t>(a.size(), b.size()));
    return static_cast<double>(inter) / static_cast<double>(denom);
}

double estimate_mean_identity_kmer(
    const std::string& consensus,
    const std::vector<std::string>& seqs,
    int32_t k) {
    if (consensus.empty() || seqs.empty()) {
        return 0.0;
    }
    const auto cons_kmers = kmer_set(consensus, k);
    if (cons_kmers.empty()) {
        return 0.0;
    }
    double sum = 0.0;
    int32_t used = 0;
    for (const auto& seq : seqs) {
        const auto seq_kmers = kmer_set(seq, k);
        if (seq_kmers.empty()) {
            continue;
        }
        sum += kmer_overlap_fraction(cons_kmers, seq_kmers);
        ++used;
    }
    if (used <= 0) {
        return 0.0;
    }
    return sum / static_cast<double>(used);
}

std::string upper_acgt(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        out.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
    }
    return out;
}

uint8_t nt_to_abpoa_code(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return 4;
    }
}

char abpoa_code_to_nt(uint8_t c) {
    switch (c) {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        default: return 'N';
    }
}

std::string abpoa_consensus(const std::vector<std::string>& seqs) {
    if (seqs.empty()) {
        return "";
    }

    abpoa_t* ab = abpoa_init();
    abpoa_para_t* abpt = abpoa_init_para();
    if (!ab || !abpt) {
        if (ab) {
            abpoa_free(ab);
        }
        if (abpt) {
            abpoa_free_para(abpt);
        }
        return "";
    }

    abpt->out_cons = 1;
    abpt->out_msa = 0;
    abpt->out_gfa = 0;
    abpt->max_n_cons = 1;
    abpoa_post_set_para(abpt);

    const int32_t n_seq = static_cast<int32_t>(seqs.size());
    std::vector<int> seq_lens(static_cast<size_t>(n_seq), 0);
    std::vector<std::vector<uint8_t>> packed(static_cast<size_t>(n_seq));
    std::vector<uint8_t*> seq_ptrs(static_cast<size_t>(n_seq), nullptr);

    for (int32_t i = 0; i < n_seq; ++i) {
        const std::string seq = upper_acgt(seqs[static_cast<size_t>(i)]);
        seq_lens[static_cast<size_t>(i)] = static_cast<int>(seq.size());
        packed[static_cast<size_t>(i)].resize(seq.size(), 4);
        for (size_t j = 0; j < seq.size(); ++j) {
            packed[static_cast<size_t>(i)][j] = nt_to_abpoa_code(seq[j]);
        }
        seq_ptrs[static_cast<size_t>(i)] = packed[static_cast<size_t>(i)].data();
    }

    const int ret = abpoa_msa(
        ab,
        abpt,
        n_seq,
        nullptr,
        seq_lens.data(),
        seq_ptrs.data(),
        nullptr,
        nullptr);

    std::string consensus;
    if (ret == 0 && ab->abc && ab->abc->n_cons > 0) {
        const int cons_len = ab->abc->cons_len[0];
        if (cons_len > 0 && ab->abc->cons_base && ab->abc->cons_base[0]) {
            consensus.reserve(static_cast<size_t>(cons_len));
            for (int i = 0; i < cons_len; ++i) {
                consensus.push_back(abpoa_code_to_nt(ab->abc->cons_base[0][i]));
            }
        }
    }

    abpoa_free(ab);
    abpoa_free_para(abpt);
    return consensus;
}

struct PostAssemblyTeDecision {
    bool pass = false;
    std::string qc = "FAIL_TE_CLASSIFICATION";
};

std::string anchor_qc_tag(const AnchorLockedReport& anchor_report) {
    if (!anchor_report.enabled) {
        return "ANCHOR_DISABLED";
    }
    if (anchor_report.fail_theta_uncertain) {
        return "ANCHOR_FAIL_THETA_UNCERTAIN";
    }
    if (!anchor_report.has_result) {
        return "ANCHOR_NO_CORE_RESULT";
    }
    return "ANCHOR_PASS";
}

bool has_prefix(const std::string& value, const char* prefix) {
    return value.rfind(prefix, 0) == 0;
}

PostAssemblyTeDecision evaluate_post_assembly_te_call(
    const ClusterTECall& te_call,
    const AssemblyCall& assembly,
    const ComponentCall& component,
    const PipelineConfig& config) {
    PostAssemblyTeDecision decision;

    if (te_call.te_name.empty()) {
        std::unordered_set<size_t> non_softclip_support;
        non_softclip_support.insert(
            component.split_sa_read_indices.begin(),
            component.split_sa_read_indices.end());
        non_softclip_support.insert(
            component.insertion_read_indices.begin(),
            component.insertion_read_indices.end());
        const int32_t non_softclip_reads = static_cast<int32_t>(non_softclip_support.size());
        const int32_t min_non_softclip_reads = std::max(2, config.te_min_fragments_for_vote);

        if (assembly.qc_pass &&
            assembly.identity_est >= std::clamp(config.assembly_min_identity_est, 0.0, 1.0) &&
            non_softclip_reads >= min_non_softclip_reads) {
            decision.pass = true;
            decision.qc = "PASS_INSERTION_TE_UNCERTAIN";
        } else {
            decision.qc = "FAIL_NO_TE_LABEL";
            return decision;
        }
    }

    if (!te_call.te_name.empty()) {
        const int32_t min_fragments = std::max(1, config.te_min_fragments_for_vote);
        if (te_call.fragment_count < min_fragments) {
            decision.qc = "FAIL_LOW_TE_FRAGMENTS";
            return decision;
        }

        const double vote_min = std::clamp(config.te_vote_fraction_min, 0.0, 1.0);
        const double identity_min = std::clamp(config.te_median_identity_min, 0.0, 1.0);
        const bool classic_pass =
            te_call.vote_fraction >= vote_min &&
            te_call.median_identity >= identity_min;
        if (classic_pass) {
            decision.pass = true;
            decision.qc = "PASS_CLASSIC";
        } else {
            const double rescue_vote_min = std::clamp(config.te_rescue_vote_fraction_min, 0.0, vote_min);
            const double rescue_identity_min = std::clamp(config.te_rescue_median_identity_min, 0.0, identity_min);
            const bool rescue_pass =
                assembly.qc_pass &&
                assembly.identity_est >= std::clamp(config.assembly_min_identity_est, 0.0, 1.0) &&
                te_call.vote_fraction >= rescue_vote_min &&
                te_call.median_identity >= rescue_identity_min;
            if (!rescue_pass) {
                decision.qc = "FAIL_LOW_TE_SUPPORT";
                return decision;
            }
            decision.pass = true;
            decision.qc = "PASS_ASM_RESCUE";
        }
    }

    const bool pure_softclip_component =
        component.split_sa_read_indices.empty() &&
        component.insertion_read_indices.empty() &&
        !component.soft_clip_read_indices.empty();
    if (!pure_softclip_component) {
        return decision;
    }

    const int32_t softclip_support_reads = static_cast<int32_t>(component.soft_clip_read_indices.size());
    if (softclip_support_reads < std::max(1, config.te_pure_softclip_min_reads)) {
        decision.pass = false;
        decision.qc = "FAIL_PURE_SOFTCLIP_LOW_READS";
        return decision;
    }
    if (te_call.fragment_count < std::max(1, config.te_pure_softclip_min_fragments)) {
        decision.pass = false;
        decision.qc = "FAIL_PURE_SOFTCLIP_LOW_FRAGMENTS";
        return decision;
    }

    const double pure_identity_min = std::clamp(config.te_pure_softclip_min_identity, 0.0, 1.0);
    if (std::max(te_call.median_identity, assembly.identity_est) < pure_identity_min) {
        decision.pass = false;
        decision.qc = "FAIL_PURE_SOFTCLIP_LOW_IDENTITY";
        return decision;
    }

    if (decision.qc == "PASS_CLASSIC") {
        decision.qc = "PASS_CLASSIC_PURE_SOFTCLIP";
    } else if (decision.qc == "PASS_ASM_RESCUE") {
        decision.qc = "PASS_ASM_RESCUE_PURE_SOFTCLIP";
    }
    return decision;
}

template <typename T>
class SafeQueue {
public:
    void push(T&& value) {
        {
            std::lock_guard<std::mutex> lock(mu_);
            queue_.push(std::move(value));
        }
        cv_.notify_one();
    }

    bool pop(T& out) {
        std::unique_lock<std::mutex> lock(mu_);
        cv_.wait(lock, [this]() { return finished_ || !queue_.empty(); });
        if (queue_.empty()) {
            return false;
        }
        out = std::move(queue_.front());
        queue_.pop();
        return true;
    }

    void close() {
        {
            std::lock_guard<std::mutex> lock(mu_);
            finished_ = true;
        }
        cv_.notify_all();
    }

private:
    std::queue<T> queue_;
    std::mutex mu_;
    std::condition_variable cv_;
    bool finished_ = false;
};

struct BinTask {
    std::vector<BamRecordPtr> records;
    int32_t tid = -1;
    int32_t bin_index = -1;
};

bool is_match_like(int op) {
    return op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF;
}

int find_first_non_hard_clip(const uint32_t* cigar, int32_t n_cigar) {
    for (int32_t i = 0; i < n_cigar; ++i) {
        if (bam_cigar_op(cigar[i]) != BAM_CHARD_CLIP) {
            return i;
        }
    }
    return -1;
}

int find_last_non_hard_clip(const uint32_t* cigar, int32_t n_cigar) {
    for (int32_t i = n_cigar - 1; i >= 0; --i) {
        if (bam_cigar_op(cigar[i]) != BAM_CHARD_CLIP) {
            return i;
        }
    }
    return -1;
}

int32_t compute_ref_end(const ReadView& read) {
    const uint32_t* cigar = read.cigar();
    const int32_t n_cigar = read.n_cigar();
    if (!cigar || n_cigar <= 0) {
        return read.pos();
    }
    int32_t ref_pos = read.pos();
    for (int32_t i = 0; i < n_cigar; ++i) {
        const int op = bam_cigar_op(cigar[i]);
        const int32_t len = static_cast<int32_t>(bam_cigar_oplen(cigar[i]));
        if ((bam_cigar_type(op) & 2) != 0) {
            ref_pos += len;
        }
    }
    return ref_pos;
}

int32_t weighted_median_position(std::vector<std::pair<int32_t, double>> pos_weights) {
    pos_weights.erase(
        std::remove_if(pos_weights.begin(), pos_weights.end(), [](const auto& row) {
            return row.second <= 0.0;
        }),
        pos_weights.end());
    if (pos_weights.empty()) {
        return -1;
    }

    std::sort(pos_weights.begin(), pos_weights.end(), [](const auto& a, const auto& b) {
        return a.first < b.first;
    });

    double total_weight = 0.0;
    for (const auto& row : pos_weights) {
        total_weight += row.second;
    }
    if (total_weight <= 0.0) {
        return pos_weights[pos_weights.size() / 2].first;
    }

    const double half_weight = 0.5 * total_weight;
    double acc_weight = 0.0;
    for (const auto& row : pos_weights) {
        acc_weight += row.second;
        if (acc_weight >= half_weight) {
            return row.first;
        }
    }
    return pos_weights.back().first;
}

std::vector<double> smooth_histogram(const std::vector<double>& hist) {
    std::vector<double> smooth(hist.size(), 0.0);
    for (size_t i = 0; i < hist.size(); ++i) {
        double value = hist[i];
        if (i >= 1) {
            value += 0.60 * hist[i - 1];
        }
        if (i + 1 < hist.size()) {
            value += 0.60 * hist[i + 1];
        }
        if (i >= 2) {
            value += 0.25 * hist[i - 2];
        }
        if (i + 2 < hist.size()) {
            value += 0.25 * hist[i + 2];
        }
        smooth[i] = value;
    }
    return smooth;
}

EvidenceBundle extract_evidence_points(
    const std::vector<BamRecordPtr>& bin_records,
    int32_t expected_tid) {
    EvidenceBundle bundle;
    bundle.read_summaries.resize(bin_records.size());

    for (size_t idx = 0; idx < bin_records.size(); ++idx) {
        if (!bin_records[idx]) {
            continue;
        }

        ReadView view(bin_records[idx].get());
        if (view.tid() != expected_tid) {
            continue;
        }

        ReadSignalSummary summary;
        summary.read_id = std::string(view.qname());
        summary.is_reverse = (view.flag() & BAM_FREVERSE) != 0;
        summary.has_sa_or_supp = view.has_sa_tag() || ((view.flag() & BAM_FSUPPLEMENTARY) != 0);

        const uint32_t* cigar = view.cigar();
        const int32_t n_cigar = view.n_cigar();
        if (!cigar || n_cigar <= 0) {
            if (summary.has_sa_or_supp) {
                summary.class_mask |= kCandidateSplitSaSupplementary;
            }
            bundle.read_summaries[idx] = std::move(summary);
            continue;
        }

        int32_t leading_soft = 0;
        int32_t trailing_soft = 0;
        const int first = find_first_non_hard_clip(cigar, n_cigar);
        const int last = find_last_non_hard_clip(cigar, n_cigar);
        if (first >= 0 && bam_cigar_op(cigar[first]) == BAM_CSOFT_CLIP) {
            leading_soft = static_cast<int32_t>(bam_cigar_oplen(cigar[first]));
        }
        if (last >= 0 && bam_cigar_op(cigar[last]) == BAM_CSOFT_CLIP) {
            trailing_soft = static_cast<int32_t>(bam_cigar_oplen(cigar[last]));
        }

        std::vector<EvidencePoint> local_points;
        local_points.reserve(8);

        int32_t ref_pos = view.pos();
        for (int32_t ci = 0; ci < n_cigar; ++ci) {
            const int op = bam_cigar_op(cigar[ci]);
            const int32_t len = static_cast<int32_t>(bam_cigar_oplen(cigar[ci]));

            if (op == BAM_CSOFT_CLIP) {
                summary.max_soft_clip = std::max(summary.max_soft_clip, len);
            } else if (op == BAM_CINS) {
                summary.max_ins = std::max(summary.max_ins, len);
            }

            if (op == BAM_CINS && len >= kLargeIndelEvidenceMin) {
                EvidencePoint point;
                point.read_index = idx;
                point.pos = ref_pos;
                point.weight = 0.80 + (static_cast<double>(std::min(len, 400)) / 200.0);
                point.kind = EvidenceKind::kIndel;
                point.signal_len = len;
                local_points.push_back(point);
            } else if (op == BAM_CDEL && len >= kLargeIndelEvidenceMin) {
                const double weight = 0.60 + (static_cast<double>(std::min(len, 400)) / 260.0);

                EvidencePoint left_point;
                left_point.read_index = idx;
                left_point.pos = ref_pos;
                left_point.weight = weight;
                left_point.kind = EvidenceKind::kIndel;
                left_point.signal_len = len;
                local_points.push_back(left_point);

                EvidencePoint right_point = left_point;
                right_point.pos = ref_pos + len;
                local_points.push_back(right_point);
            }

            if ((bam_cigar_type(op) & 2) != 0) {
                ref_pos += len;
            }
        }
        const int32_t ref_end = ref_pos;

        if (leading_soft >= kSoftClipSignalMin) {
            EvidencePoint point;
            point.read_index = idx;
            point.pos = view.pos();
            point.weight = 1.00 + (static_cast<double>(std::min(leading_soft, 500)) / 180.0);
            point.kind = EvidenceKind::kSoftClip;
            point.signal_len = leading_soft;
            local_points.push_back(point);
        }
        if (trailing_soft >= kSoftClipSignalMin) {
            EvidencePoint point;
            point.read_index = idx;
            point.pos = ref_end;
            point.weight = 1.00 + (static_cast<double>(std::min(trailing_soft, 500)) / 180.0);
            point.kind = EvidenceKind::kSoftClip;
            point.signal_len = trailing_soft;
            local_points.push_back(point);
        }

        if (summary.has_sa_or_supp) {
            EvidencePoint left_hint;
            left_hint.read_index = idx;
            left_hint.pos = view.pos();
            left_hint.weight = kSaHintWeight;
            left_hint.kind = EvidenceKind::kSAHint;
            local_points.push_back(left_hint);

            EvidencePoint right_hint = left_hint;
            right_hint.pos = ref_end;
            local_points.push_back(right_hint);
        }

        if (summary.has_sa_or_supp) {
            summary.class_mask |= kCandidateSplitSaSupplementary;
        }
        if (summary.max_soft_clip >= kSoftClipSignalMin) {
            summary.class_mask |= kCandidateSoftClip;
        }
        if (summary.max_ins >= kLongInsertionSignalMin) {
            summary.class_mask |= kCandidateLongInsertion;
        }

        for (auto& point : local_points) {
            point.class_mask = summary.class_mask;
            bundle.points.push_back(point);
        }
        bundle.read_summaries[idx] = std::move(summary);
    }

    return bundle;
}

std::vector<CandidateWindow> build_density_windows(
    const std::vector<EvidencePoint>& evidence,
    int32_t bin_start,
    int32_t bin_end) {
    if (evidence.empty()) {
        return {};
    }

    int32_t min_pos = evidence.front().pos;
    int32_t max_pos = evidence.front().pos;
    for (const auto& point : evidence) {
        min_pos = std::min(min_pos, point.pos);
        max_pos = std::max(max_pos, point.pos);
    }

    const int32_t hist_start = std::min(min_pos, bin_start) - kHistPadding;
    const int32_t hist_end = std::max(max_pos, bin_end) + kHistPadding;
    const int32_t num_bins = std::max(1, ((hist_end - hist_start) / kHistBinSize) + 1);

    std::vector<double> hist(static_cast<size_t>(num_bins), 0.0);
    for (const auto& point : evidence) {
        int32_t idx = (point.pos - hist_start) / kHistBinSize;
        idx = std::max(0, std::min(num_bins - 1, idx));
        hist[static_cast<size_t>(idx)] += point.weight;
    }
    const std::vector<double> smooth = smooth_histogram(hist);

    std::vector<DensityPeak> peaks;
    for (int32_t i = 0; i < num_bins; ++i) {
        const double center = smooth[static_cast<size_t>(i)];
        if (center < kPeakMinWeight) {
            continue;
        }
        const double left = (i > 0) ? smooth[static_cast<size_t>(i - 1)] : center;
        const double right = (i + 1 < num_bins) ? smooth[static_cast<size_t>(i + 1)] : center;
        if (center >= left && center >= right) {
            DensityPeak peak;
            peak.pos = hist_start + (i * kHistBinSize) + (kHistBinSize / 2);
            peak.weight = center;
            peaks.push_back(peak);
        }
    }

    if (peaks.empty()) {
        std::vector<std::pair<int32_t, double>> pos_weights;
        pos_weights.reserve(evidence.size());
        for (const auto& point : evidence) {
            pos_weights.push_back({point.pos, point.weight});
        }
        const int32_t fallback_center = weighted_median_position(std::move(pos_weights));
        if (fallback_center >= 0) {
            DensityPeak fallback_peak;
            fallback_peak.pos = fallback_center;
            fallback_peak.weight = 1.0;
            peaks.push_back(fallback_peak);
        }
    }
    if (peaks.empty()) {
        return {};
    }

    std::sort(peaks.begin(), peaks.end(), [](const DensityPeak& a, const DensityPeak& b) {
        return a.pos < b.pos;
    });

    std::vector<DensityPeak> merged_peaks;
    for (const auto& peak : peaks) {
        if (merged_peaks.empty() ||
            (peak.pos - merged_peaks.back().pos) > kPeakMergeDistance) {
            merged_peaks.push_back(peak);
            continue;
        }

        DensityPeak& back = merged_peaks.back();
        const double total_weight = back.weight + peak.weight;
        if (total_weight > 0.0) {
            const double weighted_center =
                ((back.weight * static_cast<double>(back.pos)) +
                 (peak.weight * static_cast<double>(peak.pos))) / total_weight;
            back.pos = static_cast<int32_t>(std::llround(weighted_center));
        }
        back.weight = total_weight;
    }

    std::vector<CandidateWindow> windows;
    windows.reserve(merged_peaks.size());

    for (const auto& peak : merged_peaks) {
        int32_t center_idx = (peak.pos - hist_start) / kHistBinSize;
        center_idx = std::max(0, std::min(num_bins - 1, center_idx));
        const double center_weight = smooth[static_cast<size_t>(center_idx)];
        const double expand_threshold = std::max(0.50, center_weight * kWindowDropRatio);

        int32_t left = center_idx;
        int32_t right = center_idx;

        while (left > 0 &&
               (center_idx - left) < kWindowMaxExpandBins &&
               smooth[static_cast<size_t>(left - 1)] >= expand_threshold) {
            --left;
        }
        while ((right + 1) < num_bins &&
               (right - center_idx) < kWindowMaxExpandBins &&
               smooth[static_cast<size_t>(right + 1)] >= expand_threshold) {
            ++right;
        }

        int32_t start = hist_start + (left * kHistBinSize);
        int32_t end = hist_start + ((right + 1) * kHistBinSize);
        if ((end - start) < kMinWindowSpan) {
            start = peak.pos - (kMinWindowSpan / 2);
            end = start + kMinWindowSpan;
        }

        start = std::max(bin_start, start);
        end = std::min(bin_end, end);
        if (end <= start) {
            continue;
        }

        CandidateWindow window;
        window.start = start;
        window.end = end;
        window.center = std::max(window.start, std::min(window.end - 1, peak.pos));
        window.peak_weight = std::max(peak.weight, center_weight);
        windows.push_back(window);
    }

    if (windows.empty()) {
        return windows;
    }

    std::sort(windows.begin(), windows.end(), [](const CandidateWindow& a, const CandidateWindow& b) {
        if (a.start != b.start) {
            return a.start < b.start;
        }
        return a.end < b.end;
    });

    std::vector<CandidateWindow> merged_windows;
    for (const auto& window : windows) {
        if (merged_windows.empty() ||
            window.start > (merged_windows.back().end + kWindowMergeGap)) {
            merged_windows.push_back(window);
            continue;
        }

        CandidateWindow& back = merged_windows.back();
        const double total_weight = back.peak_weight + window.peak_weight;
        if (total_weight > 0.0) {
            const double center =
                ((back.peak_weight * static_cast<double>(back.center)) +
                 (window.peak_weight * static_cast<double>(window.center))) / total_weight;
            back.center = static_cast<int32_t>(std::llround(center));
        }
        back.start = std::min(back.start, window.start);
        back.end = std::max(back.end, window.end);
        back.peak_weight = std::max(back.peak_weight, window.peak_weight);
    }

    return merged_windows;
}

}  // namespace

std::vector<ComponentCall> LinearBinComponentModule::build(
    const std::vector<BamRecordPtr>& bin_records,
    const std::string& chrom,
    int32_t tid,
    int32_t bin_start,
    int32_t bin_end) const {
    if (bin_records.empty()) {
        return {};
    }

    const EvidenceBundle evidence_bundle = extract_evidence_points(bin_records, tid);
    const auto windows = build_density_windows(evidence_bundle.points, bin_start, bin_end);
    if (windows.empty()) {
        return {};
    }

    std::vector<std::vector<const EvidencePoint*>> evidence_by_read(bin_records.size());
    for (const auto& point : evidence_bundle.points) {
        if (point.read_index < evidence_by_read.size()) {
            evidence_by_read[point.read_index].push_back(&point);
        }
    }

    std::vector<int32_t> assigned_window(bin_records.size(), -1);
    std::vector<int32_t> assigned_count(windows.size(), 0);
    std::vector<int32_t> ambiguous_count(windows.size(), 0);

    for (size_t read_idx = 0; read_idx < evidence_by_read.size(); ++read_idx) {
        const auto& read_points = evidence_by_read[read_idx];
        if (read_points.empty()) {
            continue;
        }

        int32_t best_window = -1;
        double best_score = 0.0;
        double second_score = 0.0;

        for (size_t wi = 0; wi < windows.size(); ++wi) {
            const auto& window = windows[wi];
            double score = 0.0;
            for (const EvidencePoint* point : read_points) {
                if (!point) {
                    continue;
                }
                const int32_t dist = std::abs(point->pos - window.center);
                const double proximity =
                    (point->pos >= window.start && point->pos <= window.end)
                    ? 1.0
                    : std::exp(-static_cast<double>(dist) / kReadAssignDecayBp);
                score += point->weight * proximity;
            }

            if (score > best_score) {
                second_score = best_score;
                best_score = score;
                best_window = static_cast<int32_t>(wi);
            } else if (score > second_score) {
                second_score = score;
            }
        }

        if (best_window < 0 || best_score < kReadAssignMinScore) {
            continue;
        }
        if (second_score >= best_score * kReadAmbiguousRatio) {
            ambiguous_count[static_cast<size_t>(best_window)] += 1;
            continue;
        }

        assigned_window[read_idx] = best_window;
        assigned_count[static_cast<size_t>(best_window)] += 1;
    }

    std::vector<ComponentCall> components;
    components.reserve(windows.size());

    for (size_t wi = 0; wi < windows.size(); ++wi) {
        if (assigned_count[wi] <= 0) {
            continue;
        }

        ComponentCall call;
        call.chrom = chrom;
        call.tid = tid;
        call.bin_start = windows[wi].start;
        call.bin_end = windows[wi].end;
        call.anchor_pos = windows[wi].center;
        call.peak_weight = windows[wi].peak_weight;

        const int32_t denom = assigned_count[wi] + ambiguous_count[wi];
        if (denom > 0) {
            call.ambiguous_frac = static_cast<double>(ambiguous_count[wi]) / static_cast<double>(denom);
        }

        std::vector<std::pair<int32_t, double>> anchor_points;
        for (size_t read_idx = 0; read_idx < assigned_window.size(); ++read_idx) {
            if (assigned_window[read_idx] != static_cast<int32_t>(wi)) {
                continue;
            }

            call.read_indices.push_back(read_idx);
            if (read_idx >= evidence_bundle.read_summaries.size()) {
                continue;
            }

            const auto& summary = evidence_bundle.read_summaries[read_idx];
            if (summary.max_soft_clip >= kSoftClipSignalMin) {
                call.soft_clip_read_indices.push_back(read_idx);
            }
            if (summary.has_sa_or_supp) {
                call.split_sa_read_indices.push_back(read_idx);
            }
            if (summary.max_ins >= kLongInsertionSignalMin) {
                call.insertion_read_indices.push_back(read_idx);
            }
        }

        for (const auto& point : evidence_bundle.points) {
            if (point.read_index >= assigned_window.size() ||
                assigned_window[point.read_index] != static_cast<int32_t>(wi)) {
                continue;
            }

            anchor_points.push_back({point.pos, point.weight});
            switch (point.kind) {
                case EvidenceKind::kSoftClip:
                    call.evidence_soft_clip_count += 1;
                    break;
                case EvidenceKind::kIndel:
                    call.evidence_indel_count += 1;
                    break;
                case EvidenceKind::kSAHint:
                    call.evidence_sa_hint_count += 1;
                    break;
            }

            BreakpointCandidate candidate;
            candidate.chrom = chrom;
            candidate.pos = point.pos;
            candidate.read_index = point.read_index;
            candidate.class_mask = point.class_mask;
            if (point.read_index < evidence_bundle.read_summaries.size()) {
                const auto& summary = evidence_bundle.read_summaries[point.read_index];
                candidate.read_id = summary.read_id;
                candidate.is_reverse = summary.is_reverse;
                candidate.anchor_len = summary.max_soft_clip;
            }
            if (point.kind == EvidenceKind::kSoftClip) {
                candidate.clip_len = point.signal_len;
            } else if (point.kind == EvidenceKind::kIndel) {
                candidate.ins_len = point.signal_len;
            }
            call.breakpoint_candidates.push_back(std::move(candidate));
        }

        if (!anchor_points.empty()) {
            const int32_t anchor = weighted_median_position(std::move(anchor_points));
            if (anchor >= 0) {
                call.anchor_pos = anchor;
            }
        }

        if (!call.read_indices.empty()) {
            components.push_back(std::move(call));
        }
    }

    return components;
}

Pipeline::Pipeline(PipelineConfig config, std::unique_ptr<BamStreamReader> bam_reader)
    : config_(std::move(config)),
      bam_reader_(std::move(bam_reader)),
      gate1_module_(),
      component_module_(),
      ins_fragment_module_(config_),
      te_classifier_module_(config_),
      anchor_locked_module_(config_) {}

PipelineResult Pipeline::run() const {
    if (!bam_reader_ || !bam_reader_->is_valid()) {
        throw std::runtime_error("BAM reader is not valid");
    }

    if (config_.enable_parallel) {
        return run_parallel();
    }
    return run_streaming();
}

PipelineResult Pipeline::run_streaming() const {
    PipelineResult result;
    StreamingState state;

    const auto progress_cb = [this](int64_t processed, int32_t tid) {
        std::cerr << "[Pipeline] processed=" << processed << " current_tid=" << tid << '\n';
        return true;
    };

    result.total_reads = bam_reader_->stream(
        [this, &state, &result](BamRecordPtr&& record) {
            consume_record(std::move(record), state, result);
        },
        progress_cb,
        config_.progress_interval);

    flush_current_bin(state, result);
    return result;
}

PipelineResult Pipeline::run_parallel() const {
    PipelineResult result;
    const int32_t worker_count = std::max(
        1,
        (config_.parallel_workers > 0) ? config_.parallel_workers : config_.bam_threads);
    SafeQueue<BinTask> bin_queue;
    std::vector<PipelineResult> worker_results(static_cast<size_t>(worker_count));
    std::vector<std::thread> workers;
    workers.reserve(static_cast<size_t>(worker_count));

    for (int32_t wi = 0; wi < worker_count; ++wi) {
        workers.emplace_back([this, wi, &bin_queue, &worker_results]() {
            BinTask task;
            while (bin_queue.pop(task)) {
                process_bin_records(
                    std::move(task.records),
                    task.tid,
                    task.bin_index,
                    worker_results[static_cast<size_t>(wi)]);
                task = BinTask();
            }
        });
    }

    int64_t gate1_passed = 0;
    int32_t current_tid = -1;
    int32_t current_bin_index = -1;
    std::vector<BamRecordPtr> current_bin_records;
    current_bin_records.reserve(config_.batch_size);

    const auto flush_bin_task = [&]() {
        if (current_bin_records.empty()) {
            return;
        }
        BinTask task;
        task.records = std::move(current_bin_records);
        task.tid = current_tid;
        task.bin_index = current_bin_index;
        bin_queue.push(std::move(task));
        current_bin_records.clear();
        current_bin_records.reserve(config_.batch_size);
    };

    const auto progress_cb = [this](int64_t processed, int32_t tid) {
        std::cerr << "[Pipeline] processed=" << processed << " current_tid=" << tid << '\n';
        return true;
    };

    const int64_t total_reads = bam_reader_->stream(
        [this, &gate1_passed, &current_tid, &current_bin_index, &current_bin_records, &flush_bin_task](
            BamRecordPtr&& record) {
            if (!record) {
                return;
            }

            ReadView view(record.get());
            if (!gate1_module_.pass_preliminary(view)) {
                return;
            }
            ++gate1_passed;

            const int32_t bin_index = (config_.bin_size > 0) ? (view.pos() / config_.bin_size) : 0;
            if (current_tid < 0) {
                current_tid = view.tid();
                current_bin_index = bin_index;
            }

            if (view.tid() != current_tid || bin_index != current_bin_index) {
                flush_bin_task();
                current_tid = view.tid();
                current_bin_index = bin_index;
            }
            current_bin_records.push_back(std::move(record));
        },
        progress_cb,
        config_.progress_interval);

    flush_bin_task();

    bin_queue.close();
    for (auto& worker : workers) {
        worker.join();
    }

    result.total_reads = total_reads;
    result.gate1_passed = gate1_passed;

    for (auto& worker_result : worker_results) {
        result.processed_bins += worker_result.processed_bins;
        result.built_components += worker_result.built_components;
        result.evidence_rows += worker_result.evidence_rows;
        result.assembled_calls += worker_result.assembled_calls;
        result.placeability_calls += worker_result.placeability_calls;
        result.genotype_calls += worker_result.genotype_calls;
        for (auto& call : worker_result.final_calls) {
            result.final_calls.push_back(std::move(call));
        }
    }

    std::sort(result.final_calls.begin(), result.final_calls.end(), [](const FinalCall& a, const FinalCall& b) {
        if (a.tid != b.tid) {
            return a.tid < b.tid;
        }
        if (a.pos != b.pos) {
            return a.pos < b.pos;
        }
        if (a.window_start != b.window_start) {
            return a.window_start < b.window_start;
        }
        if (a.window_end != b.window_end) {
            return a.window_end < b.window_end;
        }
        if (a.chrom != b.chrom) {
            return a.chrom < b.chrom;
        }
        return a.te_name < b.te_name;
    });

    return result;
}

void Pipeline::consume_record(
    BamRecordPtr&& record,
    StreamingState& state,
    PipelineResult& result) const {
    if (!record) {
        return;
    }

    ReadView view(record.get());
    if (!gate1_module_.pass_preliminary(view)) {
        return;
    }

    result.gate1_passed += 1;

    while (!state.active_window.empty()) {
        const auto& front = state.active_window.front();
        const bool different_tid = (view.tid() != front.tid);
        const bool beyond_window = (!different_tid && (view.pos() - front.pos > config_.window_size));
        if (!different_tid && !beyond_window) {
            break;
        }
        state.active_window.pop_front();
    }
    state.active_window.push_back({view.tid(), view.pos()});

    const int32_t bin_index = (config_.bin_size > 0) ? (view.pos() / config_.bin_size) : 0;

    if (state.current_tid < 0) {
        state.current_tid = view.tid();
        state.current_bin_index = bin_index;
    }

    if (view.tid() != state.current_tid || bin_index != state.current_bin_index) {
        flush_current_bin(state, result);
        state.current_tid = view.tid();
        state.current_bin_index = bin_index;
    }

    state.current_bin_records.push_back(std::move(record));
}

void Pipeline::flush_current_bin(
    StreamingState& state,
    PipelineResult& result) const {
    if (state.current_bin_records.empty()) {
        return;
    }

    process_bin_records(
        std::move(state.current_bin_records),
        state.current_tid,
        state.current_bin_index,
        result);

    state.current_bin_records.clear();
}

EvidenceFeatures Pipeline::collect_evidence(
    const ComponentCall& component,
    const std::vector<BamRecordPtr>& bin_records,
    const std::vector<InsertionFragment>& fragments,
    const ClusterTECall& te_call,
    const AnchorLockedReport& anchor_report) const {
    EvidenceFeatures features;
    features.tid = component.tid;
    features.pos = component.anchor_pos;
    features.global_cov_reads = static_cast<int32_t>(bin_records.size());
    features.local_cov_reads = static_cast<int32_t>(component.read_indices.size());
    features.evidence_point_count = static_cast<int32_t>(component.breakpoint_candidates.size());
    features.te_vote_fraction = te_call.vote_fraction;
    features.te_median_identity = te_call.median_identity;
    features.ambiguous_frac = std::clamp(component.ambiguous_frac, 0.0, 1.0);
    features.anchor_consensus_ok = anchor_report.has_result && !anchor_report.fail_theta_uncertain;

    std::unordered_set<size_t> alt_support_set;
    std::vector<int32_t> breakpoint_positions;
    breakpoint_positions.reserve(component.breakpoint_candidates.size());
    for (const auto& bp : component.breakpoint_candidates) {
        alt_support_set.insert(bp.read_index);
        if (bp.pos >= 0) {
            breakpoint_positions.push_back(bp.pos);
        }
    }
    if (alt_support_set.empty()) {
        alt_support_set.insert(
            component.split_sa_read_indices.begin(),
            component.split_sa_read_indices.end());
        alt_support_set.insert(
            component.insertion_read_indices.begin(),
            component.insertion_read_indices.end());
        alt_support_set.insert(
            component.soft_clip_read_indices.begin(),
            component.soft_clip_read_indices.end());
    }
    if (alt_support_set.empty()) {
        alt_support_set.insert(component.read_indices.begin(), component.read_indices.end());
    }
    if (breakpoint_positions.empty()) {
        breakpoint_positions.reserve(alt_support_set.size());
        for (size_t idx : alt_support_set) {
            if (idx >= bin_records.size() || !bin_records[idx]) {
                continue;
            }
            ReadView view(bin_records[idx].get());
            if (view.pos() >= 0) {
                breakpoint_positions.push_back(view.pos());
            }
        }
    }
    features.breakpoint_mad = mad_from_positions(breakpoint_positions);

    const int32_t anchor_margin = 100;
    std::unordered_set<size_t> depth_set;
    depth_set.reserve(bin_records.size());
    const int32_t anchor_start = component.anchor_pos - anchor_margin;
    const int32_t anchor_end = component.anchor_pos + anchor_margin;
    for (size_t idx = 0; idx < bin_records.size(); ++idx) {
        if (!bin_records[idx]) {
            continue;
        }
        ReadView view(bin_records[idx].get());
        if (view.tid() != component.tid) {
            continue;
        }
        const int32_t read_start = view.pos();
        const int32_t read_end = compute_ref_end(view);
        if (read_end < anchor_start || read_start > anchor_end) {
            continue;
        }
        depth_set.insert(idx);
    }

    std::unordered_set<size_t> split_sa_set(
        component.split_sa_read_indices.begin(),
        component.split_sa_read_indices.end());
    std::unordered_set<size_t> insertion_set(
        component.insertion_read_indices.begin(),
        component.insertion_read_indices.end());
    std::unordered_set<size_t> softclip_set(
        component.soft_clip_read_indices.begin(),
        component.soft_clip_read_indices.end());

    int32_t alt_overlap_support = 0;
    for (size_t idx : depth_set) {
        if (alt_support_set.find(idx) != alt_support_set.end()) {
            ++alt_overlap_support;
        }
    }

    int32_t depth_reads = static_cast<int32_t>(depth_set.size());
    if (depth_reads <= 0) {
        depth_reads = std::max(
            static_cast<int32_t>(alt_support_set.size()),
            features.local_cov_reads);
    }
    int32_t alt_support_reads = alt_overlap_support;
    if (alt_support_reads <= 0) {
        alt_support_reads = static_cast<int32_t>(alt_support_set.size());
    }
    alt_support_reads = std::max(0, std::min(alt_support_reads, depth_reads));
    const int32_t ref_support_reads = std::max(0, depth_reads - alt_support_reads);

    features.depth_reads = depth_reads;
    features.alt_support_reads = alt_support_reads;
    features.ref_support_reads = ref_support_reads;

    for (size_t read_idx : alt_support_set) {
        double weight = 0.0;
        if (split_sa_set.find(read_idx) != split_sa_set.end()) {
            weight += 1.0;
        }
        if (insertion_set.find(read_idx) != insertion_set.end()) {
            weight += 0.9;
        }
        if (softclip_set.find(read_idx) != softclip_set.end()) {
            weight += 0.6;
        }
        if (weight <= 0.0) {
            weight = 0.5;
        }
        features.effective_alt_support += weight;
    }

    int32_t softclip_total = 0;
    int32_t softclip_low_complex = 0;
    for (const auto& frag : fragments) {
        if (!is_softclip_source(frag.source)) {
            continue;
        }
        ++softclip_total;
        if (is_low_complexity_softclip(frag, config_)) {
            ++softclip_low_complex;
        }
    }
    if (softclip_total > 0) {
        features.low_complex_softclip_frac =
            static_cast<double>(softclip_low_complex) / static_cast<double>(softclip_total);
    }

    const double alpha = std::clamp(config_.evidence_min_support_alpha, 0.0, 1.0);
    const double lambda = std::clamp(config_.evidence_min_support_lambda, 0.0, 1.0);
    const double expected_cov =
        ((1.0 - lambda) * static_cast<double>(std::max(0, features.global_cov_reads))) +
        (lambda * static_cast<double>(std::max(0, features.local_cov_reads)));
    const int32_t min_support = std::max(
        1,
        static_cast<int32_t>(std::ceil(alpha * expected_cov)));

    features.min_support_required = min_support;
    features.pass_min_support = features.alt_support_reads >= min_support;
    features.pass_breakpoint_mad = features.breakpoint_mad <= std::max(0.0, config_.evidence_breakpoint_mad_max);
    features.pass_low_complexity = features.low_complex_softclip_frac <=
        std::clamp(config_.evidence_low_complex_softclip_frac_max, 0.0, 1.0);
    // Stage-A weak gate: do not require TE label; allow insertion evidence to carry candidates.
    const int32_t weak_min_fragments = std::max(1, config_.te_min_fragments_for_vote / 2);
    const bool te_fragment_support = te_call.fragment_count >= weak_min_fragments;
    const bool insertion_fragment_support =
        static_cast<int32_t>(fragments.size()) >= weak_min_fragments ||
        features.alt_support_reads >= weak_min_fragments;
    features.pass_te_consistency = te_fragment_support || insertion_fragment_support;
    features.hard_filtered =
        !(features.pass_min_support &&
          features.pass_breakpoint_mad &&
          features.pass_low_complexity &&
          features.pass_te_consistency);

    return features;
}

AssemblyCall Pipeline::assemble_component(
    const ComponentCall& component,
    const std::vector<InsertionFragment>& fragments,
    const std::vector<FragmentTEHit>& hits,
    const ClusterTECall& te_call) const {
    AssemblyCall call;
    call.tid = component.tid;
    call.pos = component.anchor_pos;
    call.assembly_mode = "NONE";
    call.qc.clear();

    (void)hits;
    (void)te_call;

    std::vector<std::string> candidates;
    candidates.reserve(fragments.size());
    for (const auto& frag : fragments) {
        if (frag.sequence.empty() || frag.length < config_.assembly_min_fragment_len) {
            continue;
        }

        candidates.push_back(upper_acgt(frag.sequence));
    }

    call.input_fragments = static_cast<int32_t>(candidates.size());
    if (candidates.empty()) {
        call.qc_pass = false;
        call.qc = "NO_CANDIDATE_FRAGMENTS";
        return call;
    }
    call.used_fragments = static_cast<int32_t>(candidates.size());
    if (call.used_fragments < config_.assembly_poa_min_reads) {
        call.qc_pass = false;
        call.qc = "INSUFFICIENT_POA_READS";
        return call;
    }
    std::sort(candidates.begin(), candidates.end(), [](const std::string& a, const std::string& b) {
        if (a.size() != b.size()) {
            return a.size() > b.size();
        }
        return a < b;
    });
    if (static_cast<int32_t>(candidates.size()) > config_.assembly_poa_max_reads) {
        candidates.resize(static_cast<size_t>(config_.assembly_poa_max_reads));
    }
    call.used_fragments = static_cast<int32_t>(candidates.size());
    call.consensus = abpoa_consensus(candidates);
    call.assembly_mode = "POA";

    call.consensus = upper_acgt(call.consensus);
    call.consensus_len = static_cast<int32_t>(call.consensus.size());
    call.identity_est = estimate_mean_identity_kmer(call.consensus, candidates, config_.assembly_kmer_size);

    if (call.consensus.empty()) {
        call.qc_pass = false;
        call.qc = "EMPTY_CONSENSUS";
    } else if (call.consensus_len < config_.assembly_min_consensus_len) {
        call.qc_pass = false;
        call.qc = "CONSENSUS_TOO_SHORT";
    } else if (call.used_fragments >= 2 && call.identity_est < config_.assembly_min_identity_est) {
        call.qc_pass = false;
        call.qc = "LOW_IDENTITY";
    } else {
        call.qc_pass = true;
        call.qc = "PASS";
    }

    return call;
}

PlaceabilityReport Pipeline::score_placeability(
    const AssemblyCall& assembly,
    const EvidenceFeatures& evidence) const {
    PlaceabilityReport report;
    report.tid = assembly.tid;
    report.pos = assembly.pos;
    report.depth_reads = evidence.depth_reads;
    report.support_reads = evidence.alt_support_reads;
    report.ref_support_reads = evidence.ref_support_reads;
    report.min_support_required = evidence.min_support_required;
    report.breakpoint_mad = evidence.breakpoint_mad;
    report.low_complex_softclip_frac = evidence.low_complex_softclip_frac;
    report.hard_filtered = evidence.hard_filtered;

    if (evidence.hard_filtered) {
        report.delta_score = 0.0;
        report.evidence_p = 0.0;
        report.evidence_q = 0.0;
        report.tier = 3;
        return report;
    }

    const double anchor_bonus = evidence.anchor_consensus_ok ? 1.0 : 0.0;
    const double mad_term = std::exp(-std::max(0.0, evidence.breakpoint_mad) / 25.0);
    const double z =
        (1.8 * std::log1p(std::max(0.0, evidence.effective_alt_support))) +
        (1.2 * mad_term) +
        (1.2 * evidence.te_vote_fraction) +
        (1.0 * evidence.te_median_identity) +
        (0.7 * anchor_bonus) -
        (1.2 * evidence.ambiguous_frac) -
        (0.8 * evidence.low_complex_softclip_frac);

    const double p = std::clamp(sigmoid(z - config_.evidence_logit_bias), 0.0, 1.0);
    const double q = -10.0 * std::log10(std::max(1e-9, 1.0 - p));
    report.evidence_p = p;
    report.evidence_q = q;
    report.delta_score = q;

    const double tier1 = std::clamp(config_.evidence_tier1_prob, 0.0, 1.0);
    const double tier2 = std::clamp(config_.evidence_tier2_prob, 0.0, tier1);
    if (p >= tier1) {
        report.tier = 1;
    } else if (p >= tier2) {
        report.tier = 2;
    } else {
        report.tier = 3;
    }
    return report;
}

GenotypeCall Pipeline::genotype_call(
    const AssemblyCall& assembly,
    const PlaceabilityReport& placeability) const {
    GenotypeCall call;
    call.tid = assembly.tid;
    call.pos = assembly.pos;

    int32_t depth = std::max(0, placeability.depth_reads);
    int32_t alt = std::max(0, placeability.support_reads);
    int32_t ref = std::max(0, placeability.ref_support_reads);

    if (depth <= 0) {
        depth = alt + ref;
    }
    if (depth <= 0) {
        call.genotype = "./.";
        call.af = 0.0;
        call.gq = 0;
        return call;
    }
    if ((alt + ref) <= 0) {
        alt = std::min(alt, depth);
        ref = std::max(0, depth - alt);
    } else if ((alt + ref) != depth) {
        depth = alt + ref;
    }

    alt = std::max(0, std::min(alt, depth));
    ref = std::max(0, depth - alt);
    call.af = std::clamp(
        static_cast<double>(alt) / static_cast<double>(std::max(1, depth)),
        0.0,
        1.0);

    if (depth < std::max(1, config_.genotype_min_depth)) {
        call.genotype = "./.";
        call.gq = 0;
        return call;
    }

    const double error_rate = std::clamp(config_.genotype_error_rate, 1e-4, 0.25);
    const double ll_00 = log_binomial_likelihood(alt, ref, error_rate);
    const double ll_01 = log_binomial_likelihood(alt, ref, 0.5);
    const double ll_11 = log_binomial_likelihood(alt, ref, 1.0 - error_rate);

    std::string best_gt = "0/0";
    double best_ll = ll_00;
    double second_ll = -std::numeric_limits<double>::infinity();

    auto update_best = [&](const char* gt, double ll) {
        if (ll > best_ll) {
            second_ll = best_ll;
            best_ll = ll;
            best_gt = gt;
        } else if (ll > second_ll) {
            second_ll = ll;
        }
    };

    update_best("0/1", ll_01);
    update_best("1/1", ll_11);
    call.genotype = best_gt;

    const double delta_ll = std::max(0.0, best_ll - second_ll);
    const double gq = 4.3429448190325175 * delta_ll;  // -10*log10(exp(-delta_ll))
    call.gq = std::max(0, std::min(99, static_cast<int32_t>(std::lround(gq))));

    return call;
}

void Pipeline::process_bin_records(
    std::vector<BamRecordPtr>&& bin_records,
    int32_t tid,
    int32_t bin_index,
    PipelineResult& result) const {
    if (bin_records.empty()) {
        return;
    }

    result.processed_bins += 1;

    const int32_t bin_start = bin_index * config_.bin_size;
    const int32_t bin_end = bin_start + config_.bin_size;
    const std::string chrom = bam_reader_->chromosome_name(tid);

    auto components = component_module_.build(bin_records, chrom, tid, bin_start, bin_end);
    result.built_components += static_cast<int64_t>(components.size());

    for (const auto& component : components) {
        std::vector<InsertionFragment> fragments = ins_fragment_module_.extract(component, bin_records);
        std::vector<FragmentTEHit> hits;
        ClusterTECall te_call;

        if (te_classifier_module_.is_enabled()) {
            hits = te_classifier_module_.classify(fragments);
            te_call = te_classifier_module_.vote_cluster(hits);
        }

        AnchorLockedReport anchor_report;
        if (anchor_locked_module_.is_enabled()) {
            anchor_report = anchor_locked_module_.resolve(component, fragments, hits, te_call);
        }

        auto evidence = collect_evidence(component, bin_records, fragments, te_call, anchor_report);
        result.evidence_rows += static_cast<int64_t>(evidence.evidence_point_count);
        if (evidence.hard_filtered) {
            continue;
        }

        AssemblyCall assembly = assemble_component(component, fragments, hits, te_call);
        result.assembled_calls += 1;
        if (!assembly.qc_pass) {
            continue;
        }

        const PostAssemblyTeDecision te_decision =
            evaluate_post_assembly_te_call(te_call, assembly, component, config_);
        if (!te_decision.pass) {
            continue;
        }

        PlaceabilityReport placeability = score_placeability(assembly, evidence);
        result.placeability_calls += 1;

        GenotypeCall genotype = genotype_call(assembly, placeability);
        result.genotype_calls += 1;

        FinalCall call;
        call.chrom = component.chrom;
        call.tid = assembly.tid;
        call.pos = assembly.pos;
        call.window_start = component.bin_start;
        call.window_end = component.bin_end;
        call.te_name = te_call.te_name.empty() ? "UNK" : te_call.te_name;
        call.te_vote_fraction = te_call.vote_fraction;
        call.te_median_identity = te_call.median_identity;
        call.te_fragment_count = te_call.fragment_count;
        call.te_top1_name = te_call.top1_te_name;
        call.te_top2_name = te_call.top2_te_name;
        call.te_posterior_top1 = te_call.posterior_top1;
        call.te_posterior_top2 = te_call.posterior_top2;
        call.te_posterior_margin = te_call.posterior_margin;
        call.te_theta = theta_tag(anchor_report.theta0);
        call.te_mad_fwd = anchor_report.mad_fwd;
        call.te_mad_rev = anchor_report.mad_rev;
        call.te_breakpoint_core = anchor_report.te_breakpoint_core;
        call.te_breakpoint_window_start = anchor_report.te_breakpoint_window_start;
        call.te_breakpoint_window_end = anchor_report.te_breakpoint_window_end;
        call.te_core_candidates = anchor_report.core_candidate_count;
        call.te_core_set = anchor_report.core_set_count;
        call.te_split_sa_core_frac = anchor_report.split_sa_core_frac;
        call.te_ref_junc_pos_min = anchor_report.ref_junc_pos_min;
        call.te_ref_junc_pos_max = anchor_report.ref_junc_pos_max;
        call.te_qc = te_decision.qc + "|" + anchor_qc_tag(anchor_report);
        call.te_status = has_prefix(te_decision.qc, "PASS_INSERTION_TE_UNCERTAIN")
            ? "TE_UNCERTAIN"
            : "TE_CERTAIN";

        call.tier = placeability.tier;
        call.support_reads = placeability.support_reads;
        call.genotype = genotype.genotype;
        call.af = genotype.af;
        call.gq = genotype.gq;
        call.asm_mode = assembly.assembly_mode;
        call.asm_input_fragments = assembly.input_fragments;
        call.asm_used_fragments = assembly.used_fragments;
        call.asm_consensus_len = assembly.consensus_len;
        call.asm_identity_est = assembly.identity_est;
        call.asm_qc = assembly.qc;
        result.final_calls.push_back(std::move(call));
    }
}

std::unique_ptr<Pipeline> build_default_pipeline(const PipelineConfig& config) {
    auto reader = make_bam_reader(config.bam_path, config.bam_threads);
    return std::make_unique<Pipeline>(config, std::move(reader));
}

}  // namespace placer
