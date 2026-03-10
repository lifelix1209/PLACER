#include "denovo.h"

#include "bam_io.h"
#include "pipeline.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <memory>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <htslib/sam.h>

namespace placer {
namespace {

struct ClipInfo {
    int32_t leading = 0;
    int32_t trailing = 0;
    int32_t ref_end = 0;
};

struct InsOp {
    int32_t start = 0;
    int32_t len = 0;
    int32_t ref_pos = -1;
};

enum class EvidenceMatchClass : uint8_t {
    kNone = 0,
    kAmbiguous = 1,
    kFamily = 2,
    kExact = 3
};

struct ReadEvidenceDecision {
    EvidenceMatchClass match_class = EvidenceMatchClass::kNone;
    ParentReadEvidence evidence;
};

struct CandidateScanWindow {
    int32_t fetch_start = -1;
    int32_t fetch_end = -1;
    int32_t match_start = -1;
    int32_t match_end = -1;
};

bool is_match_like(int op) {
    return op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF;
}

std::string upper_copy(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::toupper(ch));
    });
    return value;
}

std::string sanitize_token(std::string value) {
    for (char& ch : value) {
        const unsigned char uch = static_cast<unsigned char>(ch);
        if (!(std::isalnum(uch) || ch == '.' || ch == '_' || ch == '-')) {
            ch = '_';
        }
    }
    return value;
}

std::string infer_te_family(const std::string& te_name) {
    const std::string upper = upper_copy(te_name);
    if (upper.empty() || upper == "NA") {
        return "NA";
    }
    if (upper.rfind("ALU", 0) == 0) {
        return "ALU";
    }
    if (upper.rfind("SVA", 0) == 0) {
        return "SVA";
    }
    if (upper.rfind("L1", 0) == 0 || upper.rfind("LINE1", 0) == 0 || upper.rfind("LINE_1", 0) == 0) {
        return "LINE1";
    }
    if (upper.rfind("ERV", 0) == 0 || upper.find("HERV") != std::string::npos) {
        return "ERV";
    }
    const std::string delimiters = "|#:/ _-";
    const size_t cut = upper.find_first_of(delimiters);
    if (cut == std::string::npos || cut == 0) {
        return upper;
    }
    return upper.substr(0, cut);
}

std::string evidence_type_for_source(InsertionFragmentSource source) {
    switch (source) {
        case InsertionFragmentSource::kClipRefLeft:
        case InsertionFragmentSource::kClipRefRight:
            return "soft_clip";
        case InsertionFragmentSource::kCigarInsertion:
            return "long_insertion";
        case InsertionFragmentSource::kSplitSa:
            return "split_sa";
        default:
            break;
    }
    return "unknown";
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
    if (cigar == nullptr || n_cigar <= 0) {
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

ClipInfo analyze_clip_info(const ReadView& read) {
    ClipInfo info;
    const uint32_t* cigar = read.cigar();
    const int32_t n_cigar = read.n_cigar();
    if (cigar == nullptr || n_cigar <= 0) {
        info.ref_end = read.pos();
        return info;
    }

    const int first = find_first_non_hard_clip(cigar, n_cigar);
    const int last = find_last_non_hard_clip(cigar, n_cigar);
    if (first >= 0 && bam_cigar_op(cigar[first]) == BAM_CSOFT_CLIP) {
        info.leading = static_cast<int32_t>(bam_cigar_oplen(cigar[first]));
    }
    if (last >= 0 && bam_cigar_op(cigar[last]) == BAM_CSOFT_CLIP) {
        info.trailing = static_cast<int32_t>(bam_cigar_oplen(cigar[last]));
    }
    info.ref_end = compute_ref_end(read);
    return info;
}

std::vector<InsOp> find_long_insertions(const ReadView& read, int32_t min_long_ins) {
    std::vector<InsOp> out;
    const uint32_t* cigar = read.cigar();
    const int32_t n_cigar = read.n_cigar();
    if (cigar == nullptr || n_cigar <= 0) {
        return out;
    }

    int32_t qpos = 0;
    int32_t rpos = read.pos();
    for (int32_t i = 0; i < n_cigar; ++i) {
        const int op = bam_cigar_op(cigar[i]);
        const int32_t len = static_cast<int32_t>(bam_cigar_oplen(cigar[i]));

        if (is_match_like(op)) {
            qpos += len;
            rpos += len;
            continue;
        }

        if (op == BAM_CINS) {
            if (len >= min_long_ins) {
                InsOp ins;
                ins.start = qpos;
                ins.len = len;
                ins.ref_pos = rpos;
                out.push_back(ins);
            }
            qpos += len;
            continue;
        }

        if ((bam_cigar_type(op) & 1) != 0) {
            qpos += len;
        }
        if ((bam_cigar_type(op) & 2) != 0) {
            rpos += len;
        }
    }

    return out;
}

CandidateScanWindow build_scan_window(const DenovoChildCandidate& candidate, const DenovoConfig& config) {
    const int32_t event_start = std::max(0, candidate.event_start);
    const int32_t event_end = std::max(event_start, candidate.event_end);
    const int32_t span = std::max(1, event_end - event_start + 1);
    const int32_t adaptive_padding = std::clamp(
        std::max(config.default_match_window, span / 2),
        config.default_match_window,
        config.max_match_window);

    CandidateScanWindow window;
    window.match_start = std::max(0, event_start - adaptive_padding);
    window.match_end = event_end + adaptive_padding;
    window.fetch_start = std::max(0, event_start - config.fetch_window);
    window.fetch_end = event_end + config.fetch_window;
    return window;
}

bool position_in_window(int32_t pos, const CandidateScanWindow& window) {
    return pos >= window.match_start && pos <= window.match_end;
}

std::optional<int32_t> nearest_breakpoint_in_window(
    int32_t a,
    int32_t b,
    const CandidateScanWindow& window) {
    const bool a_ok = position_in_window(a, window);
    const bool b_ok = position_in_window(b, window);
    if (!a_ok && !b_ok) {
        return std::nullopt;
    }
    if (a_ok && !b_ok) {
        return a;
    }
    if (!a_ok && b_ok) {
        return b;
    }
    const int32_t mid = (window.match_start + window.match_end) / 2;
    return (std::abs(a - mid) <= std::abs(b - mid)) ? a : b;
}

bool parse_sa_field(const std::string& field, std::string& chrom, int32_t& pos_1based) {
    chrom.clear();
    pos_1based = -1;

    std::stringstream ss(field);
    std::string token;
    if (!std::getline(ss, chrom, ',')) {
        return false;
    }
    if (!std::getline(ss, token, ',')) {
        return false;
    }

    try {
        pos_1based = std::stoi(token);
        return !chrom.empty() && pos_1based > 0;
    } catch (...) {
        return false;
    }
}

bool has_sa_breakpoint_in_window(const ReadView& read, const std::string& chrom, const CandidateScanWindow& window) {
    std::string sa;
    if (!read.get_string_tag("SA", sa) || sa.empty()) {
        return false;
    }

    std::stringstream ss(sa);
    std::string field;
    while (std::getline(ss, field, ';')) {
        if (field.empty()) {
            continue;
        }
        std::string sa_chrom;
        int32_t sa_pos_1based = -1;
        if (!parse_sa_field(field, sa_chrom, sa_pos_1based)) {
            continue;
        }
        const int32_t sa_pos_0based = std::max(0, sa_pos_1based - 1);
        if (sa_chrom == chrom && position_in_window(sa_pos_0based, window)) {
            return true;
        }
    }

    return false;
}

InsertionFragment make_fragment(
    const DenovoChildCandidate& child,
    const std::string& read_name,
    const std::string& tag,
    int32_t breakpoint_pos,
    InsertionFragmentSource source,
    int32_t start,
    const std::string& sequence) {
    InsertionFragment fragment;
    fragment.fragment_id = child.chrom + ":" + std::to_string(child.pos) + "|" +
        sanitize_token(read_name) + "|" + tag + "|" + std::to_string(breakpoint_pos);
    fragment.chrom = child.chrom;
    fragment.anchor_pos = breakpoint_pos;
    fragment.read_id = read_name;
    fragment.start = start;
    fragment.length = static_cast<int32_t>(sequence.size());
    fragment.read_len = static_cast<int32_t>(sequence.size());
    fragment.source = source;
    fragment.sequence = sequence;
    return fragment;
}

ReadEvidenceDecision choose_better_decision(
    const ReadEvidenceDecision& lhs,
    const ReadEvidenceDecision& rhs) {
    if (static_cast<int>(rhs.match_class) > static_cast<int>(lhs.match_class)) {
        return rhs;
    }
    if (static_cast<int>(rhs.match_class) < static_cast<int>(lhs.match_class)) {
        return lhs;
    }
    if (rhs.evidence.fragment_len > lhs.evidence.fragment_len) {
        return rhs;
    }
    return lhs;
}

}  // namespace

struct ParentPoolScanner::Impl {
    explicit Impl(const DenovoConfig& config)
        : pipeline_config(make_pipeline_config(config)),
          classifier(pipeline_config),
          split_sa_module(pipeline_config),
          min_softclip_len(std::max(1, config.min_softclip_len)),
          min_insertion_len(std::max(1, config.min_insertion_len)),
          parent_mapq_min(std::max(0, config.parent_mapq_min)) {}

    static PipelineConfig make_pipeline_config(const DenovoConfig& config) {
        PipelineConfig pipeline_config;
        pipeline_config.te_fasta_path = config.te_fasta_path;
        pipeline_config.ins_fragment_hits_tsv_path.clear();
        pipeline_config.ins_fragments_fasta_path.clear();
        pipeline_config.te_kmer_size = 13;
        pipeline_config.te_kmer_sizes_csv = "9,11,13";
        pipeline_config.te_median_identity_min = 0.20;
        pipeline_config.te_low_kmer_rescue_enable = true;
        pipeline_config.te_low_kmer_rescue_topn = 3;
        pipeline_config.te_low_kmer_rescue_min_frag_len = 40;
        pipeline_config.min_soft_clip_for_seq_extract = std::max(1, config.min_softclip_len);
        pipeline_config.min_long_ins_for_seq_extract = std::max(1, config.min_insertion_len);
        pipeline_config.min_sa_aln_len_for_seq_extract =
            std::max(1, std::min(config.min_softclip_len, config.min_insertion_len));
        pipeline_config.max_sa_per_read = 3;
        pipeline_config.short_ins_enable = false;
        return pipeline_config;
    }

    std::vector<ReadEvidenceDecision> classify_fragments(
        const DenovoChildCandidate& child,
        const std::string& bam_path,
        const std::vector<InsertionFragment>& fragments) const {
        std::vector<ReadEvidenceDecision> decisions;
        if (fragments.empty()) {
            return decisions;
        }

        std::vector<FragmentTEHit> hits;
        if (classifier.is_enabled()) {
            hits = classifier.classify(fragments);
        }

        decisions.reserve(fragments.size());
        for (size_t i = 0; i < fragments.size(); ++i) {
            const auto& fragment = fragments[i];
            ReadEvidenceDecision decision;
            decision.evidence.parent_bam_path = bam_path;
            decision.evidence.read_name =
                fragment.read_id.empty() ? sanitize_token(fragment.fragment_id) : fragment.read_id;
            decision.evidence.evidence_type = evidence_type_for_source(fragment.source);
            decision.evidence.breakpoint_pos =
                fragment.ref_junc_pos >= 0 ? fragment.ref_junc_pos : fragment.anchor_pos;
            decision.evidence.fragment_len = static_cast<int32_t>(fragment.sequence.size());

            if (i < hits.size() && !hits[i].te_name.empty()) {
                const std::string child_te_upper = upper_copy(child.te_name);
                const std::string hit_te_upper = upper_copy(hits[i].te_name);
                const std::string child_family = infer_te_family(child.te_name);
                const std::string hit_family = infer_te_family(hits[i].te_name);

                decision.evidence.matched_te = hits[i].te_name;
                decision.evidence.matched_family = hit_family;
                if (hit_te_upper == child_te_upper) {
                    decision.match_class = EvidenceMatchClass::kExact;
                    decision.evidence.veto_reason = "EXACT_TE_MATCH";
                } else if (child_family != "NA" && hit_family != "NA" && child_family == hit_family) {
                    decision.match_class = EvidenceMatchClass::kFamily;
                    decision.evidence.veto_reason = "FAMILY_MATCH";
                } else {
                    decision.match_class = EvidenceMatchClass::kAmbiguous;
                    decision.evidence.veto_reason = "NON_MATCHING_TE_FRAGMENT";
                }
            } else {
                decision.match_class = EvidenceMatchClass::kAmbiguous;
                decision.evidence.veto_reason = "UNCLASSIFIED_INSERTION_FRAGMENT";
            }

            decisions.push_back(std::move(decision));
        }

        return decisions;
    }

    ReadEvidenceDecision classify_read(
        const DenovoChildCandidate& child,
        const std::string& bam_path,
        const bam1_t* record,
        const CandidateScanWindow& window,
        bool family_match_veto) const {
        (void)family_match_veto;
        ReadEvidenceDecision best;
        if (record == nullptr) {
            return best;
        }

        ReadView read(record);
        const uint16_t flag = read.flag();
        if ((flag & BAM_FUNMAP) != 0 || (flag & BAM_FSECONDARY) != 0) {
            return best;
        }
        if (read.mapq() < parent_mapq_min) {
            return best;
        }

        const std::string read_name = sanitize_token(std::string(read.qname()));
        const ClipInfo clip = analyze_clip_info(read);
        std::vector<InsertionFragment> fragments;
        fragments.reserve(4);

        if (clip.leading >= min_softclip_len && position_in_window(read.pos(), window)) {
            const std::string sequence = read.decode_subsequence(0, clip.leading);
            if (!sequence.empty()) {
                fragments.push_back(make_fragment(
                    child,
                    read_name,
                    "softclip_left",
                    read.pos(),
                    InsertionFragmentSource::kClipRefLeft,
                    0,
                    sequence));
            }
        }

        if (clip.trailing >= min_softclip_len && position_in_window(clip.ref_end, window)) {
            const int32_t start = std::max(0, read.seq_len() - clip.trailing);
            const std::string sequence = read.decode_subsequence(start, clip.trailing);
            if (!sequence.empty()) {
                fragments.push_back(make_fragment(
                    child,
                    read_name,
                    "softclip_right",
                    clip.ref_end,
                    InsertionFragmentSource::kClipRefRight,
                    start,
                    sequence));
            }
        }

        for (const auto& ins : find_long_insertions(read, min_insertion_len)) {
            if (!position_in_window(ins.ref_pos, window)) {
                continue;
            }
            const std::string sequence = read.decode_subsequence(ins.start, ins.len);
            if (sequence.empty()) {
                continue;
            }
            fragments.push_back(make_fragment(
                child,
                read_name,
                "ins",
                ins.ref_pos,
                InsertionFragmentSource::kCigarInsertion,
                ins.start,
                sequence));
        }

        if (!fragments.empty()) {
            const auto decisions = classify_fragments(child, bam_path, fragments);
            for (const auto& decision : decisions) {
                best = choose_better_decision(best, decision);
            }
        }

        if (best.match_class == EvidenceMatchClass::kNone) {
            const auto split_bp = nearest_breakpoint_in_window(read.pos(), clip.ref_end, window);
            if (((flag & BAM_FSUPPLEMENTARY) != 0 || has_sa_breakpoint_in_window(read, child.chrom, window)) &&
                split_bp.has_value()) {
                ReadEvidenceDecision ambiguous;
                ambiguous.match_class = EvidenceMatchClass::kAmbiguous;
                ambiguous.evidence.parent_bam_path = bam_path;
                ambiguous.evidence.read_name = read_name;
                ambiguous.evidence.evidence_type = "split_sa";
                ambiguous.evidence.breakpoint_pos = *split_bp;
                ambiguous.evidence.fragment_len = 0;
                ambiguous.evidence.veto_reason = "SPLIT_ALIGNMENT_NEAR_LOCUS";
                best = choose_better_decision(best, ambiguous);
            }
        }

        return best;
    }

    PipelineConfig pipeline_config;
    TEKmerQuickClassifierModule classifier;
    SplitSAFragmentModule split_sa_module;
    int32_t min_softclip_len = 50;
    int32_t min_insertion_len = 50;
    int32_t parent_mapq_min = 0;
};

ParentPoolScanner::ParentPoolScanner(DenovoConfig config) : config_(std::move(config)) {
    readers_.reserve(config_.parent_bam_paths.size());
    for (const auto& bam_path : config_.parent_bam_paths) {
        auto reader = std::make_unique<IndexedBamReader>(bam_path, config_.bam_threads);
        if (!reader->is_valid()) {
            throw std::runtime_error("Failed to open parent BAM: " + bam_path);
        }
        if (!reader->has_index()) {
            throw std::runtime_error("Missing BAM index for parent BAM: " + bam_path);
        }
        readers_.push_back(std::move(reader));
    }

    impl_ = std::make_unique<Impl>(config_);
    implementation_status_ =
        impl_->classifier.is_enabled()
        ? "TARGETED_PARENT_SCAN_WITH_SPLIT_FRAGMENT_CLASSIFICATION"
        : "TARGETED_PARENT_SCAN_WITH_SPLIT_NO_TEKMER_CLASSIFIER";
}

ParentPoolScanner::~ParentPoolScanner() = default;

ParentEvidenceSummary ParentPoolScanner::scan_candidate(const DenovoChildCandidate& candidate) const {
    ParentEvidenceSummary summary;
    summary.scanner_status = implementation_status_;
    if (readers_.empty()) {
        return summary;
    }

    const CandidateScanWindow window = build_scan_window(candidate, config_);
    std::unordered_map<std::string, ReadEvidenceDecision> best_by_read_key;

    for (size_t i = 0; i < readers_.size(); ++i) {
        const auto& reader = readers_[i];
        const std::string& bam_path = config_.parent_bam_paths[i];
        std::vector<BamRecordPtr> fetched_records;
        int32_t fetched_tid = -1;
        reader->fetch(candidate.chrom, window.fetch_start, window.fetch_end, [&](const bam1_t* record) {
            if (record != nullptr &&
                (record->core.flag & BAM_FUNMAP) == 0 &&
                (record->core.flag & BAM_FSECONDARY) == 0 &&
                record->core.qual >= impl_->parent_mapq_min) {
                BamRecordPtr copy(bam_dup1(record));
                if (copy) {
                    fetched_tid = record->core.tid;
                    fetched_records.push_back(std::move(copy));
                }
            }

            const ReadEvidenceDecision decision = impl_->classify_read(
                candidate,
                bam_path,
                record,
                window,
                config_.family_match_veto);
            if (decision.match_class == EvidenceMatchClass::kNone) {
                return true;
            }

            const std::string key = bam_path + "\n" + decision.evidence.read_name;
            const auto it = best_by_read_key.find(key);
            if (it == best_by_read_key.end()) {
                best_by_read_key.emplace(key, decision);
            } else {
                it->second = choose_better_decision(it->second, decision);
            }
            return true;
        });

        if (!fetched_records.empty()) {
            std::vector<const bam1_t*> raw_records;
            raw_records.reserve(fetched_records.size());
            for (const auto& record : fetched_records) {
                raw_records.push_back(record.get());
            }

            ComponentCall component;
            component.chrom = candidate.chrom;
            component.tid = fetched_tid;
            component.anchor_pos = candidate.pos;
            component.bin_start = window.fetch_start;
            component.bin_end = window.fetch_end;
            component.read_indices.resize(raw_records.size());
            std::iota(component.read_indices.begin(), component.read_indices.end(), 0);

            const auto split_fragments = impl_->split_sa_module.extract(component, raw_records);
            std::vector<InsertionFragment> filtered_split_fragments;
            filtered_split_fragments.reserve(split_fragments.size());
            for (const auto& fragment : split_fragments) {
                if (fragment.sequence.empty()) {
                    continue;
                }
                if (fragment.ref_junc_pos < 0 || !position_in_window(fragment.ref_junc_pos, window)) {
                    continue;
                }
                filtered_split_fragments.push_back(fragment);
            }

            const auto split_decisions = impl_->classify_fragments(
                candidate,
                bam_path,
                filtered_split_fragments);
            for (const auto& decision : split_decisions) {
                const std::string key = bam_path + "\n" + decision.evidence.read_name;
                const auto it = best_by_read_key.find(key);
                if (it == best_by_read_key.end()) {
                    best_by_read_key.emplace(key, decision);
                } else {
                    it->second = choose_better_decision(it->second, decision);
                }
            }
        }
    }

    std::unordered_set<std::string> support_bams;
    for (const auto& kv : best_by_read_key) {
        const auto& decision = kv.second;
        if (decision.match_class == EvidenceMatchClass::kExact) {
            summary.total_support_reads += 1;
            summary.exact_te_reads += 1;
            support_bams.insert(decision.evidence.parent_bam_path);
            summary.veto_reads.push_back(decision.evidence);
        } else if (decision.match_class == EvidenceMatchClass::kFamily) {
            summary.total_support_reads += 1;
            summary.family_te_reads += 1;
            support_bams.insert(decision.evidence.parent_bam_path);
            summary.veto_reads.push_back(decision.evidence);
        } else if (decision.match_class == EvidenceMatchClass::kAmbiguous) {
            summary.ambiguous_signal_reads += 1;
            summary.veto_reads.push_back(decision.evidence);
        }
    }

    summary.support_bams.assign(support_bams.begin(), support_bams.end());
    std::sort(summary.support_bams.begin(), summary.support_bams.end());
    return summary;
}

}  // namespace placer
