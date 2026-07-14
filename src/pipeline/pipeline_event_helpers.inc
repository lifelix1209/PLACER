int bool_as_int(bool value) {
    return value ? 1 : 0;
}

const char* final_hypothesis_kind_name(FinalHypothesisKind kind) {
    switch (kind) {
        case FinalHypothesisKind::kReference:
            return "REFERENCE";
        case FinalHypothesisKind::kInsertionNonTe:
            return "NON_TE_INSERTION";
        case FinalHypothesisKind::kTeUnknown:
            return "TE_UNKNOWN";
        case FinalHypothesisKind::kTeResolved:
            return "TE_RESOLVED";
    }
    return "UNKNOWN";
}

int32_t evidence_specificity_rank(EvidenceKind kind) {
    switch (kind) {
        case EvidenceKind::kIndel:
            return 0;
        case EvidenceKind::kSoftClip:
            return 1;
        case EvidenceKind::kSAHint:
            return 2;
    }
    return std::numeric_limits<int32_t>::max();
}

void accumulate_read_window_assignment_score(
    const EvidencePoint& point,
    double point_score,
    ReadWindowAssignmentScore& score) {
    score.valid = true;
    score.total_score += point_score;

    const int32_t specificity_rank = evidence_specificity_rank(point.kind);
    if (specificity_rank < score.specificity_rank) {
        score.specificity_rank = specificity_rank;
        score.specificity_score = point_score;
        return;
    }
    if (specificity_rank == score.specificity_rank &&
        point_score > score.specificity_score) {
        score.specificity_score = point_score;
    }
}

bool better_read_window_assignment_score(
    const ReadWindowAssignmentScore& lhs,
    const ReadWindowAssignmentScore& rhs) {
    if (!lhs.valid) {
        return false;
    }
    if (!rhs.valid) {
        return true;
    }
    if (lhs.specificity_rank != rhs.specificity_rank) {
        return lhs.specificity_rank < rhs.specificity_rank;
    }
    if (lhs.specificity_score != rhs.specificity_score) {
        return lhs.specificity_score > rhs.specificity_score;
    }
    return lhs.total_score > rhs.total_score;
}

bool better_event_flank_placement(
    const EventFlankPlacement& lhs,
    const EventFlankPlacement& rhs) {
    if (lhs.identity != rhs.identity) {
        return lhs.identity > rhs.identity;
    }
    if (lhs.breakpoint_delta != rhs.breakpoint_delta) {
        return lhs.breakpoint_delta < rhs.breakpoint_delta;
    }
    if (lhs.endpoint_offset != rhs.endpoint_offset) {
        return lhs.endpoint_offset < rhs.endpoint_offset;
    }
    if (lhs.align_len != rhs.align_len) {
        return lhs.align_len > rhs.align_len;
    }
    if (lhs.ref_start != rhs.ref_start) {
        return lhs.ref_start < rhs.ref_start;
    }
    return lhs.ref_end < rhs.ref_end;
}

std::string display_or_na(const std::string& value) {
    return value.empty() ? "NA" : value;
}

std::string event_segmentation_cache_key(
    const ComponentCall& component,
    const EventReadEvidence& event_evidence,
    const EventConsensus& event_consensus) {
    const int32_t bp_left = std::min(event_evidence.bp_left, event_evidence.bp_right);
    const int32_t bp_right = std::max(event_evidence.bp_left, event_evidence.bp_right);
    std::string key;
    key.reserve(
        component.chrom.size() +
        event_consensus.consensus_seq.size() +
        64);
    key.append(std::to_string(component.chrom.size()));
    key.push_back(':');
    key.append(component.chrom);
    key.push_back('|');
    key.append(std::to_string(bp_left));
    key.push_back('|');
    key.append(std::to_string(bp_right));
    key.push_back('|');
    key.append(std::to_string(event_consensus.consensus_seq.size()));
    key.push_back(':');
    key.append(event_consensus.consensus_seq);
    return key;
}

std::string summarize_seq_ends(const std::string& seq, size_t keep = 24) {
    if (seq.empty()) {
        return "NA";
    }
    if (seq.size() <= (2 * keep)) {
        return seq;
    }
    return seq.substr(0, keep) + "..." + seq.substr(seq.size() - keep);
}

bool is_primary_alignment(const bam1_t* record) {
    if (!record) {
        return false;
    }
    const uint16_t flag = record->core.flag;
    return (flag & BAM_FUNMAP) == 0 &&
           (flag & BAM_FSECONDARY) == 0 &&
           (flag & BAM_FSUPPLEMENTARY) == 0;
}

std::string record_qname(const bam1_t* record) {
    if (!record) {
        return {};
    }
    ReadView view(record);
    return std::string(view.qname());
}

struct LocalEventSignal {
    bool split = false;
    bool indel = false;
    bool left_clip = false;
    bool right_clip = false;
    int32_t split_left_pos = -1;
    int32_t split_right_pos = -1;
    int32_t indel_pos = -1;
    int32_t max_raw_cigar_insert_len = 0;
    int32_t left_clip_pos = -1;
    int32_t right_clip_pos = -1;

    bool any() const {
        return split || indel || left_clip || right_clip;
    }
};

LocalEventSignal classify_local_event_signal(
    const bam1_t* record,
    const std::string& chrom,
    int32_t window_start,
    int32_t window_end) {
    LocalEventSignal signal;
    if (!record) {
        return signal;
    }

    ReadView read(record);
    const std::vector<int32_t> split_positions =
        robust_local_split_insertion_positions(read, chrom);

    const uint32_t* cigar = read.cigar();
    const int32_t n_cigar = read.n_cigar();
    if (!cigar || n_cigar <= 0) {
        return signal;
    }

    int first = -1;
    int last = -1;
    for (int32_t i = 0; i < n_cigar; ++i) {
        if (bam_cigar_op(cigar[i]) != BAM_CHARD_CLIP) {
            first = i;
            break;
        }
    }
    for (int32_t i = n_cigar - 1; i >= 0; --i) {
        if (bam_cigar_op(cigar[i]) != BAM_CHARD_CLIP) {
            last = i;
            break;
        }
    }

    int32_t ref_pos = read.pos();
    const int32_t window_center = window_start + ((window_end - window_start) / 2);
    int32_t best_indel_dist = std::numeric_limits<int32_t>::max();
    for (int32_t i = 0; i < n_cigar; ++i) {
        const int op = bam_cigar_op(cigar[i]);
        const int32_t len = static_cast<int32_t>(bam_cigar_oplen(cigar[i]));

        if (i == first && op == BAM_CSOFT_CLIP && len >= kSoftClipSignalMin &&
            ref_pos >= window_start && ref_pos <= window_end) {
            signal.left_clip = true;
            signal.left_clip_pos = ref_pos;
        }
        if (i == last && op == BAM_CSOFT_CLIP && len >= kSoftClipSignalMin &&
            ref_pos >= window_start && ref_pos <= window_end) {
            signal.right_clip = true;
            signal.right_clip_pos = ref_pos;
        }
        if (op == BAM_CINS && ref_pos >= window_start && ref_pos <= window_end) {
            signal.max_raw_cigar_insert_len = std::max(
                signal.max_raw_cigar_insert_len,
                len);
            if (len < kLongInsertionSignalMin) {
                continue;
            }
            signal.indel = true;
            const int32_t dist = std::abs(ref_pos - window_center);
            if (dist < best_indel_dist) {
                best_indel_dist = dist;
                signal.indel_pos = ref_pos;
            }
        }

        if ((bam_cigar_type(op) & 2) != 0) {
            ref_pos += len;
        }
    }

    for (int32_t pos : split_positions) {
        if (pos < window_start || pos > window_end) {
            continue;
        }
        signal.split = true;
        if (signal.split_left_pos < 0 || pos < signal.split_left_pos) {
            signal.split_left_pos = pos;
        }
        if (signal.split_right_pos < 0 || pos > signal.split_right_pos) {
            signal.split_right_pos = pos;
        }
    }
    return signal;
}

bool read_has_local_event_signal(
    const bam1_t* record,
    const std::string& chrom,
    int32_t window_start,
    int32_t window_end) {
    return classify_local_event_signal(record, chrom, window_start, window_end).any();
}

std::vector<std::string> sorted_string_set(const std::unordered_set<std::string>& values) {
    std::vector<std::string> out;
    out.reserve(values.size());
    for (const auto& value : values) {
        out.push_back(value);
    }
    std::sort(out.begin(), out.end());
    return out;
}

std::string upper_acgt(const std::string& s);

std::string collect_aligned_query_bases_before(
    const ReadView& read,
    int32_t query_limit,
    int32_t max_bases) {
    if (query_limit <= 0 || max_bases <= 0) {
        return {};
    }
    const int32_t take_start = std::max(0, query_limit - max_bases);
    const int32_t take_len = query_limit - take_start;
    return upper_acgt(read.decode_subsequence(take_start, take_len));
}

std::string collect_aligned_query_bases_after(
    const ReadView& read,
    int32_t query_start,
    int32_t max_bases) {
    if (query_start < 0 || max_bases <= 0) {
        return {};
    }
    return upper_acgt(read.decode_subsequence(query_start, max_bases));
}

std::unordered_map<std::string, const bam1_t*> best_records_by_qname(
    const std::vector<const bam1_t*>& local_records) {
    std::unordered_map<std::string, const bam1_t*> records_by_qname;
    records_by_qname.reserve(local_records.size());
    for (const bam1_t* record : local_records) {
        if (!record) {
            continue;
        }
        const std::string qname = record_qname(record);
        if (qname.empty()) {
            continue;
        }

        auto it = records_by_qname.find(qname);
        if (it == records_by_qname.end()) {
            records_by_qname.emplace(qname, record);
            continue;
        }

        ReadView candidate(record);
        ReadView incumbent(it->second);
        const int32_t candidate_seq_len = candidate.seq_len();
        const int32_t incumbent_seq_len = incumbent.seq_len();
        if (candidate_seq_len > incumbent_seq_len ||
            (candidate_seq_len == incumbent_seq_len &&
             is_primary_alignment(record) &&
             !is_primary_alignment(it->second))) {
            it->second = record;
        }
    }
    return records_by_qname;
}

bool can_build_event_string_from_fragment(
    const bam1_t* record,
    const InsertionFragment& fragment) {
    if (!record || fragment.sequence.empty() || fragment.start < 0 || fragment.length <= 0) {
        return false;
    }

    ReadView read(record);
    const int32_t insert_start = fragment.start;
    const int32_t insert_end = fragment.start + fragment.length;
    if (insert_end > read.seq_len()) {
        return false;
    }

    const int32_t left_flank_len = insert_start - std::max(0, insert_start - kEventConsensusFlankBp);
    const int32_t right_flank_len = std::max(0, std::min(read.seq_len(), insert_end + kEventConsensusFlankBp) - insert_end);

    switch (fragment.source) {
        case InsertionFragmentSource::kClipRefLeft:
            return right_flank_len >= kEventConsensusMinAnchorBp;
        case InsertionFragmentSource::kClipRefRight:
            return left_flank_len >= kEventConsensusMinAnchorBp;
        case InsertionFragmentSource::kCigarInsertion:
        case InsertionFragmentSource::kSplitSa:
            return left_flank_len >= kEventConsensusMinAnchorBp &&
                right_flank_len >= kEventConsensusMinAnchorBp;
        case InsertionFragmentSource::kUnknown:
        default:
            return false;
    }
}

std::string build_event_string_from_fragment(
    const bam1_t* record,
    const InsertionFragment& fragment) {
    if (!can_build_event_string_from_fragment(record, fragment)) {
        return {};
    }

    ReadView read(record);
    const int32_t insert_start = fragment.start;
    const int32_t insert_end = fragment.start + fragment.length;
    if (insert_end > read.seq_len()) {
        return {};
    }

    const std::string insert_seq = upper_acgt(fragment.sequence);
    if (insert_seq.empty()) {
        return {};
    }

    const auto fetch_left_flank = [&](int32_t flank_len) {
        if (flank_len <= 0) {
            return std::string();
        }
        return collect_aligned_query_bases_before(read, insert_start, flank_len);
    };
    const auto fetch_right_flank = [&](int32_t flank_len) {
        if (flank_len <= 0) {
            return std::string();
        }
        return collect_aligned_query_bases_after(read, insert_end, flank_len);
    };

    switch (fragment.source) {
        case InsertionFragmentSource::kClipRefLeft: {
            const std::string right_flank = fetch_right_flank(kEventConsensusFlankBp);
            if (static_cast<int32_t>(right_flank.size()) < kEventConsensusMinAnchorBp) {
                return {};
            }
            return insert_seq + right_flank;
        }
        case InsertionFragmentSource::kClipRefRight: {
            const std::string left_flank = fetch_left_flank(kEventConsensusFlankBp);
            if (static_cast<int32_t>(left_flank.size()) < kEventConsensusMinAnchorBp) {
                return {};
            }
            return left_flank + insert_seq;
        }
        case InsertionFragmentSource::kCigarInsertion:
        case InsertionFragmentSource::kSplitSa: {
            const std::string left_flank = fetch_left_flank(kEventConsensusFlankBp);
            const std::string right_flank = fetch_right_flank(kEventConsensusFlankBp);
            if (static_cast<int32_t>(left_flank.size()) < kEventConsensusMinAnchorBp ||
                static_cast<int32_t>(right_flank.size()) < kEventConsensusMinAnchorBp) {
                return {};
            }
            return left_flank + insert_seq + right_flank;
        }
        case InsertionFragmentSource::kUnknown:
        default:
            return {};
    }
}

bool fragment_supports_event_consensus(
    const InsertionFragment& fragment) {
    switch (fragment.source) {
        case InsertionFragmentSource::kClipRefLeft:
        case InsertionFragmentSource::kClipRefRight:
        case InsertionFragmentSource::kCigarInsertion:
        case InsertionFragmentSource::kSplitSa:
            return true;
        case InsertionFragmentSource::kUnknown:
        default:
            return false;
    }
}

bool fragment_has_full_event_context(
    const InsertionFragment& fragment) {
    switch (fragment.source) {
        case InsertionFragmentSource::kCigarInsertion:
        case InsertionFragmentSource::kSplitSa:
            return true;
        case InsertionFragmentSource::kClipRefLeft:
        case InsertionFragmentSource::kClipRefRight:
        case InsertionFragmentSource::kUnknown:
        default:
            return false;
    }
}

bool fragment_is_local_to_event(
    const InsertionFragment& fragment,
    int32_t bp_left,
    int32_t bp_right) {
    if (fragment.ref_junc_pos < 0 || bp_left < 0 || bp_right < 0) {
        return false;
    }
    constexpr int32_t kEventFragmentJunctionSlackBp = 25;
    const int32_t span_start = std::max(0, std::min(bp_left, bp_right) - kEventFragmentJunctionSlackBp);
    const int32_t span_end = std::max(span_start, std::max(bp_left, bp_right) + kEventFragmentJunctionSlackBp);
    return fragment.ref_junc_pos >= span_start && fragment.ref_junc_pos <= span_end;
}

void emit_pipeline_log_line(const std::string& line) {
    static std::mutex log_mutex;
    const std::lock_guard<std::mutex> lock(log_mutex);
    std::cerr << line << '\n';
}

using SharedRecordBatch = std::shared_ptr<const std::vector<BufferedRecord>>;

int32_t exact_bin_index_for_record(const bam1_t* record, int32_t bin_size) {
    if (!record || bin_size <= 0) {
        return 0;
    }
    ReadView view(record);
    if (view.pos() < 0) {
        return 0;
    }
    return view.pos() / bin_size;
}

std::vector<const bam1_t*> collect_task_record_ptrs(const ExactBinTask& task) {
    std::vector<const bam1_t*> bin_records;
    bin_records.reserve(task.records.size());
    for (const auto& record_ref : task.records) {
        if (record_ref.record) {
            bin_records.push_back(record_ref.record);
        }
    }
    return bin_records;
}

void append_pipeline_performance(
    PipelinePerformanceSummary& dst,
    const PipelinePerformanceSummary& src) {
    dst.pipeline_wall_seconds += src.pipeline_wall_seconds;
    dst.pipeline_construction_seconds += src.pipeline_construction_seconds;
    dst.bam_stream_wall_seconds += src.bam_stream_wall_seconds;
    dst.bin_processing_wall_seconds += src.bin_processing_wall_seconds;
    dst.finalization_seconds += src.finalization_seconds;
    dst.gate1_preliminary_seconds += src.gate1_preliminary_seconds;
    dst.exact_bin_scan_seconds += src.exact_bin_scan_seconds;
    dst.component_build_seconds += src.component_build_seconds;
    dst.local_fetch_seconds += src.local_fetch_seconds;
    dst.projection_seconds += src.projection_seconds;
    dst.local_component_refresh_seconds += src.local_component_refresh_seconds;
    dst.fragment_extract_seconds += src.fragment_extract_seconds;
    dst.signal_cache_build_seconds += src.signal_cache_build_seconds;
    dst.hypothesis_summary_seconds += src.hypothesis_summary_seconds;
    dst.validator_input_count_seconds += src.validator_input_count_seconds;
    dst.event_consensus_seconds += src.event_consensus_seconds;
    dst.event_segmentation_seconds += src.event_segmentation_seconds;
    dst.te_alignment_seconds += src.te_alignment_seconds;
    dst.te_blast_batches += src.te_blast_batches;
    dst.te_blast_queries += src.te_blast_queries;
    dst.te_blast_cache_hits += src.te_blast_cache_hits;
    dst.te_blast_cache_misses += src.te_blast_cache_misses;
    dst.te_blast_deduplicated_queries += src.te_blast_deduplicated_queries;
    dst.segmentation_edit_distance_calls += src.segmentation_edit_distance_calls;
    dst.segmentation_edit_distance_cache_hits +=
        src.segmentation_edit_distance_cache_hits;
    dst.segmentation_edit_distance_cache_misses +=
        src.segmentation_edit_distance_cache_misses;
    dst.segmentation_seed_bins_total += src.segmentation_seed_bins_total;
    dst.segmentation_paired_searches += src.segmentation_paired_searches;
    dst.segmentation_one_sided_searches += src.segmentation_one_sided_searches;
    dst.segmentation_endpoint_slack_searches +=
        src.segmentation_endpoint_slack_searches;
    dst.segmentation_reverse_complement_retries +=
        src.segmentation_reverse_complement_retries;
    dst.segmentation_cache_hits += src.segmentation_cache_hits;
    dst.segmentation_cache_misses += src.segmentation_cache_misses;
}

void append_partial_pipeline_result(PipelineResult& dst, PipelineResult&& src) {
    dst.processed_bins += src.processed_bins;
    dst.built_components += src.built_components;
    dst.event_consensus_calls += src.event_consensus_calls;
    dst.genotype_calls += src.genotype_calls;
    append_pipeline_performance(dst.performance, src.performance);
    if (!src.final_calls.empty()) {
        dst.final_calls.insert(
            dst.final_calls.end(),
            std::make_move_iterator(src.final_calls.begin()),
            std::make_move_iterator(src.final_calls.end()));
    }
    if (!src.evidence_ledger.empty()) {
        dst.evidence_ledger.insert(
            dst.evidence_ledger.end(),
            std::make_move_iterator(src.evidence_ledger.begin()),
            std::make_move_iterator(src.evidence_ledger.end()));
    }
}

bool final_call_sort_less(const FinalCall& a, const FinalCall& b) {
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
}

bool same_call_locus(const FinalCall& a, const FinalCall& b) {
    if (a.tid != b.tid) {
        return false;
    }
    if (!a.chrom.empty() && !b.chrom.empty() && a.chrom != b.chrom) {
        return false;
    }
    const auto call_interval = [](const FinalCall& call) {
        int32_t left = -1;
        int32_t right = -1;
        if (call.bp_left >= 0 && call.bp_right >= 0) {
            left = std::min(call.bp_left, call.bp_right);
            right = std::max(call.bp_left, call.bp_right);
        } else if (call.pos >= 0) {
            left = call.pos;
            right = call.pos;
        }
        return std::pair<int32_t, int32_t>{left, right};
    };
    const auto [a_left, a_right] = call_interval(a);
    const auto [b_left, b_right] = call_interval(b);
    if (a_left >= 0 && a_right >= 0 && b_left >= 0 && b_right >= 0) {
        const bool a_has_interval = a_right > a_left;
        const bool b_has_interval = b_right > b_left;
        if (a_has_interval && b_has_interval) {
            return a_right >= b_left && b_right >= a_left;
        }
        if (a_right < b_left) {
            return b_left - a_right <= kFinalCallDedupDistanceBp;
        }
        if (b_right < a_left) {
            return a_left - b_right <= kFinalCallDedupDistanceBp;
        }
        return true;
    }
    return a.pos >= 0 && b.pos >= 0 &&
           std::abs(a.pos - b.pos) <= kFinalCallDedupDistanceBp;
}

bool final_call_has_robust_te_lfdr(const FinalCall& call) {
    return call.robust_mechanistic_qc == "PASS_TE_LFDR" &&
           call.robust_mechanistic_worst_case_lfdr <= 0.10;
}

int32_t final_call_precision_rank(const FinalCall& call) {
    if (call.final_qc.rfind("PASS_TE_CLOSED", 0) == 0) {
        return 4;
    }
    if (call.final_qc.rfind("PASS_STRUCTURAL_INSERTION", 0) == 0 &&
        call.final_qc.find("EVENT_CLUSTER_PROMOTED") == std::string::npos) {
        return 3;
    }
    if (call.final_qc.rfind("PASS_TE_IMPRECISE", 0) == 0 &&
        call.final_qc.find("EVENT_CLUSTER_PROMOTED") == std::string::npos) {
        return 2;
    }
    if (call.final_qc.find("EVENT_CLUSTER_PROMOTED") != std::string::npos) {
        return 1;
    }
    return 0;
}

bool prefer_new_call(const FinalCall& cur, const FinalCall& incumbent) {
    const bool cur_stable_event_community =
        cur.final_qc.find("EVENT_COMMUNITY_STABLE") != std::string::npos;
    const bool incumbent_stable_event_community =
        incumbent.final_qc.find("EVENT_COMMUNITY_STABLE") != std::string::npos;
    if (cur_stable_event_community != incumbent_stable_event_community) {
        return cur_stable_event_community;
    }
    const bool cur_robust_te = final_call_has_robust_te_lfdr(cur);
    const bool incumbent_robust_te = final_call_has_robust_te_lfdr(incumbent);
    if (cur_robust_te != incumbent_robust_te) {
        return cur_robust_te;
    }
    if (cur_robust_te && incumbent_robust_te &&
        cur.robust_mechanistic_worst_case_lfdr !=
            incumbent.robust_mechanistic_worst_case_lfdr) {
        return cur.robust_mechanistic_worst_case_lfdr <
               incumbent.robust_mechanistic_worst_case_lfdr;
    }
    const int32_t cur_precision = final_call_precision_rank(cur);
    const int32_t incumbent_precision = final_call_precision_rank(incumbent);
    if (cur_precision != incumbent_precision) {
        return cur_precision > incumbent_precision;
    }
    if (cur.support_reads != incumbent.support_reads) {
        return cur.support_reads > incumbent.support_reads;
    }
    if (cur.gq != incumbent.gq) {
        return cur.gq > incumbent.gq;
    }
    if (cur.cross_family_margin != incumbent.cross_family_margin) {
        return cur.cross_family_margin > incumbent.cross_family_margin;
    }
    if (cur.best_te_identity != incumbent.best_te_identity) {
        return cur.best_te_identity > incumbent.best_te_identity;
    }
    if (cur.te_name != incumbent.te_name) {
        if (incumbent.te_name.empty() || incumbent.te_name == "UNK") {
            return true;
        }
        if (cur.te_name.empty() || cur.te_name == "UNK") {
            return false;
        }
    }
    if (cur.event_consensus_len != incumbent.event_consensus_len) {
        return cur.event_consensus_len > incumbent.event_consensus_len;
    }
    if (cur.best_te_query_coverage != incumbent.best_te_query_coverage) {
        return cur.best_te_query_coverage > incumbent.best_te_query_coverage;
    }
    return cur.pos < incumbent.pos;
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
