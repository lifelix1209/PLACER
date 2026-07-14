EventConsensus Pipeline::build_event_consensus(
    const ComponentCall& component,
    const std::vector<const bam1_t*>& local_records,
    const std::vector<InsertionFragment>& fragments,
    const EventReadEvidence& event_evidence,
    ConsensusContextMode mode) const {
    EventConsensus consensus;
    (void)component;

    const ConsensusInputSummary inputs = collect_event_consensus_inputs(
        local_records,
        fragments,
        event_evidence);
    bool use_full_context = false;
    switch (mode) {
        case ConsensusContextMode::kPartialOnly:
            use_full_context = false;
            break;
        case ConsensusContextMode::kPreferFull:
        default:
            use_full_context = !inputs.full_event_by_qname.empty();
            break;
    }
    consensus.full_context_input_reads = inputs.full_context_input_reads;
    consensus.partial_context_input_reads = inputs.partial_context_input_reads;
    consensus.left_anchor_input_reads = inputs.left_anchor_input_reads;
    consensus.right_anchor_input_reads = inputs.right_anchor_input_reads;
    consensus.used_full_context = use_full_context;
    const auto& event_by_qname = use_full_context
        ? inputs.full_event_by_qname
        : inputs.partial_event_by_qname;
    std::vector<std::string> event_strings;
    event_strings.reserve(event_by_qname.size());
    for (const auto& row : event_by_qname) {
        event_strings.push_back(row.second);
    }

    consensus.input_event_reads = static_cast<int32_t>(event_strings.size());
    if (consensus.input_event_reads <= 0) {
        consensus.qc_reason = "NO_EVENT_STRING_READS";
        return consensus;
    }
    const int32_t min_event_reads = use_full_context
        ? 1
        : std::max(1, config_.event_consensus_poa_min_reads);
    if (consensus.input_event_reads < min_event_reads) {
        consensus.qc_reason = "INSUFFICIENT_EVENT_READS";
        return consensus;
    }

    std::sort(event_strings.begin(), event_strings.end(), [](const std::string& a, const std::string& b) {
        if (a.size() != b.size()) {
            return a.size() > b.size();
        }
        return a < b;
    });
    if (static_cast<int32_t>(event_strings.size()) > config_.event_consensus_poa_max_reads) {
        event_strings.resize(static_cast<size_t>(config_.event_consensus_poa_max_reads));
    }

    consensus.consensus_seq = upper_acgt(abpoa_consensus(event_strings));
    consensus.consensus_len = static_cast<int32_t>(consensus.consensus_seq.size());
    if (consensus.consensus_seq.empty()) {
        consensus.qc_reason = "EMPTY_EVENT_CONSENSUS";
        return consensus;
    }

    consensus.qc_pass = true;
    consensus.qc_reason = mode == ConsensusContextMode::kPartialOnly
        ? "PASS_EVENT_CONSENSUS_PARTIAL_RETRY"
        : "PASS_EVENT_CONSENSUS";
    return consensus;
}

Pipeline::ConsensusInputSummary Pipeline::collect_event_consensus_inputs(
    const std::vector<const bam1_t*>& local_records,
    const std::vector<InsertionFragment>& fragments,
    const EventReadEvidence& event_evidence) const {
    ConsensusInputSummary summary;

    std::unordered_set<std::string> support_qnames(
        event_evidence.support_qnames.begin(),
        event_evidence.support_qnames.end());
    const std::unordered_map<std::string, const bam1_t*> records_by_qname =
        best_records_by_qname(local_records);

    std::unordered_set<std::string> full_context_qnames;
    std::unordered_set<std::string> partial_context_qnames;
    std::unordered_set<std::string> left_anchor_qnames;
    std::unordered_set<std::string> right_anchor_qnames;
    summary.full_event_by_qname.reserve(fragments.size());
    summary.partial_event_by_qname.reserve(fragments.size());
    full_context_qnames.reserve(fragments.size());
    partial_context_qnames.reserve(fragments.size());
    left_anchor_qnames.reserve(fragments.size());
    right_anchor_qnames.reserve(fragments.size());

    for (const auto& fragment : fragments) {
        if (fragment.read_id.empty() ||
            support_qnames.find(fragment.read_id) == support_qnames.end() ||
            !fragment_supports_event_consensus(fragment) ||
            !fragment_is_local_to_event(fragment, event_evidence.bp_left, event_evidence.bp_right)) {
            continue;
        }
        const auto record_it = records_by_qname.find(fragment.read_id);
        if (record_it == records_by_qname.end()) {
            continue;
        }
        if (!can_build_event_string_from_fragment(record_it->second, fragment)) {
            continue;
        }
        const std::string event_string = build_event_string_from_fragment(
            record_it->second,
            fragment);
        if (event_string.empty()) {
            continue;
        }

        switch (fragment.source) {
            case InsertionFragmentSource::kCigarInsertion:
            case InsertionFragmentSource::kSplitSa:
                full_context_qnames.insert(fragment.read_id);
                left_anchor_qnames.insert(fragment.read_id);
                right_anchor_qnames.insert(fragment.read_id);
                break;
            case InsertionFragmentSource::kClipRefLeft:
                partial_context_qnames.insert(fragment.read_id);
                right_anchor_qnames.insert(fragment.read_id);
                break;
            case InsertionFragmentSource::kClipRefRight:
                partial_context_qnames.insert(fragment.read_id);
                left_anchor_qnames.insert(fragment.read_id);
                break;
            case InsertionFragmentSource::kUnknown:
            default:
                break;
        }

        auto& event_by_qname = fragment_has_full_event_context(fragment)
            ? summary.full_event_by_qname
            : summary.partial_event_by_qname;
        auto it = event_by_qname.find(fragment.read_id);
        if (it == event_by_qname.end() || event_string.size() > it->second.size()) {
            event_by_qname[fragment.read_id] = event_string;
        }
    }

    summary.full_context_input_reads = static_cast<int32_t>(full_context_qnames.size());
    summary.partial_context_input_reads = static_cast<int32_t>(partial_context_qnames.size());
    summary.left_anchor_input_reads = static_cast<int32_t>(left_anchor_qnames.size());
    summary.right_anchor_input_reads = static_cast<int32_t>(right_anchor_qnames.size());
    return summary;
}

Pipeline::ConsensusInputCounts Pipeline::collect_event_consensus_input_counts(
    const std::vector<const bam1_t*>& local_records,
    const std::vector<InsertionFragment>& fragments,
    const EventReadEvidence& event_evidence) const {
    ConsensusInputCounts counts;

    std::unordered_set<std::string> support_qnames(
        event_evidence.support_qnames.begin(),
        event_evidence.support_qnames.end());
    const std::unordered_map<std::string, const bam1_t*> records_by_qname =
        best_records_by_qname(local_records);

    std::unordered_set<std::string> full_context_qnames;
    std::unordered_set<std::string> partial_context_qnames;
    std::unordered_set<std::string> left_anchor_qnames;
    std::unordered_set<std::string> right_anchor_qnames;
    full_context_qnames.reserve(fragments.size());
    partial_context_qnames.reserve(fragments.size());
    left_anchor_qnames.reserve(fragments.size());
    right_anchor_qnames.reserve(fragments.size());

    for (const auto& fragment : fragments) {
        if (fragment.read_id.empty() ||
            support_qnames.find(fragment.read_id) == support_qnames.end() ||
            !fragment_supports_event_consensus(fragment) ||
            !fragment_is_local_to_event(fragment, event_evidence.bp_left, event_evidence.bp_right)) {
            continue;
        }
        const auto record_it = records_by_qname.find(fragment.read_id);
        if (record_it == records_by_qname.end() ||
            !can_build_event_string_from_fragment(record_it->second, fragment)) {
            continue;
        }

        switch (fragment.source) {
            case InsertionFragmentSource::kCigarInsertion:
            case InsertionFragmentSource::kSplitSa:
                full_context_qnames.insert(fragment.read_id);
                left_anchor_qnames.insert(fragment.read_id);
                right_anchor_qnames.insert(fragment.read_id);
                break;
            case InsertionFragmentSource::kClipRefLeft:
                partial_context_qnames.insert(fragment.read_id);
                right_anchor_qnames.insert(fragment.read_id);
                break;
            case InsertionFragmentSource::kClipRefRight:
                partial_context_qnames.insert(fragment.read_id);
                left_anchor_qnames.insert(fragment.read_id);
                break;
            case InsertionFragmentSource::kUnknown:
            default:
                break;
        }
    }

    counts.full_context_input_reads = static_cast<int32_t>(full_context_qnames.size());
    counts.partial_context_input_reads = static_cast<int32_t>(partial_context_qnames.size());
    counts.left_anchor_input_reads = static_cast<int32_t>(left_anchor_qnames.size());
    counts.right_anchor_input_reads = static_cast<int32_t>(right_anchor_qnames.size());
    counts.input_event_reads =
        counts.full_context_input_reads > 0
            ? counts.full_context_input_reads
            : counts.partial_context_input_reads;
    return counts;
}

Pipeline::HypothesisValidatorEvidence Pipeline::collect_hypothesis_validator_evidence(
    const HypothesisSummary& summary,
    const ConsensusInputCounts& inputs,
    int32_t anchor_pos) const {
    HypothesisValidatorEvidence out;
    out.summary = summary;
    out.event_evidence = summary.event_evidence;
    out.precise_support = summary.alt_split_reads + summary.alt_indel_reads;
    out.breakpoint_width = std::max(0, summary.bp_right - summary.bp_left);
    out.anchor_distance = std::min(
        std::abs(summary.bp_left - anchor_pos),
        std::abs(summary.bp_right - anchor_pos));
    out.full_context_input_reads = inputs.full_context_input_reads;
    out.partial_context_input_reads = inputs.partial_context_input_reads;
    out.left_anchor_input_reads = inputs.left_anchor_input_reads;
    out.right_anchor_input_reads = inputs.right_anchor_input_reads;
    out.input_event_reads = inputs.input_event_reads;

    const bool has_bilateral_anchor =
        out.left_anchor_input_reads > 0 && out.right_anchor_input_reads > 0;
    const bool has_precise_or_full_context =
        out.full_context_input_reads > 0 || out.precise_support > 0;
    const bool has_identifiable_bilateral_partial_context =
        has_identifiable_bilateral_partial_context_support(
            out.left_anchor_input_reads,
            out.right_anchor_input_reads,
            out.partial_context_input_reads,
            out.input_event_reads);
    const bool has_validator_support =
        has_precise_or_full_context || has_identifiable_bilateral_partial_context;

    if (!has_bilateral_anchor && !has_precise_or_full_context) {
        out.qc_reason = "VALIDATOR_NO_BILATERAL_ANCHOR";
        return out;
    }
    if (!has_validator_support) {
        out.qc_reason = "VALIDATOR_NO_PRECISE_OR_FULL_CONTEXT";
        return out;
    }

    out.feasible_for_expensive_stage = true;
    out.qc_reason = "VALIDATOR_PASS";
    return out;
}

bool Pipeline::has_identifiable_bilateral_partial_context_support(
    int32_t left_anchor_input_reads,
    int32_t right_anchor_input_reads,
    int32_t partial_context_input_reads,
    int32_t input_event_reads) const {
    const bool has_bilateral_anchor =
        left_anchor_input_reads > 0 && right_anchor_input_reads > 0;
    return has_bilateral_anchor &&
        partial_context_input_reads >= 2 &&
        input_event_reads >= 2;
}

bool Pipeline::compare_hypothesis_validator_priority(
    const HypothesisValidatorEvidence& lhs,
    const HypothesisValidatorEvidence& rhs) const {
    if (lhs.precise_support != rhs.precise_support) {
        return lhs.precise_support > rhs.precise_support;
    }
    if (lhs.full_context_input_reads != rhs.full_context_input_reads) {
        return lhs.full_context_input_reads > rhs.full_context_input_reads;
    }
    const int32_t lhs_anchor_balance =
        std::min(lhs.left_anchor_input_reads, lhs.right_anchor_input_reads);
    const int32_t rhs_anchor_balance =
        std::min(rhs.left_anchor_input_reads, rhs.right_anchor_input_reads);
    if (lhs_anchor_balance != rhs_anchor_balance) {
        return lhs_anchor_balance > rhs_anchor_balance;
    }
    if (lhs.breakpoint_width != rhs.breakpoint_width) {
        return lhs.breakpoint_width < rhs.breakpoint_width;
    }
    if (lhs.anchor_distance != rhs.anchor_distance) {
        return lhs.anchor_distance < rhs.anchor_distance;
    }
    if (lhs.summary.alt_struct_reads != rhs.summary.alt_struct_reads) {
        return lhs.summary.alt_struct_reads > rhs.summary.alt_struct_reads;
    }
    return lhs.summary.original_index < rhs.summary.original_index;
}

std::vector<Pipeline::ShortlistedHypothesis> Pipeline::build_expensive_stage_shortlist(
    const std::vector<HypothesisValidatorEvidence>& candidates) const {
    constexpr int32_t kRescuePreciseAnchorMinDistanceBp = 500;
    std::vector<HypothesisValidatorEvidence> feasible;
    feasible.reserve(candidates.size());
    for (const auto& candidate : candidates) {
        if (candidate.feasible_for_expensive_stage) {
            feasible.push_back(candidate);
        }
    }
    if (feasible.empty()) {
        return {};
    }

    std::sort(feasible.begin(), feasible.end(), [&](const auto& lhs, const auto& rhs) {
        return compare_hypothesis_validator_priority(lhs, rhs);
    });

    std::vector<ShortlistedHypothesis> out;
    out.push_back(ShortlistedHypothesis{feasible.front(), true});

    const auto support_jaccard = [](
                                     const std::vector<std::string>& lhs,
                                     const std::vector<std::string>& rhs) {
        size_t i = 0;
        size_t j = 0;
        size_t intersect = 0;
        while (i < lhs.size() && j < rhs.size()) {
            if (lhs[i] == rhs[j]) {
                ++intersect;
                ++i;
                ++j;
            } else if (lhs[i] < rhs[j]) {
                ++i;
            } else {
                ++j;
            }
        }
        const size_t union_size = lhs.size() + rhs.size() - intersect;
        return union_size == 0
            ? 0.0
            : static_cast<double>(intersect) / static_cast<double>(union_size);
    };

    for (size_t i = 1; i < feasible.size(); ++i) {
        const auto& candidate = feasible[i];
        bool redundant = false;
        for (const auto& selected : out) {
            const auto& incumbent = selected.validator;
            const bool spatially_distinct =
                std::abs(candidate.summary.bp_left - incumbent.summary.bp_left) > 60 ||
                std::abs(candidate.summary.bp_right - incumbent.summary.bp_right) > 60;
            const bool support_distinct =
                support_jaccard(
                    candidate.summary.support_qnames,
                    incumbent.summary.support_qnames) < 0.5;
            if (!spatially_distinct && !support_distinct) {
                redundant = true;
                break;
            }
        }
        if (redundant) {
            continue;
        }
        out.push_back(ShortlistedHypothesis{candidate, false});
    }

    const HypothesisValidatorEvidence* rescue_precise_anchor = nullptr;
    for (const auto& candidate : feasible) {
        const auto& summary = candidate.summary;
        if (summary.hypothesis_support != 1 ||
            (summary.hypothesis_priority != 0 &&
             summary.hypothesis_priority != 1) ||
            summary.bp_left != summary.bp_right) {
            continue;
        }
        bool already_selected = false;
        bool spatially_distinct = true;
        for (const auto& selected : out) {
            const auto& incumbent = selected.validator.summary;
            if (incumbent.original_index == summary.original_index ||
                (incumbent.bp_left == summary.bp_left &&
                 incumbent.bp_right == summary.bp_right)) {
                already_selected = true;
                break;
            }
            const int32_t incumbent_center =
                incumbent.bp_left +
                ((incumbent.bp_right - incumbent.bp_left) / 2);
            if (std::abs(summary.bp_left - incumbent_center) <
                kRescuePreciseAnchorMinDistanceBp) {
                spatially_distinct = false;
                break;
            }
        }
        if (already_selected || !spatially_distinct) {
            continue;
        }
        if (rescue_precise_anchor == nullptr ||
            summary.hypothesis_priority <
                rescue_precise_anchor->summary.hypothesis_priority ||
            (summary.hypothesis_priority ==
                 rescue_precise_anchor->summary.hypothesis_priority &&
             summary.hypothesis_score >
                 rescue_precise_anchor->summary.hypothesis_score) ||
            (summary.hypothesis_priority ==
                 rescue_precise_anchor->summary.hypothesis_priority &&
             summary.hypothesis_score ==
                 rescue_precise_anchor->summary.hypothesis_score &&
             summary.bp_left < rescue_precise_anchor->summary.bp_left)) {
            rescue_precise_anchor = &candidate;
        }
    }
    if (rescue_precise_anchor != nullptr) {
        out.push_back(ShortlistedHypothesis{*rescue_precise_anchor, false});
    }
    return out;
}
