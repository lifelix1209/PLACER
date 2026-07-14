std::vector<int32_t> collect_alt_observed_lengths(
    const ComponentCall& component,
    const EventReadEvidence& event_evidence) {
    std::unordered_set<std::string> support_qnames(
        event_evidence.support_qnames.begin(),
        event_evidence.support_qnames.end());
    std::vector<int32_t> lengths;
    lengths.reserve(component.breakpoint_candidates.size());
    for (const auto& bp : component.breakpoint_candidates) {
        if (bp.read_id.empty() || support_qnames.find(bp.read_id) == support_qnames.end()) {
            continue;
        }
        if (bp.ins_len > 0) {
            lengths.push_back(bp.ins_len);
        } else if (bp.clip_len > 0) {
            lengths.push_back(bp.clip_len);
        }
    }
    return lengths;
}

int32_t infer_event_length_from_alt_support(const std::vector<int32_t>& lengths) {
    return median_i32(lengths);
}

double support_jaccard(
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
}

Pipeline::HypothesisSummary Pipeline::build_hypothesis_summary(
    const ComponentCall& component,
    const ComponentSignalCache& signal_cache,
    const std::vector<InsertionFragment>& fragments,
    size_t original_index,
    int32_t bp_left,
    int32_t bp_right,
    int32_t hypothesis_support,
    int32_t hypothesis_priority) const {
    const EventReadEvidence event_evidence = collect_event_read_evidence_for_bounds_cached(
        component,
        signal_cache,
        fragments,
        bp_left,
        bp_right,
        bp_left,
        bp_right);

    HypothesisSummary summary;
    summary.original_index = original_index;
    summary.event_evidence = event_evidence;
    summary.bp_left = event_evidence.bp_left;
    summary.bp_right = event_evidence.bp_right;
    summary.alt_split_reads = event_evidence.alt_split_reads;
    summary.alt_indel_reads = event_evidence.alt_indel_reads;
    summary.alt_left_clip_reads = event_evidence.alt_left_clip_reads;
    summary.alt_right_clip_reads = event_evidence.alt_right_clip_reads;
    summary.alt_struct_reads = event_evidence.alt_struct_reads;
    summary.raw_cigar_insert_reads = event_evidence.raw_cigar_insert_reads;
    summary.max_raw_cigar_insert_len = event_evidence.max_raw_cigar_insert_len;
    summary.ref_span_reads = event_evidence.ref_span_reads;
    summary.inferred_event_length = infer_event_length_from_alt_support(
        collect_alt_observed_lengths(component, event_evidence));
    summary.hypothesis_support = hypothesis_support;
    summary.hypothesis_priority = hypothesis_priority;
    summary.hypothesis_score =
        static_cast<double>(hypothesis_support) *
        static_cast<double>(breakpoint_hypothesis_support_weight(hypothesis_priority));
    summary.support_qnames = event_evidence.support_qnames;
    return summary;
}

std::vector<Pipeline::HypothesisSummary> Pipeline::select_hypothesis_summaries_for_expensive_stage(
    const std::vector<HypothesisSummary>& summaries) const {
    const auto collapsed = collapse_hypothesis_summaries(summaries);
    std::vector<HypothesisSummary> selected;
    selected.reserve(collapsed.size());
    for (size_t i = 0; i < collapsed.size(); ++i) {
        if (should_keep_hypothesis_for_expensive_stage(collapsed[i], i == 0)) {
            selected.push_back(collapsed[i]);
        }
    }
    return selected;
}

std::vector<Pipeline::HypothesisSummary> Pipeline::collapse_hypothesis_summaries(
    const std::vector<HypothesisSummary>& summaries) const {
    constexpr int32_t kHypothesisCollapseSlackBp = 30;
    constexpr double kHypothesisCollapseMinJaccard = 0.8;

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

    const auto better_summary_representative = [](
                                                   const HypothesisSummary& lhs,
                                                   const HypothesisSummary& rhs) {
        if (lhs.alt_struct_reads != rhs.alt_struct_reads) {
            return lhs.alt_struct_reads > rhs.alt_struct_reads;
        }
        const int32_t lhs_precise = lhs.alt_split_reads + lhs.alt_indel_reads;
        const int32_t rhs_precise = rhs.alt_split_reads + rhs.alt_indel_reads;
        if (lhs_precise != rhs_precise) {
            return lhs_precise > rhs_precise;
        }
        const int32_t lhs_width = std::abs(lhs.bp_right - lhs.bp_left);
        const int32_t rhs_width = std::abs(rhs.bp_right - rhs.bp_left);
        if (lhs_width != rhs_width) {
            return lhs_width < rhs_width;
        }
        return lhs.original_index < rhs.original_index;
    };

    std::vector<HypothesisSummary> kept;
    kept.reserve(summaries.size());
    for (const auto& summary : summaries) {
        bool merged = false;
        for (auto& incumbent : kept) {
            const bool close_left =
                std::abs(summary.bp_left - incumbent.bp_left) <= kHypothesisCollapseSlackBp;
            const bool close_right =
                std::abs(summary.bp_right - incumbent.bp_right) <= kHypothesisCollapseSlackBp;
            if (!close_left || !close_right) {
                continue;
            }
            if (support_jaccard(summary.support_qnames, incumbent.support_qnames) <
                kHypothesisCollapseMinJaccard) {
                continue;
            }
            if (better_summary_representative(summary, incumbent)) {
                incumbent = summary;
            }
            merged = true;
            break;
        }
        if (!merged) {
            kept.push_back(summary);
        }
    }
    return kept;
}

bool Pipeline::should_keep_hypothesis_for_expensive_stage(
    const HypothesisSummary& summary,
    bool is_top_ranked_survivor) const {
    if (is_top_ranked_survivor) {
        return true;
    }
    if (summary.alt_split_reads > 0 || summary.alt_indel_reads > 0) {
        return true;
    }
    const bool bilateral_clip =
        summary.alt_left_clip_reads > 0 && summary.alt_right_clip_reads > 0;
    if (bilateral_clip && summary.alt_struct_reads >= 2) {
        return true;
    }
    const bool unilateral_clip =
        (summary.alt_left_clip_reads > 0) != (summary.alt_right_clip_reads > 0);
    const int32_t unilateral_clip_reads = std::max(
        std::max(0, summary.alt_left_clip_reads),
        std::max(0, summary.alt_right_clip_reads));
    const int32_t ref = std::max(0, summary.ref_span_reads);
    const bool low_reference_conflict =
        summary.alt_struct_reads >= (2 * std::max(1, ref));
    return unilateral_clip &&
           unilateral_clip_reads >= 3 &&
           summary.alt_struct_reads >= 3 &&
           low_reference_conflict;
}

bool Pipeline::should_record_hypothesis_in_evidence_ledger(
    const HypothesisSummary& summary) const {
    return summary.alt_struct_reads > 0 ||
           summary.alt_split_reads > 0 ||
           summary.alt_indel_reads > 0 ||
           summary.alt_left_clip_reads > 0 ||
           summary.alt_right_clip_reads > 0 ||
           summary.ref_span_reads > 0;
}

GenotypeCall Pipeline::genotype_call(
    const ComponentCall& component,
    const EventReadEvidence& event_evidence) const {
    GenotypeCall call;
    call.tid = component.tid;
    call.pos = component.anchor_pos;

    EventGenotypeInput input;
    input.alt_struct_reads = event_evidence.alt_struct_reads;
    input.alt_split_reads = event_evidence.alt_split_reads;
    input.alt_indel_reads = event_evidence.alt_indel_reads;
    input.alt_left_clip_reads = event_evidence.alt_left_clip_reads;
    input.alt_right_clip_reads = event_evidence.alt_right_clip_reads;
    input.ref_span_reads = event_evidence.ref_span_reads;
    input.min_depth = config_.genotype_min_depth;
    input.error_rate = config_.genotype_error_rate;
    input.overdispersion = config_.genotype_overdispersion;
    input.alt_observed_lengths = collect_alt_observed_lengths(component, event_evidence);
    input.event_length = infer_event_length_from_alt_support(input.alt_observed_lengths);
    const EventGenotypeDecision decision = genotype_event_from_alt_vs_ref(input);

    call.genotype = decision.best_gt;
    call.af = decision.allele_fraction;
    call.gq = decision.gq;
    return call;
}

BoundaryEvidence Pipeline::analyze_boundary(
    const EventSegmentation& event_segmentation,
    int32_t breakpoint_envelope_width) const {
    FinalBoundaryInput boundary_input;
    boundary_input.left_ref_start = event_segmentation.left_ref_start;
    boundary_input.left_ref_end = event_segmentation.left_ref_end;
    boundary_input.right_ref_start = event_segmentation.right_ref_start;
    boundary_input.right_ref_end = event_segmentation.right_ref_end;
    boundary_input.tsd_min_len = std::max(1, config_.tsd_min_len);
    boundary_input.tsd_max_len = std::max(boundary_input.tsd_min_len, config_.tsd_max_len);
    return evaluate_boundary_evidence(boundary_input, breakpoint_envelope_width);
}

namespace {

struct BreakpointPosteriorSummary {
    double ci_width = 0.0;
    double entropy = 0.0;
};

// L1: breakpoint-position posterior as a product/mixture of per-signal kernels.
// Each supporting signal is a noisy measurement of the latent breakpoint: a
// split/SA supplementary position is precise (small sigma), a CIGAR insertion
// position slightly less so, a soft-clip position least. The posterior over
// position is the (normalised) sum of Gaussian kernels; from it we report the
// 90% credible-interval width and the normalised entropy -- an honest breakpoint
// uncertainty rather than a single coordinate.
BreakpointPosteriorSummary compute_breakpoint_position_posterior(
    const std::vector<BreakpointCandidate>& candidates) {
    BreakpointPosteriorSummary summary;
    std::vector<std::pair<int32_t, double>> signals;
    signals.reserve(candidates.size());
    int32_t lo = std::numeric_limits<int32_t>::max();
    int32_t hi = std::numeric_limits<int32_t>::min();
    for (const auto& bp : candidates) {
        if (bp.pos < 0) {
            continue;
        }
        double sigma = 8.0;
        if ((bp.class_mask & kCandidateSplitSaSupplementary) != 0) {
            sigma = 2.0;
        } else if ((bp.class_mask & kCandidateLongInsertion) != 0) {
            sigma = 3.0;
        } else if ((bp.class_mask & kCandidateSoftClip) != 0) {
            sigma = 6.0;
        }
        signals.push_back({bp.pos, sigma});
        lo = std::min(lo, bp.pos);
        hi = std::max(hi, bp.pos);
    }
    if (signals.size() < 2) {
        return summary;  // not enough signals to form a posterior
    }

    constexpr int32_t kPad = 20;
    constexpr int64_t kMaxSpan = 2000;
    int64_t grid_lo = static_cast<int64_t>(lo) - kPad;
    int64_t grid_hi = static_cast<int64_t>(hi) + kPad;
    if (grid_hi - grid_lo > kMaxSpan) {
        grid_hi = grid_lo + kMaxSpan;
    }
    const int64_t n = grid_hi - grid_lo + 1;
    std::vector<double> dens(static_cast<size_t>(n), 0.0);
    double total = 0.0;
    for (int64_t gi = 0; gi < n; ++gi) {
        const double x = static_cast<double>(grid_lo + gi);
        double f = 0.0;
        for (const auto& signal : signals) {
            const double z = (x - static_cast<double>(signal.first)) / signal.second;
            f += std::exp(-0.5 * z * z) / signal.second;
        }
        dens[static_cast<size_t>(gi)] = f;
        total += f;
    }
    if (total <= 0.0) {
        return summary;
    }

    double entropy = 0.0;
    for (double& d : dens) {
        d /= total;
        if (d > 0.0) {
            entropy -= d * std::log(d);
        }
    }
    summary.entropy = entropy / std::log(static_cast<double>(n));

    double cum = 0.0;
    int64_t q05 = 0;
    int64_t q95 = n - 1;
    for (int64_t gi = 0; gi < n; ++gi) {
        cum += dens[static_cast<size_t>(gi)];
        if (cum >= 0.05) {
            q05 = gi;
            break;
        }
    }
    cum = 0.0;
    for (int64_t gi = 0; gi < n; ++gi) {
        cum += dens[static_cast<size_t>(gi)];
        if (cum >= 0.95) {
            q95 = gi;
            break;
        }
    }
    summary.ci_width = static_cast<double>(q95 - q05);
    return summary;
}

}  // namespace

FinalCall Pipeline::emit_final_te_call(
    const ComponentCall& component,
    const EventReadEvidence& event_evidence,
    const EventConsensus& event_consensus,
    const EventSegmentation& event_segmentation,
    const TEAlignmentEvidence& te_alignment,
    const GenotypeCall& genotype,
    const FinalBoundaryDecision& boundary,
    const JointDecisionResult& joint) const {
    FinalCall call;
    call.chrom = component.chrom;
    call.tid = component.tid;
    const bool one_sided_left =
        event_segmentation.pass &&
        !event_segmentation.insert_seq.empty() &&
        event_segmentation.left_flank_align_len > 0 &&
        event_segmentation.right_flank_align_len == 0 &&
        event_segmentation.left_ref_end >= 0;
    const bool one_sided_right =
        event_segmentation.pass &&
        !event_segmentation.insert_seq.empty() &&
        event_segmentation.right_flank_align_len > 0 &&
        event_segmentation.left_flank_align_len == 0 &&
        event_segmentation.right_ref_start >= 0;
    const bool one_sided_flank_anchored = one_sided_left || one_sided_right;
    const bool evidence_tier_call =
        joint.emit_structural_event_call ||
        joint.emit_evidence_te_call ||
        joint.final_qc == "PASS_TE_IMPRECISE";
    if (one_sided_flank_anchored) {
        const int32_t anchored_pos = one_sided_left
            ? event_segmentation.left_ref_end
            : event_segmentation.right_ref_start;
        call.bp_left = anchored_pos;
        call.bp_right = anchored_pos;
    } else if (evidence_tier_call &&
        event_evidence.bp_left >= 0 &&
        event_evidence.bp_right >= 0) {
        call.bp_left = std::min(event_evidence.bp_left, event_evidence.bp_right);
        call.bp_right = std::max(event_evidence.bp_left, event_evidence.bp_right);
    } else {
        call.bp_left = std::min(
            event_segmentation.left_ref_end,
            event_segmentation.right_ref_start);
        call.bp_right = std::max(
            event_segmentation.left_ref_end,
            event_segmentation.right_ref_start);
    }
    call.pos = (call.bp_left >= 0 && call.bp_right >= 0)
        ? (call.bp_left + ((call.bp_right - call.bp_left) / 2))
        : component.anchor_pos;
    call.window_start = component.bin_start;
    call.window_end = component.bin_end;

    const BreakpointPosteriorSummary bp_posterior =
        compute_breakpoint_position_posterior(component.breakpoint_candidates);
    call.bp_ci_width = bp_posterior.ci_width;
    call.bp_posterior_entropy = bp_posterior.entropy;

    const bool emit_unknown_te =
        joint.emit_unknown_te || joint.emit_structural_event_call;
    const TeSequenceExplanation& sequence_explanation =
        te_alignment.te_sequence_explanation;
    const bool sequence_family_certificate =
        emit_unknown_te &&
        !joint.emit_structural_event_call &&
        joint.emit_te_call &&
        te_alignment.pass &&
        te_alignment.qc_reason == "PASS_INSERT_TE_ALIGNMENT" &&
        te_alignment.annotation_confidence == "HIGH" &&
        !te_alignment.best_family.empty() &&
        !te_alignment.best_subfamily.empty() &&
        te_alignment.cross_family_margin >= config_.te_family_margin_min &&
        sequence_explanation.status == TeAnnotationStatus::kResolved &&
        sequence_explanation.family == te_alignment.best_family &&
        sequence_explanation.subfamily == te_alignment.best_subfamily &&
        joint.te_structure_path_confidence >= (1.0 - config_.final_fdr_q);
    call.family = emit_unknown_te
        ? "UNKNOWN"
        : (te_alignment.best_family.empty() ? "NA" : te_alignment.best_family);
    call.subfamily = emit_unknown_te
        ? "UNKNOWN"
        : (te_alignment.best_subfamily.empty() ? "NA" : te_alignment.best_subfamily);
    call.family_committed =
        !emit_unknown_te &&
        te_alignment.pass &&
        (te_alignment.qc_reason == "PASS_INSERT_TE_ALIGNMENT" ||
         te_alignment.qc_reason == "PASS_INSERT_TE_ALIGNMENT_FAMILY_ONLY") &&
        !te_alignment.best_family.empty();
    if (sequence_family_certificate) {
        call.sequence_family_candidate = te_alignment.best_family;
        call.sequence_subfamily_candidate = te_alignment.best_subfamily;
        call.sequence_family_commit_eligible = true;
    }
    call.te_name = (call.subfamily != "NA") ? call.subfamily : call.family;
    call.strand = "NA";
    call.insert_len = static_cast<int32_t>(event_segmentation.insert_seq.size());
    call.best_te_identity = te_alignment.best_identity;
    call.te_consensus_start = te_alignment.te_consensus_start;
    call.te_consensus_end = te_alignment.te_consensus_end;
    call.best_te_query_coverage = te_alignment.best_query_coverage;
    call.cross_family_margin = te_alignment.cross_family_margin;
    call.te_sequence_model_label = te_alignment.sequence_model_label;
    call.te_sequence_model_score = te_alignment.sequence_model_score;
    call.te_sequence_model_gc = te_alignment.sequence_model_gc;
    call.te_sequence_model_entropy = te_alignment.sequence_model_entropy;
    call.te_sequence_model_tandem_fraction = te_alignment.sequence_model_tandem_fraction;
    call.te_sequence_model_low_complexity_fraction =
        te_alignment.sequence_model_low_complexity_fraction;
    call.te_sequence_model_jsd_k5 = te_alignment.sequence_model_jsd_k5;
    call.te_sequence_model_jsd_k6 = te_alignment.sequence_model_jsd_k6;
    call.te_sequence_model_k9_containment = te_alignment.sequence_model_k9_containment;
    call.te_annotation_confidence = te_alignment.annotation_confidence;
    call.te_annotation_class = te_alignment.annotation_class;
    call.te_annotation_order = te_alignment.annotation_order;
    call.te_annotation_intervals = te_alignment.annotation_intervals;
    call.te_annotation_residual_fraction = te_alignment.annotation_residual_fraction;
    call.te_annotation_masked_fraction = te_alignment.annotation_masked_fraction;
    call.left_flank_align_len = event_segmentation.left_flank_align_len;
    call.right_flank_align_len = event_segmentation.right_flank_align_len;
    call.event_consensus_len = event_consensus.consensus_len;
    call.insert_seq = event_segmentation.insert_seq;

    call.te_qc = te_alignment.qc_reason;
    call.final_qc = joint.final_qc;
    call.best_explanation = joint.best_explanation;
    call.explanation_residual = joint.explanation_residual;
    call.explanation_path = joint.explanation_path;
    call.te_structure_path = joint.te_structure_path;
    call.te_structure_log_evidence = joint.te_structure_log_evidence;
    call.nonte_structure_log_evidence = joint.nonte_structure_log_evidence;
    call.artifact_structure_log_evidence = joint.artifact_structure_log_evidence;
    call.te_structure_path_confidence = joint.te_structure_path_confidence;
    call.polyA_posterior = joint.polyA_posterior;
    call.transduction_posterior = joint.transduction_posterior;
    call.te_core_coverage = joint.te_core_coverage;
    call.unexplained_high_complexity_bp = joint.unexplained_high_complexity_bp;
    call.te_posterior = joint.te_posterior;
    call.non_te_posterior = joint.non_te_posterior;
    call.artifact_posterior = joint.artifact_posterior;
    call.te_vs_artifact_log_odds = joint.te_vs_artifact_log_odds;
    call.te_vs_non_te_log_odds = joint.te_vs_non_te_log_odds;
    call.posterior_qc = joint.posterior_qc;
    call.latent_mechanism = joint.latent_mechanism;
    call.family_activity_prior = joint.family_activity_prior;
    call.lfdr = joint.lfdr;
    call.worst_case_lfdr = joint.worst_case_lfdr;
    call.lfdr_qc = joint.lfdr_qc;
    call.mechanistic_lower_log_bf_te_vs_artifact =
        joint.mechanistic_lower_log_bf_te_vs_artifact;
    call.mechanistic_lower_log_bf_te_vs_non_te =
        joint.mechanistic_lower_log_bf_te_vs_non_te;
    call.mechanistic_ref_conflict_signal =
        joint.mechanistic_ref_conflict_signal;
    call.mechanistic_ambiguity_width =
        joint.mechanistic_ambiguity_width;
    call.mechanistic_blocks = joint.mechanistic_blocks;
    call.robust_mechanistic_lfdr = joint.robust_mechanistic_lfdr;
    call.robust_mechanistic_worst_case_lfdr =
        joint.robust_mechanistic_worst_case_lfdr;
    call.robust_mechanistic_qc = joint.robust_mechanistic_qc;
    append_qc_token(call.final_qc, te_alignment.qc_reason.c_str());
    append_qc_token(call.final_qc, boundary.qc.c_str());
    append_qc_token(call.final_qc, event_consensus.qc_reason.c_str());
    append_qc_token(call.final_qc, event_segmentation.qc_reason.c_str());
    if (one_sided_flank_anchored) {
        append_qc_token(call.final_qc, "CONSENSUS_FLANK_ANCHORED_BREAKPOINT");
    }

    call.tsd_type = boundary.boundary_type;
    call.tsd_len = boundary.boundary_len;
    call.tsd_seq = "NA";
    call.tsd_bg_p = 1.0;

    call.support_reads = event_evidence.alt_struct_reads;
    call.alt_struct_reads = event_evidence.alt_struct_reads;
    call.raw_cigar_insert_reads = event_evidence.raw_cigar_insert_reads;
    call.max_raw_cigar_insert_len = event_evidence.max_raw_cigar_insert_len;
    call.ref_span_reads = event_evidence.ref_span_reads;
    call.low_mapq_ref_span_reads = event_evidence.low_mapq_ref_span_reads;
    call.support_qnames = event_evidence.support_qnames;

    call.genotype = genotype.genotype;
    call.af = genotype.af;
    call.gq = genotype.gq;
    return call;
}
