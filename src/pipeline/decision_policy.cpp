#include "decision_policy.h"
#include "pipeline.h"

#include <algorithm>
#include <cmath>
#include <unordered_set>
#include <vector>

namespace placer {

namespace {

int32_t median_i32(std::vector<int32_t> values) {
    if (values.empty()) {
        return -1;
    }
    std::sort(values.begin(), values.end());
    const size_t mid = values.size() / 2;
    if (values.size() % 2 == 1) {
        return values[mid];
    }
    return static_cast<int32_t>(
        (static_cast<int64_t>(values[mid - 1]) + static_cast<int64_t>(values[mid])) / 2);
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

}  // namespace

InsertionAcceptanceDecision evaluate_insertion_acceptance(
    const InsertionEvidence& evidence,
    const InsertionAcceptanceParams& params) {
    InsertionAcceptanceDecision decision;

    if (evidence.tier <= 2) {
        decision.pass = true;
        return decision;
    }

    if (!params.emit_low_confidence_calls) {
        return decision;
    }

    const int32_t min_support = std::max(1, params.low_conf_min_support_reads);
    const int32_t max_tier = std::clamp(params.low_conf_max_tier, 1, 3);
    const bool low_conf_ok =
        evidence.support_reads >= min_support &&
        evidence.tier <= max_tier;
    if (low_conf_ok) {
        decision.pass = true;
        decision.low_confidence = true;
    }
    return decision;
}

TeOpenSetDecision classify_te_open_set(
    const TeOpenSetInput& input,
    const TeOpenSetParams& params) {
    TeOpenSetDecision decision;

    const double conf_certain_min = std::clamp(params.conf_certain_min, 0.0, 1.0);
    const double conf_uncertain_min = std::clamp(
        std::min(params.conf_uncertain_min, conf_certain_min),
        0.0,
        1.0);
    const bool prob_certain = input.te_confidence_prob >= conf_certain_min;
    const bool prob_uncertain = input.te_confidence_prob >= conf_uncertain_min;

    const bool uncertain_path = !input.te_gate_pass || input.te_gate_uncertain_path;
    if (!uncertain_path) {
        decision.status = "TE_CERTAIN";
        return decision;
    }

    if (!input.te_gate_pass) {
        decision.add_te_gate_fail_tag = true;
    }
    if (input.has_proxy_signal && (input.proxy_certain || prob_certain)) {
        decision.status = "TE_CERTAIN";
        decision.promote_top1_name = true;
        decision.qc_proxy_tag = "TE_PROXY_CERTAIN";
    } else if (input.has_proxy_signal && prob_uncertain) {
        decision.status = "TE_UNCERTAIN";
        decision.qc_proxy_tag = "TE_PROXY_WEAK";
    } else {
        decision.status = "NON_TE";
        decision.qc_proxy_tag = "TE_PROXY_NONE";
    }

    return decision;
}

PostAssemblyTeDecision evaluate_post_assembly_te_decision(
    const ClusterTECall& te_call,
    const AssemblyCall& assembly,
    const ComponentCall& component,
    const PipelineConfig& config) {
    PostAssemblyTeDecision decision;

    std::unordered_set<size_t> non_softclip_support;
    non_softclip_support.insert(
        component.split_sa_read_indices.begin(),
        component.split_sa_read_indices.end());
    non_softclip_support.insert(
        component.insertion_read_indices.begin(),
        component.insertion_read_indices.end());
    const int32_t non_softclip_reads = static_cast<int32_t>(non_softclip_support.size());
    const int32_t softclip_reads = static_cast<int32_t>(component.soft_clip_read_indices.size());

    std::unordered_set<size_t> all_support = non_softclip_support;
    all_support.insert(
        component.soft_clip_read_indices.begin(),
        component.soft_clip_read_indices.end());
    const int32_t total_support_reads = static_cast<int32_t>(all_support.size());

    std::vector<int32_t> bp_positions;
    bp_positions.reserve(component.breakpoint_candidates.size());
    for (const auto& bp : component.breakpoint_candidates) {
        if (bp.pos >= 0) {
            bp_positions.push_back(bp.pos);
        }
    }
    const double bp_mad = mad_from_positions(bp_positions);
    const bool one_sided_geom_ok =
        bp_positions.size() < 2 ||
        bp_mad <= std::max(0.0, config.te_one_sided_breakpoint_mad_max);

    const bool no_softclip_component =
        softclip_reads <= 0 && non_softclip_reads > 0;
    const bool split_indel_dominant =
        non_softclip_reads > 0 && non_softclip_reads >= softclip_reads;

    const int32_t short_ins_len_min = std::max(1, config.short_ins_min_len);
    const int32_t short_ins_len_max = std::max(short_ins_len_min, config.short_ins_max_len);
    const bool short_ins_candidate =
        config.short_ins_enable &&
        assembly.consensus_len >= short_ins_len_min &&
        assembly.consensus_len <= short_ins_len_max &&
        total_support_reads >= std::max(1, config.short_ins_min_reads);
    const double short_ins_identity_min = std::clamp(
        std::max(0.0, config.assembly_min_identity_est - config.short_ins_kmer_relax_identity),
        0.0,
        1.0);

    if (te_call.te_name.empty()) {
        const int32_t min_non_softclip_reads = std::max(
            2,
            no_softclip_component ? config.te_no_softclip_min_reads : config.te_min_fragments_for_vote);

        if (assembly.qc_pass &&
            assembly.identity_est >= std::clamp(config.assembly_min_identity_est, 0.0, 1.0) &&
            non_softclip_reads >= min_non_softclip_reads) {
            decision.pass = true;
            if (no_softclip_component &&
                one_sided_geom_ok &&
                non_softclip_reads >= std::max(1, config.te_no_softclip_min_reads) &&
                assembly.identity_est >= std::clamp(config.te_no_softclip_identity_min, 0.0, 1.0)) {
                decision.qc = "PASS_SPLIT_INDEL_RESCUE";
            } else {
                decision.qc = "PASS_INSERTION_TE_UNCERTAIN";
            }
        } else if (short_ins_candidate &&
                   split_indel_dominant &&
                   assembly.qc_pass &&
                   assembly.identity_est >= short_ins_identity_min &&
                   one_sided_geom_ok) {
            decision.pass = true;
            decision.qc = "PASS_SHORT_INS_RESCUE";
        } else {
            if (no_softclip_component && !one_sided_geom_ok) {
                decision.qc = "FAIL_SPLIT_INDEL_INCONSISTENT";
                return decision;
            }
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
                const bool split_indel_rescue =
                    split_indel_dominant &&
                    one_sided_geom_ok &&
                    assembly.qc_pass &&
                    non_softclip_reads >= std::max(1, config.te_no_softclip_min_reads) &&
                    te_call.fragment_count >= std::max(1, config.te_no_softclip_min_fragments) &&
                    assembly.identity_est >= std::clamp(config.te_no_softclip_identity_min, 0.0, 1.0);
                if (split_indel_rescue) {
                    decision.pass = true;
                    decision.qc = "PASS_SPLIT_INDEL_RESCUE";
                    return decision;
                }

                const bool short_ins_rescue =
                    short_ins_candidate &&
                    split_indel_dominant &&
                    one_sided_geom_ok &&
                    assembly.qc_pass &&
                    assembly.identity_est >= short_ins_identity_min &&
                    (te_call.median_identity + std::clamp(config.short_ins_kmer_relax_identity, 0.0, 1.0)) >= rescue_identity_min;
                if (short_ins_rescue) {
                    decision.pass = true;
                    decision.qc = "PASS_SHORT_INS_RESCUE";
                    return decision;
                }

                if (no_softclip_component && !one_sided_geom_ok) {
                    decision.qc = "FAIL_SPLIT_INDEL_INCONSISTENT";
                    return decision;
                }
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

}  // namespace placer
