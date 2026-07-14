#ifdef NDEBUG
#undef NDEBUG
#endif
#include "pipeline.h"

#include <cassert>
#include <vector>

int main() {
    using namespace placer;

    {
        std::vector<ComponentFinalCallCandidate> candidates;

        ComponentFinalCallCandidate truth_near;
        truth_near.pos = 17680845;
        truth_near.score = 4.19922;
        truth_near.emit_te = true;
        candidates.push_back(truth_near);

        ComponentFinalCallCandidate distant_resolved;
        distant_resolved.pos = 17682680;
        distant_resolved.score = 5.95314;
        distant_resolved.emit_te = true;
        candidates.push_back(distant_resolved);

        const auto selected = select_component_final_call_indices(candidates);
        assert(selected.size() == 2);
        assert(selected[0] == 0);
        assert(selected[1] == 1);
    }

    {
        std::vector<ComponentFinalCallCandidate> candidates;

        ComponentFinalCallCandidate weaker;
        weaker.pos = 1000;
        weaker.score = 2.0;
        weaker.emit_te = true;
        candidates.push_back(weaker);

        ComponentFinalCallCandidate stronger_same_locus;
        stronger_same_locus.pos = 1040;
        stronger_same_locus.score = 3.0;
        stronger_same_locus.emit_te = true;
        candidates.push_back(stronger_same_locus);

        ComponentFinalCallCandidate rejected;
        rejected.pos = 1400;
        rejected.score = 9.0;
        rejected.emit_te = false;
        candidates.push_back(rejected);

        const auto selected = select_component_final_call_indices(candidates);
        assert(selected.size() == 1);
        assert(selected[0] == 1);
    }

    {
        std::vector<ComponentFinalCallCandidate> candidates;

        ComponentFinalCallCandidate left_event;
        left_event.pos = 1000;
        left_event.score = 5.0;
        left_event.emit_te = true;
        candidates.push_back(left_event);

        ComponentFinalCallCandidate anchor_only;
        anchor_only.pos = 2000;
        anchor_only.score = 100.0;
        anchor_only.emit_te = false;
        candidates.push_back(anchor_only);

        ComponentFinalCallCandidate right_event;
        right_event.pos = 3000;
        right_event.score = 4.0;
        right_event.emit_te = true;
        candidates.push_back(right_event);

        const auto selected = select_component_final_call_indices(candidates);
        assert(selected.size() == 2);
        assert(selected[0] == 0);
        assert(selected[1] == 2);
    }

    {
        std::vector<ComponentFinalCallCandidate> candidates;

        ComponentFinalCallCandidate left_truth;
        left_truth.pos = 1000;
        left_truth.score = 2.0;
        left_truth.emit_te = true;
        candidates.push_back(left_truth);

        ComponentFinalCallCandidate adjacent_truth;
        adjacent_truth.pos = 1070;
        adjacent_truth.score = 3.0;
        adjacent_truth.emit_te = true;
        candidates.push_back(adjacent_truth);

        const auto selected = select_component_final_call_indices(candidates);
        assert(selected.size() == 2);
        assert(selected[0] == 0);
        assert(selected[1] == 1);
    }

    {
        std::vector<ComponentFinalCallCandidate> candidates;

        ComponentFinalCallCandidate te_evidence;
        te_evidence.pos = 28237107;
        te_evidence.anchor_pos = 28235587;
        te_evidence.score = 4.39;
        te_evidence.emit_te = true;
        te_evidence.evidence_te = true;
        te_evidence.one_sided_segmentation = true;
        te_evidence.anchor_support = 6;
        te_evidence.anchor_ref_span_reads = 8;
        candidates.push_back(te_evidence);

        retether_evidence_supported_final_call_positions(candidates);

        assert(candidates[0].pos == 28237107);
    }

    {
        std::vector<ComponentFinalCallCandidate> candidates;

        ComponentFinalCallCandidate high_confidence;
        high_confidence.pos = 38193873;
        high_confidence.anchor_pos = 38195168;
        high_confidence.score = 4.37;
        high_confidence.emit_te = true;
        high_confidence.evidence_te = false;
        high_confidence.one_sided_segmentation = false;
        high_confidence.anchor_support = 12;
        candidates.push_back(high_confidence);

        retether_evidence_supported_final_call_positions(candidates);

        assert(candidates[0].pos == 38193873);
    }

    {
        std::vector<ComponentFinalCallCandidate> candidates;

        ComponentFinalCallCandidate distant_indel_anchor;
        distant_indel_anchor.pos = 28234790;
        distant_indel_anchor.emit_te = false;
        distant_indel_anchor.anchor_support = 3;
        distant_indel_anchor.anchor_priority = 1;
        distant_indel_anchor.anchor_hypothesis_score = 24.0;
        distant_indel_anchor.component_anchor_pos = 28236188;
        candidates.push_back(distant_indel_anchor);

        ComponentFinalCallCandidate truth_near_indel_anchor;
        truth_near_indel_anchor.pos = 28235587;
        truth_near_indel_anchor.emit_te = false;
        truth_near_indel_anchor.anchor_support = 1;
        truth_near_indel_anchor.anchor_priority = 1;
        truth_near_indel_anchor.anchor_hypothesis_score = 8.0;
        truth_near_indel_anchor.component_anchor_pos = 28236188;
        candidates.push_back(truth_near_indel_anchor);

        ComponentFinalCallCandidate component_proximal_clip_anchor;
        component_proximal_clip_anchor.pos = 28236410;
        component_proximal_clip_anchor.emit_te = false;
        component_proximal_clip_anchor.anchor_support = 8;
        component_proximal_clip_anchor.anchor_priority = 3;
        component_proximal_clip_anchor.anchor_hypothesis_score = 8.0;
        component_proximal_clip_anchor.component_anchor_pos = 28236188;
        candidates.push_back(component_proximal_clip_anchor);

        ComponentFinalCallCandidate te_evidence;
        te_evidence.pos = 28237107;
        te_evidence.anchor_pos = 28237127;
        te_evidence.score = 4.39;
        te_evidence.emit_te = true;
        te_evidence.evidence_te = true;
        te_evidence.one_sided_segmentation = true;
        te_evidence.anchor_support = 2;
        te_evidence.anchor_ref_span_reads = 7;
        te_evidence.anchor_priority = 2;
        te_evidence.anchor_hypothesis_score = 14.0;
        te_evidence.component_anchor_pos = 28236188;
        candidates.push_back(te_evidence);

        retether_evidence_supported_final_call_positions(candidates);

        assert(candidates[3].pos == 28237107);
    }

    {
        std::vector<ComponentFinalCallCandidate> candidates;

        ComponentFinalCallCandidate far_strong_anchor;
        far_strong_anchor.pos = 17680845;
        far_strong_anchor.emit_te = false;
        far_strong_anchor.anchor_support = 6;
        far_strong_anchor.anchor_priority = 1;
        far_strong_anchor.anchor_hypothesis_score = 48.0;
        far_strong_anchor.component_anchor_pos = 17681616;
        candidates.push_back(far_strong_anchor);

        ComponentFinalCallCandidate component_proximal_clip_anchor;
        component_proximal_clip_anchor.pos = 17681611;
        component_proximal_clip_anchor.emit_te = false;
        component_proximal_clip_anchor.anchor_support = 7;
        component_proximal_clip_anchor.anchor_priority = 3;
        component_proximal_clip_anchor.anchor_hypothesis_score = 7.0;
        component_proximal_clip_anchor.component_anchor_pos = 17681616;
        candidates.push_back(component_proximal_clip_anchor);

        ComponentFinalCallCandidate te_evidence;
        te_evidence.pos = 17681850;
        te_evidence.anchor_pos = 17681850;
        te_evidence.score = 4.14;
        te_evidence.emit_te = true;
        te_evidence.evidence_te = true;
        te_evidence.one_sided_segmentation = true;
        te_evidence.anchor_support = 1;
        te_evidence.anchor_priority = 0;
        te_evidence.anchor_hypothesis_score = 8.0;
        te_evidence.component_anchor_pos = 17681616;
        candidates.push_back(te_evidence);

        retether_evidence_supported_final_call_positions(candidates);

        assert(candidates[2].pos == 17681850);
    }

    {
        std::vector<ComponentFinalCallCandidate> candidates;

        ComponentFinalCallCandidate component_proximal_anchor;
        component_proximal_anchor.pos = 10556971;
        component_proximal_anchor.emit_te = false;
        component_proximal_anchor.anchor_support = 8;
        component_proximal_anchor.anchor_priority = 3;
        component_proximal_anchor.anchor_hypothesis_score = 8.0;
        component_proximal_anchor.component_anchor_pos = 10556989;
        candidates.push_back(component_proximal_anchor);

        ComponentFinalCallCandidate resolved_te_evidence;
        resolved_te_evidence.pos = 10557659;
        resolved_te_evidence.anchor_pos = 10557659;
        resolved_te_evidence.score = 3.48;
        resolved_te_evidence.emit_te = true;
        resolved_te_evidence.evidence_te = true;
        resolved_te_evidence.resolved_te = true;
        resolved_te_evidence.one_sided_segmentation = true;
        resolved_te_evidence.anchor_support = 2;
        resolved_te_evidence.anchor_ref_span_reads = 8;
        resolved_te_evidence.anchor_priority = 1;
        resolved_te_evidence.anchor_hypothesis_score = 16.0;
        resolved_te_evidence.component_anchor_pos = 10556989;
        resolved_te_evidence.retether_direction = -1;
        candidates.push_back(resolved_te_evidence);

        retether_evidence_supported_final_call_positions(candidates);

        assert(candidates[1].pos == 10557659);
    }

    {
        std::vector<ComponentFinalCallCandidate> candidates;

        ComponentFinalCallCandidate current_position_anchor;
        current_position_anchor.pos = 17261244;
        current_position_anchor.emit_te = false;
        current_position_anchor.anchor_support = 1;
        current_position_anchor.anchor_priority = 2;
        current_position_anchor.anchor_hypothesis_score = 7.0;
        current_position_anchor.component_anchor_pos = 17261846;
        candidates.push_back(current_position_anchor);

        ComponentFinalCallCandidate local_precise_anchor;
        local_precise_anchor.pos = 17261978;
        local_precise_anchor.emit_te = false;
        local_precise_anchor.anchor_support = 2;
        local_precise_anchor.anchor_priority = 2;
        local_precise_anchor.anchor_hypothesis_score = 14.0;
        local_precise_anchor.component_anchor_pos = 17261846;
        candidates.push_back(local_precise_anchor);

        ComponentFinalCallCandidate unknown_te_evidence;
        unknown_te_evidence.pos = 17261846;
        unknown_te_evidence.anchor_pos = 17261978;
        unknown_te_evidence.score = 2.98;
        unknown_te_evidence.emit_te = true;
        unknown_te_evidence.evidence_te = true;
        unknown_te_evidence.one_sided_segmentation = true;
        unknown_te_evidence.anchor_support = 2;
        unknown_te_evidence.anchor_ref_span_reads = 17;
        unknown_te_evidence.anchor_priority = 5;
        unknown_te_evidence.anchor_hypothesis_score = 10.0;
        unknown_te_evidence.component_anchor_pos = 17261846;
        unknown_te_evidence.retether_direction = -1;
        candidates.push_back(unknown_te_evidence);

        retether_evidence_supported_final_call_positions(candidates);

        assert(candidates[2].pos == 17261846);
    }

    {
        std::vector<ComponentFinalCallCandidate> candidates;

        ComponentFinalCallCandidate left_precise_anchor;
        left_precise_anchor.pos = 17261244;
        left_precise_anchor.emit_te = false;
        left_precise_anchor.anchor_support = 1;
        left_precise_anchor.anchor_priority = 2;
        left_precise_anchor.anchor_hypothesis_score = 7.0;
        left_precise_anchor.component_anchor_pos = 17261846;
        candidates.push_back(left_precise_anchor);

        ComponentFinalCallCandidate local_precise_anchor;
        local_precise_anchor.pos = 17261978;
        local_precise_anchor.emit_te = false;
        local_precise_anchor.anchor_support = 2;
        local_precise_anchor.anchor_priority = 2;
        local_precise_anchor.anchor_hypothesis_score = 14.0;
        local_precise_anchor.component_anchor_pos = 17261846;
        candidates.push_back(local_precise_anchor);

        ComponentFinalCallCandidate unknown_te_evidence;
        unknown_te_evidence.pos = 17261978;
        unknown_te_evidence.anchor_pos = 17261978;
        unknown_te_evidence.score = 2.98;
        unknown_te_evidence.emit_te = true;
        unknown_te_evidence.evidence_te = true;
        unknown_te_evidence.one_sided_segmentation = true;
        unknown_te_evidence.anchor_support = 2;
        unknown_te_evidence.anchor_ref_span_reads = 17;
        unknown_te_evidence.anchor_priority = 2;
        unknown_te_evidence.anchor_hypothesis_score = 14.0;
        unknown_te_evidence.component_anchor_pos = 17261846;
        unknown_te_evidence.retether_direction = -1;
        candidates.push_back(unknown_te_evidence);

        retether_evidence_supported_final_call_positions(candidates);

        assert(candidates[2].pos == 17261978);
    }

    {
        std::vector<ComponentFinalCallCandidate> candidates;

        ComponentFinalCallCandidate precise_anchor;
        precise_anchor.pos = 9500393;
        precise_anchor.emit_te = false;
        precise_anchor.anchor_support = 1;
        precise_anchor.anchor_priority = 1;
        precise_anchor.anchor_hypothesis_score = 8.0;
        precise_anchor.component_anchor_pos = 9501801;
        candidates.push_back(precise_anchor);

        ComponentFinalCallCandidate complete_te_evidence;
        complete_te_evidence.pos = 9499393;
        complete_te_evidence.anchor_pos = 9500393;
        complete_te_evidence.score = 3.70;
        complete_te_evidence.emit_te = true;
        complete_te_evidence.evidence_te = true;
        complete_te_evidence.one_sided_segmentation = false;
        complete_te_evidence.anchor_support = 5;
        complete_te_evidence.anchor_ref_span_reads = 4;
        complete_te_evidence.anchor_priority = 1;
        complete_te_evidence.anchor_hypothesis_score = 8.0;
        complete_te_evidence.component_anchor_pos = 9501801;
        candidates.push_back(complete_te_evidence);

        retether_evidence_supported_final_call_positions(candidates);

        assert(candidates[1].pos == 9500393);
    }

    {
        std::vector<ComponentFinalCallCandidate> candidates;

        ComponentFinalCallCandidate component_proximal_anchor;
        component_proximal_anchor.pos = 67144249;
        component_proximal_anchor.emit_te = false;
        component_proximal_anchor.anchor_support = 14;
        component_proximal_anchor.anchor_priority = 5;
        component_proximal_anchor.anchor_hypothesis_score = 14.0;
        component_proximal_anchor.component_anchor_pos = 67145040;
        candidates.push_back(component_proximal_anchor);

        ComponentFinalCallCandidate complete_unknown_evidence;
        complete_unknown_evidence.pos = 67144113;
        complete_unknown_evidence.anchor_pos = 67144113;
        complete_unknown_evidence.score = 4.67;
        complete_unknown_evidence.emit_te = true;
        complete_unknown_evidence.evidence_te = true;
        complete_unknown_evidence.resolved_te = false;
        complete_unknown_evidence.one_sided_segmentation = false;
        complete_unknown_evidence.anchor_support = 10;
        complete_unknown_evidence.anchor_ref_span_reads = 4;
        complete_unknown_evidence.anchor_priority = 3;
        complete_unknown_evidence.anchor_hypothesis_score = 10.0;
        complete_unknown_evidence.component_anchor_pos = 67145040;
        candidates.push_back(complete_unknown_evidence);

        retether_evidence_supported_final_call_positions(candidates);

        assert(candidates[1].pos == 67144113);
    }

    {
        std::vector<ComponentFinalCallCandidate> candidates;

        ComponentFinalCallCandidate truth_precise_anchor;
        truth_precise_anchor.pos = 31965349;
        truth_precise_anchor.emit_te = false;
        truth_precise_anchor.anchor_support = 1;
        truth_precise_anchor.anchor_priority = 4;
        truth_precise_anchor.anchor_hypothesis_score = 6.0;
        truth_precise_anchor.component_anchor_pos = 31965756;
        candidates.push_back(truth_precise_anchor);

        ComponentFinalCallCandidate complete_unknown_evidence;
        complete_unknown_evidence.pos = 31965688;
        complete_unknown_evidence.anchor_pos = 31965688;
        complete_unknown_evidence.score = 4.00;
        complete_unknown_evidence.emit_te = true;
        complete_unknown_evidence.evidence_te = true;
        complete_unknown_evidence.resolved_te = false;
        complete_unknown_evidence.one_sided_segmentation = false;
        complete_unknown_evidence.anchor_support = 7;
        complete_unknown_evidence.anchor_ref_span_reads = 6;
        complete_unknown_evidence.anchor_priority = 5;
        complete_unknown_evidence.anchor_hypothesis_score = 7.0;
        complete_unknown_evidence.component_anchor_pos = 31965756;
        candidates.push_back(complete_unknown_evidence);

        retether_evidence_supported_final_call_positions(candidates);

        assert(candidates[1].pos == 31965349);
    }

    return 0;
}
