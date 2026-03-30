#include "decision_policy.h"
#include "pipeline.h"

#include <cassert>
#include <string>
#include <vector>

int main() {
    using namespace placer;

    {
        auto make_call = [](
                             int32_t pos,
                             const std::string& te_name,
                             int32_t support_reads,
                             double cross_family_margin,
                             int32_t gq) {
            FinalCall call;
            call.chrom = "1";
            call.tid = 0;
            call.pos = pos;
            call.te_name = te_name;
            call.support_reads = support_reads;
            call.cross_family_margin = cross_family_margin;
            call.gq = gq;
            return call;
        };

        PipelineResult result;
        result.final_calls.push_back(make_call(1000, "ALU:Ya5", 8, 0.20, 18));
        result.final_calls.push_back(make_call(1040, "L1:L1HS", 12, 0.45, 27));
        finalize_final_calls(result);
        assert(result.final_calls.size() == 1);
        assert(result.final_calls.front().pos == 1040);
        assert(result.final_calls.front().te_name == "L1:L1HS");
        assert(result.final_pass_calls == 1);
    }

    {
        auto make_call = [](
                             int32_t pos,
                             int32_t support_reads,
                             int32_t gq,
                             double cross_family_margin) {
            FinalCall call;
            call.chrom = "1";
            call.tid = 0;
            call.pos = pos;
            call.te_name = "ALU:Ya5";
            call.support_reads = support_reads;
            call.gq = gq;
            call.cross_family_margin = cross_family_margin;
            return call;
        };

        PipelineResult result;
        result.final_calls.push_back(make_call(1000, 5, 12, 0.30));
        result.final_calls.push_back(make_call(1080, 9, 20, 0.55));
        result.final_calls.push_back(make_call(1185, 7, 16, 0.40));
        finalize_final_calls(result);
        assert(result.final_calls.size() == 2);
        assert(result.final_calls[0].pos == 1080);
        assert(result.final_calls[1].pos == 1185);
        assert(result.final_pass_calls == 2);
    }

    {
        bam1_t* record = bam_init1();
        assert(record != nullptr);

        const std::string qname = "read_softclip_with_sa";
        const std::string seq(100, 'A');
        const std::string qual(100, 'I');
        const uint32_t cigar[] = {
            bam_cigar_gen(30, BAM_CSOFT_CLIP),
            bam_cigar_gen(70, BAM_CMATCH),
        };
        const int set_ret = bam_set1(
            record,
            qname.size() + 1,
            qname.c_str(),
            0,
            0,
            1000,
            60,
            sizeof(cigar) / sizeof(cigar[0]),
            cigar,
            -1,
            -1,
            0,
            seq.size(),
            seq.c_str(),
            qual.c_str(),
            0);
        assert(set_ret >= 0);

        const char sa_value[] = "1,2000,+,70M30S,60,1;";
        const int aux_ret = bam_aux_append(
            record,
            "SA",
            'Z',
            static_cast<int>(sizeof(sa_value)),
            reinterpret_cast<const uint8_t*>(sa_value));
        assert(aux_ret == 0);

        const std::vector<const bam1_t*> records = {record};
        LinearBinComponentModule module;
        const auto components = module.build(records, "1", 0, 900, 1200);
        assert(components.size() == 1);

        bool saw_softclip_candidate = false;
        int split_hint_count = 0;
        for (const auto& bp : components.front().breakpoint_candidates) {
            if (bp.clip_len > 0) {
                saw_softclip_candidate = true;
                assert((bp.class_mask & kCandidateSoftClip) != 0);
                assert((bp.class_mask & kCandidateSplitSaSupplementary) == 0);
                assert((bp.class_mask & kCandidateLongInsertion) == 0);
            } else if ((bp.class_mask & kCandidateSplitSaSupplementary) != 0) {
                split_hint_count += 1;
                assert(bp.pos == 1000);
            }
        }
        assert(saw_softclip_candidate);
        assert(split_hint_count == 1);

        bam_destroy1(record);
    }

    {
        bam1_t* record = bam_init1();
        assert(record != nullptr);

        const std::string qname = "read_trailing_softclip_with_sa";
        const std::string seq(100, 'A');
        const std::string qual(100, 'I');
        const uint32_t cigar[] = {
            bam_cigar_gen(70, BAM_CMATCH),
            bam_cigar_gen(30, BAM_CSOFT_CLIP),
        };
        const int set_ret = bam_set1(
            record,
            qname.size() + 1,
            qname.c_str(),
            0,
            0,
            1000,
            60,
            sizeof(cigar) / sizeof(cigar[0]),
            cigar,
            -1,
            -1,
            0,
            seq.size(),
            seq.c_str(),
            qual.c_str(),
            0);
        assert(set_ret >= 0);

        const char sa_value[] = "1,2000,+,30S70M,60,1;";
        const int aux_ret = bam_aux_append(
            record,
            "SA",
            'Z',
            static_cast<int>(sizeof(sa_value)),
            reinterpret_cast<const uint8_t*>(sa_value));
        assert(aux_ret == 0);

        const std::vector<const bam1_t*> records = {record};
        LinearBinComponentModule module;
        const auto components = module.build(records, "1", 0, 900, 1200);
        assert(components.size() == 1);

        int split_hint_count = 0;
        for (const auto& bp : components.front().breakpoint_candidates) {
            if ((bp.class_mask & kCandidateSplitSaSupplementary) != 0 && bp.clip_len == 0) {
                split_hint_count += 1;
                assert(bp.pos == 1070);
            }
        }
        assert(split_hint_count == 1);

        bam_destroy1(record);
    }

    {
        bam1_t* record = bam_init1();
        assert(record != nullptr);

        const std::string qname = "read_double_softclip";
        const std::string seq(100, 'A');
        const std::string qual(100, 'I');
        const uint32_t cigar[] = {
            bam_cigar_gen(30, BAM_CSOFT_CLIP),
            bam_cigar_gen(40, BAM_CMATCH),
            bam_cigar_gen(30, BAM_CSOFT_CLIP),
        };
        const int set_ret = bam_set1(
            record,
            qname.size() + 1,
            qname.c_str(),
            0,
            0,
            1000,
            60,
            sizeof(cigar) / sizeof(cigar[0]),
            cigar,
            -1,
            -1,
            0,
            seq.size(),
            seq.c_str(),
            qual.c_str(),
            0);
        assert(set_ret >= 0);

        const std::vector<const bam1_t*> records = {record};
        LinearBinComponentModule module;
        const auto components = module.build(records, "1", 0, 900, 1200);
        assert(components.size() == 1);
        assert(components.front().evidence_soft_clip_count == 1);

        int softclip_candidates = 0;
        for (const auto& bp : components.front().breakpoint_candidates) {
            if (bp.clip_len > 0) {
                softclip_candidates += 1;
            }
        }
        assert(softclip_candidates == 1);

        bam_destroy1(record);
    }

    {
        bam1_t* record = bam_init1();
        assert(record != nullptr);

        const std::string qname = "read_large_deletion_only";
        const std::string seq(100, 'A');
        const std::string qual(100, 'I');
        const uint32_t cigar[] = {
            bam_cigar_gen(50, BAM_CMATCH),
            bam_cigar_gen(60, BAM_CDEL),
            bam_cigar_gen(50, BAM_CMATCH),
        };
        const int set_ret = bam_set1(
            record,
            qname.size() + 1,
            qname.c_str(),
            0,
            0,
            1000,
            60,
            sizeof(cigar) / sizeof(cigar[0]),
            cigar,
            -1,
            -1,
            0,
            seq.size(),
            seq.c_str(),
            qual.c_str(),
            0);
        assert(set_ret >= 0);

        const std::vector<const bam1_t*> records = {record};
        LinearBinComponentModule module;
        const auto components = module.build(records, "1", 0, 900, 1200);
        assert(components.empty());

        bam_destroy1(record);
    }

    {
        EventGenotypeInput in;
        in.alt_struct_reads = 6;
        in.ref_span_reads = 4;
        in.min_depth = 3;
        in.min_gq = 20;
        in.error_rate = 0.02;

        const auto d = genotype_event_from_alt_vs_ref(in);
        assert(d.best_gt == "0/1");
        assert(d.allele_fraction > 0.55 && d.allele_fraction < 0.65);
        assert(d.gq >= 20);
        assert(d.pass);
    }

    {
        EventGenotypeInput in;
        in.alt_struct_reads = 5;
        in.ref_span_reads = 1;
        in.min_depth = 3;
        in.min_gq = 20;
        in.error_rate = 0.02;

        const auto d = genotype_event_from_alt_vs_ref(in);
        assert(d.best_gt == "1/1");
        assert(d.allele_fraction > 0.80);
        assert(d.gq >= 20);
        assert(d.pass);
    }

    {
        EventGenotypeInput in;
        in.alt_struct_reads = 2;
        in.ref_span_reads = 0;
        in.min_depth = 3;
        in.min_gq = 20;
        in.error_rate = 0.02;

        const auto d = genotype_event_from_alt_vs_ref(in);
        assert(d.best_gt == "1/1");
        assert(d.allele_fraction > 0.95);
        assert(d.gq >= 20);
        assert(d.pass);
    }

    {
        EventGenotypeInput in;
        in.alt_struct_reads = 0;
        in.ref_span_reads = 8;
        in.min_depth = 3;
        in.min_gq = 20;
        in.error_rate = 0.02;

        const auto d = genotype_event_from_alt_vs_ref(in);
        assert(d.best_gt == "0/0");
        assert(!d.pass);
    }

    {
        FinalBoundaryInput in;
        in.left_ref_start = 100;
        in.left_ref_end = 220;
        in.right_ref_start = 211;
        in.right_ref_end = 330;
        in.tsd_min_len = 3;
        in.tsd_max_len = 50;

        const auto d = check_boundary_consistency(in);
        assert(d.pass);
        assert(d.boundary_type == "TSD");
        assert(d.boundary_len == 9);
        assert(d.qc == "PASS_BOUNDARY_TSD");
    }

    {
        FinalBoundaryInput in;
        in.left_ref_start = 100;
        in.left_ref_end = 220;
        in.right_ref_start = 220;
        in.right_ref_end = 340;

        const auto d = check_boundary_consistency(in);
        assert(d.pass);
        assert(d.boundary_type == "BLUNT");
        assert(d.qc == "PASS_BOUNDARY_BLUNT");
    }

    {
        FinalBoundaryInput in;
        in.left_ref_start = 100;
        in.left_ref_end = 220;
        in.right_ref_start = 227;
        in.right_ref_end = 340;

        const auto d = check_boundary_consistency(in);
        assert(d.pass);
        assert(d.boundary_type == "SMALL_DEL");
        assert(d.boundary_len == 7);
        assert(d.qc == "PASS_BOUNDARY_SMALL_DEL");
    }

    {
        FinalBoundaryInput in;
        in.left_ref_start = 100;
        in.left_ref_end = 220;
        in.right_ref_start = 160;
        in.right_ref_end = 340;
        in.tsd_min_len = 3;
        in.tsd_max_len = 50;

        const auto d = check_boundary_consistency(in);
        assert(!d.pass);
        assert(d.qc == "REJECT_BOUNDARY_TSD_RANGE");
    }

    {
        FinalTeAcceptanceInput in;
        in.event_existence_pass = true;
        in.event_closure_pass = true;
        in.te_sequence_pass = true;
        in.boundary_pass = true;

        const auto d = evaluate_final_te_acceptance(in);
        assert(d.pass);
        assert(d.qc == "PASS_FINAL_TE_CALL");
    }

    {
        FinalTeAcceptanceInput in;
        in.event_existence_pass = true;
        in.event_closure_pass = false;
        in.te_sequence_pass = true;
        in.boundary_pass = true;

        const auto d = evaluate_final_te_acceptance(in);
        assert(!d.pass);
        assert(d.qc == "REJECT_EVENT_CLOSURE");
    }

    {
        EventGenotypeInput in;
        in.alt_struct_reads = 6;
        in.ref_span_reads = 4;
        in.min_depth = 3;
        in.min_gq = 20;
        in.error_rate = 0.02;

        const auto e = build_event_existence_evidence(in);
        assert(e.best_gt == "0/1");
        assert(e.depth == 10);
        assert(e.gq >= 20);
        assert(e.score > 0.0);
    }

    {
        EventSegmentationEvidence seg;
        seg.has_consensus = true;
        seg.has_left_flank = true;
        seg.has_right_flank = true;
        seg.has_insert_seq = true;
        seg.pair_valid = true;
        seg.left_align_len = 95;
        seg.right_align_len = 90;
        seg.left_identity = 0.98;
        seg.right_identity = 0.97;
        seg.insert_len = 210;
        seg.score = 1.9;
        seg.qc = "PASS_EVENT_SEGMENTATION";

        TEAlignmentEvidence te;
        te.best_family = "L1";
        te.best_subfamily = "L1-3";
        te.best_identity = 0.95;
        te.best_query_coverage = 0.94;
        te.cross_family_margin = 0.30;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT";

        BoundaryEvidence boundary;
        boundary.geometry_defined = true;
        boundary.canonical_pass = false;
        boundary.evidence_consistent = true;
        boundary.boundary_type = "NONCANONICAL";
        boundary.boundary_len = 64;
        boundary.score = 0.25;
        boundary.qc = "PASS_BOUNDARY_NONCANONICAL_CONSISTENT";

        EventGenotypeInput in;
        in.alt_struct_reads = 4;
        in.ref_span_reads = 5;
        in.min_depth = 3;
        in.min_gq = 20;
        in.error_rate = 0.02;

        const auto existence = build_event_existence_evidence(in);
        const auto joint = evaluate_joint_hypotheses(existence, seg, te, boundary);
        assert(joint.emit_te_call);
        assert(joint.best.kind == FinalHypothesisKind::kTeResolved);
    }

    {
        EventSegmentationEvidence seg;
        seg.has_consensus = true;
        seg.has_left_flank = true;
        seg.has_right_flank = true;
        seg.has_insert_seq = true;
        seg.pair_valid = true;
        seg.left_align_len = 98;
        seg.right_align_len = 58;
        seg.left_identity = 0.96;
        seg.right_identity = 0.93;
        seg.insert_len = 472;
        seg.score = 1.3;
        seg.qc = "PASS_EVENT_SEGMENTATION";

        TEAlignmentEvidence te;
        te.best_family = "Unknown";
        te.best_subfamily = "Unknown-62";
        te.best_identity = 0.853526;
        te.best_query_coverage = 1.0;
        te.cross_family_margin = 0.324078;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT";

        BoundaryEvidence boundary;
        boundary.geometry_defined = true;
        boundary.canonical_pass = true;
        boundary.evidence_consistent = true;
        boundary.boundary_type = "TSD";
        boundary.boundary_len = 18;
        boundary.score = 1.0;
        boundary.qc = "PASS_BOUNDARY_TSD";

        EventGenotypeInput in;
        in.alt_struct_reads = 7;
        in.ref_span_reads = 1;
        in.min_depth = 3;
        in.min_gq = 20;
        in.error_rate = 0.02;

        const auto existence = build_event_existence_evidence(in);
        const auto joint = evaluate_joint_hypotheses(existence, seg, te, boundary);
        assert(joint.best.kind == FinalHypothesisKind::kTeResolved);
        assert(joint.runner_up.kind == FinalHypothesisKind::kTeUnknown);
        assert(joint.emit_te_call);
        assert(!joint.emit_unknown_te);
        assert(joint.final_qc == "PASS_FINAL_TE_CALL");
    }

    {
        EventSegmentationEvidence seg;
        seg.has_consensus = true;
        seg.has_left_flank = true;
        seg.has_right_flank = true;
        seg.has_insert_seq = false;
        seg.pair_valid = false;
        seg.score = -2.0;
        seg.qc = "NO_EVENT_SEGMENTATION";

        TEAlignmentEvidence te;

        BoundaryEvidence boundary;
        boundary.geometry_defined = false;
        boundary.score = -2.0;
        boundary.qc = "REJECT_BOUNDARY_UNSET";

        EventGenotypeInput in;
        in.alt_struct_reads = 8;
        in.ref_span_reads = 0;
        in.min_depth = 3;
        in.min_gq = 20;
        in.error_rate = 0.02;

        const auto existence = build_event_existence_evidence(in);
        const auto joint = evaluate_joint_hypotheses(existence, seg, te, boundary);
        assert(!joint.emit_te_call);
        assert(joint.best.kind != FinalHypothesisKind::kTeResolved);
        assert(joint.best.kind != FinalHypothesisKind::kTeUnknown);
    }

    {
        EventSegmentationEvidence seg;
        seg.has_consensus = true;
        seg.has_left_flank = true;
        seg.has_right_flank = true;
        seg.has_insert_seq = true;
        seg.pair_valid = true;
        seg.left_align_len = 88;
        seg.right_align_len = 91;
        seg.left_identity = 0.95;
        seg.right_identity = 0.96;
        seg.insert_len = 180;
        seg.score = 1.5;
        seg.qc = "PASS_EVENT_SEGMENTATION";

        TEAlignmentEvidence te;
        te.best_family = "UNKNOWN";
        te.best_subfamily = "UNKNOWN";
        te.best_identity = 0.79;
        te.best_query_coverage = 0.90;
        te.cross_family_margin = 0.01;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";

        BoundaryEvidence boundary;
        boundary.geometry_defined = true;
        boundary.canonical_pass = true;
        boundary.evidence_consistent = true;
        boundary.boundary_type = "TSD";
        boundary.boundary_len = 11;
        boundary.score = 1.0;
        boundary.qc = "PASS_BOUNDARY_TSD";

        EventGenotypeInput in;
        in.alt_struct_reads = 5;
        in.ref_span_reads = 2;
        in.min_depth = 3;
        in.min_gq = 20;
        in.error_rate = 0.02;

        const auto existence = build_event_existence_evidence(in);
        const auto joint = evaluate_joint_hypotheses(existence, seg, te, boundary);
        assert(joint.emit_te_call);
        assert(joint.best.kind == FinalHypothesisKind::kTeUnknown);
    }

    {
        EventSegmentationEvidence seg;
        seg.has_consensus = true;
        seg.has_left_flank = true;
        seg.has_right_flank = true;
        seg.has_insert_seq = true;
        seg.pair_valid = true;
        seg.left_align_len = 110;
        seg.right_align_len = 110;
        seg.left_identity = 0.95;
        seg.right_identity = 0.95;
        seg.insert_len = 248;
        seg.score = 1.0;
        seg.qc = "PASS_EVENT_SEGMENTATION";

        TEAlignmentEvidence te;
        te.best_family = "UNKNOWN";
        te.best_subfamily = "UNKNOWN";
        te.best_identity = 0.65942;
        te.best_query_coverage = 0.979839;
        te.cross_family_margin = 0.00597882;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";

        BoundaryEvidence boundary;
        boundary.geometry_defined = true;
        boundary.canonical_pass = false;
        boundary.evidence_consistent = true;
        boundary.boundary_type = "NONCANONICAL";
        boundary.boundary_len = 55;
        boundary.score = 0.25;
        boundary.qc = "PASS_BOUNDARY_NONCANONICAL_CONSISTENT";

        EventGenotypeInput in;
        in.alt_struct_reads = 7;
        in.ref_span_reads = 4;
        in.min_depth = 3;
        in.min_gq = 20;
        in.error_rate = 0.02;

        const auto existence = build_event_existence_evidence(in);
        const auto joint = evaluate_joint_hypotheses(existence, seg, te, boundary);
        assert(joint.emit_te_call);
        assert(joint.emit_unknown_te);
        assert(joint.best.kind == FinalHypothesisKind::kTeUnknown);
    }

    {
        EventSegmentationEvidence seg;
        seg.has_consensus = true;
        seg.has_left_flank = true;
        seg.has_right_flank = true;
        seg.has_insert_seq = true;
        seg.pair_valid = true;
        seg.left_align_len = 64;
        seg.right_align_len = 91;
        seg.left_identity = 0.92;
        seg.right_identity = 0.96;
        seg.insert_len = 1624;
        seg.score = 0.291;
        seg.qc = "PASS_EVENT_SEGMENTATION";

        TEAlignmentEvidence te;
        te.best_family = "UNKNOWN";
        te.best_subfamily = "UNKNOWN";
        te.best_identity = 0.858844;
        te.best_query_coverage = 0.700123;
        te.cross_family_margin = 0.0478658;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";

        BoundaryEvidence boundary;
        boundary.geometry_defined = true;
        boundary.canonical_pass = true;
        boundary.evidence_consistent = true;
        boundary.boundary_type = "TSD";
        boundary.boundary_len = 17;
        boundary.score = 1.0;
        boundary.qc = "PASS_BOUNDARY_TSD";

        EventGenotypeInput in;
        in.alt_struct_reads = 6;
        in.ref_span_reads = 9;
        in.min_depth = 3;
        in.min_gq = 20;
        in.error_rate = 0.02;

        const auto existence = build_event_existence_evidence(in);
        const auto joint = evaluate_joint_hypotheses(existence, seg, te, boundary);
        assert(joint.emit_te_call);
        assert(joint.emit_unknown_te);
        assert(joint.best.kind == FinalHypothesisKind::kTeUnknown);
    }

    {
        JointHypothesisScore best_te;
        best_te.kind = FinalHypothesisKind::kTeResolved;
        best_te.total = 1.60;

        JointHypothesisScore best_non_te;
        best_non_te.kind = FinalHypothesisKind::kInsertionNonTe;
        best_non_te.total = 0.70;

        assert(should_emit_te_call(best_te, best_non_te));
    }

    {
        JointHypothesisScore best_te;
        best_te.kind = FinalHypothesisKind::kTeUnknown;
        best_te.total = 1.49;

        JointHypothesisScore best_non_te;
        best_non_te.kind = FinalHypothesisKind::kInsertionNonTe;
        best_non_te.total = 0.20;

        assert(!should_emit_te_call(best_te, best_non_te));
    }

    {
        JointHypothesisScore best_te;
        best_te.kind = FinalHypothesisKind::kTeResolved;
        best_te.total = 1.80;

        JointHypothesisScore best_non_te;
        best_non_te.kind = FinalHypothesisKind::kInsertionNonTe;
        best_non_te.total = 1.10;

        assert(!should_emit_te_call(best_te, best_non_te));
    }

    {
        EventSegmentation seg;
        seg.pass = true;
        seg.qc_reason = "PASS_EVENT_SEGMENTATION";
        seg.left_flank_align_len = 94;
        seg.right_flank_align_len = 97;
        seg.left_flank_identity = 0.95;
        seg.right_flank_identity = 0.96;
        seg.insert_seq = std::string(233, 'A');

        const auto evidence = analyze_event_segmentation_for_test(true, seg);
        assert(evidence.has_insert_seq);
        assert(evidence.pair_valid);
        assert(evidence.score > 0.0);
    }

    {
        FinalBoundaryInput in;
        in.left_ref_start = 100;
        in.left_ref_end = 220;
        in.right_ref_start = 156;
        in.right_ref_end = 340;

        const auto evidence = evaluate_boundary_evidence(in, 90);
        assert(evidence.geometry_defined);
        assert(!evidence.canonical_pass);
        assert(evidence.evidence_consistent);
        assert(evidence.score > 0.0);
    }

    return 0;
}
