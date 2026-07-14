#ifdef NDEBUG
#undef NDEBUG
#endif
#include "decision_policy.h"
#include "pipeline.h"

#include <cassert>
#include <string>
#include <vector>

namespace {

void set_bam_breakdown(
    placer::EventGenotypeInput& input,
    int32_t split_reads,
    int32_t indel_reads,
    int32_t left_clip_reads,
    int32_t right_clip_reads) {
    input.alt_split_reads = split_reads;
    input.alt_indel_reads = indel_reads;
    input.alt_left_clip_reads = left_clip_reads;
    input.alt_right_clip_reads = right_clip_reads;
}

void set_bam_breakdown(
    placer::EventExistenceEvidence& evidence,
    int32_t split_reads,
    int32_t indel_reads,
    int32_t left_clip_reads,
    int32_t right_clip_reads) {
    evidence.alt_split_reads = split_reads;
    evidence.alt_indel_reads = indel_reads;
    evidence.alt_left_clip_reads = left_clip_reads;
    evidence.alt_right_clip_reads = right_clip_reads;
}

placer::EventExistenceEvidence make_igv_existence(
    int32_t alt_struct_reads,
    int32_t alt_split_reads,
    int32_t alt_indel_reads,
    int32_t alt_left_clip_reads,
    int32_t alt_right_clip_reads,
    int32_t ref_span_reads,
    int32_t gq,
    double score) {
    placer::EventExistenceEvidence existence;
    existence.best_gt = "0/1";
    existence.alt_struct_reads = alt_struct_reads;
    set_bam_breakdown(
        existence,
        alt_split_reads,
        alt_indel_reads,
        alt_left_clip_reads,
        alt_right_clip_reads);
    existence.ref_span_reads = ref_span_reads;
    existence.depth = alt_struct_reads + ref_span_reads;
    existence.af = existence.depth > 0
        ? static_cast<double>(alt_struct_reads) / static_cast<double>(existence.depth)
        : 0.0;
    existence.gq = gq;
    existence.score = score;
    return existence;
}

placer::EventSegmentationEvidence make_igv_complete_segmentation(
    int32_t insert_len,
    double score,
    int32_t left_align_len = 80,
    int32_t right_align_len = 80) {
    placer::EventSegmentationEvidence segmentation;
    segmentation.has_consensus = true;
    segmentation.has_left_flank = true;
    segmentation.has_right_flank = true;
    segmentation.has_insert_seq = true;
    segmentation.pair_valid = true;
    segmentation.left_align_len = left_align_len;
    segmentation.right_align_len = right_align_len;
    segmentation.left_identity = 0.95;
    segmentation.right_identity = 0.95;
    segmentation.insert_len = insert_len;
    segmentation.score = score;
    segmentation.qc = "PASS_EVENT_SEGMENTATION";
    return segmentation;
}

placer::BoundaryEvidence make_igv_boundary(
    const std::string& boundary_type,
    int32_t boundary_len) {
    placer::BoundaryEvidence boundary;
    boundary.geometry_defined = true;
    boundary.canonical_pass = true;
    boundary.evidence_consistent = true;
    boundary.boundary_type = boundary_type;
    boundary.boundary_len = boundary_len;
    boundary.score = 1.0;
    boundary.qc = boundary_type == "TSD"
        ? "PASS_BOUNDARY_TSD"
        : (boundary_type == "BLUNT"
               ? "PASS_BOUNDARY_BLUNT"
               : "PASS_BOUNDARY_SMALL_DEL");
    return boundary;
}

placer::TEAlignmentEvidence make_igv_te(
    const std::string& family,
    const std::string& subfamily,
    double identity,
    double query_coverage,
    double cross_family_margin,
    bool pass,
    const std::string& qc_reason,
    const std::string& sequence_model_label,
    double sequence_model_score) {
    placer::TEAlignmentEvidence te;
    te.best_family = family;
    te.best_subfamily = subfamily;
    te.best_identity = identity;
    te.best_query_coverage = query_coverage;
    te.cross_family_margin = cross_family_margin;
    te.pass = pass;
    te.qc_reason = qc_reason;
    te.sequence_model_label = sequence_model_label;
    te.sequence_model_score = sequence_model_score;
    return te;
}

void assert_igv_recall_emit(
    const placer::EventExistenceEvidence& existence,
    const placer::EventSegmentationEvidence& segmentation,
    const placer::TEAlignmentEvidence& te,
    const placer::BoundaryEvidence& boundary,
    bool expect_unknown) {
    const auto joint = placer::evaluate_joint_hypotheses(
        existence,
        segmentation,
        te,
        boundary);
    assert(joint.final_qc.find("RESCUE") == std::string::npos);
    if (joint.emit_te_call) {
        assert(joint.emit_unknown_te == expect_unknown);
        assert(joint.final_qc == "PASS_TE_CLOSED" ||
               joint.final_qc == "PASS_TE_IMPRECISE");
    } else {
        assert(joint.final_qc == "TE_AMBIGUOUS" ||
               joint.final_qc == "TE_LFDR_HIGH" ||
               joint.final_qc == "PASS_NONTE_INSERTION" ||
               joint.final_qc == "PASS_STRUCTURAL_INSERTION" ||
               joint.final_qc == "REFERENCE_OR_ARTIFACT" ||
               joint.final_qc == "NO_CALL_INCOMPLETE");
    }
}

}  // namespace

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
        result.final_calls.push_back(make_call(1030, 5, 12, 0.30));
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

        int softclip_candidates = 0;
        for (const auto& bp : components.front().breakpoint_candidates) {
            if (bp.clip_len > 0) {
                softclip_candidates += 1;
            }
            if ((bp.class_mask & kCandidateSplitSaSupplementary) != 0) {
                assert(false && "SA without query-gap insertion should not be split support");
            }
        }
        assert(components.front().split_sa_read_indices.empty());
        assert(components.front().soft_clip_read_indices.size() == 1);
        assert(softclip_candidates == 1);

        bam_destroy1(record);
    }

    {
        bam1_t* record = bam_init1();
        assert(record != nullptr);

        const std::string qname = "read_trailing_softclip_with_sa";
        const std::string seq(200, 'A');
        const std::string qual(200, 'I');
        const uint32_t cigar[] = {
            bam_cigar_gen(80, BAM_CMATCH),
            bam_cigar_gen(120, BAM_CSOFT_CLIP),
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

        const char sa_value[] = "1,1081,+,140S60M,60,1;";
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
                assert(bp.pos == 1080);
            }
        }
        assert(split_hint_count == 1);
        assert(components.front().split_sa_read_indices.size() == 1);
        assert(components.front().soft_clip_read_indices.empty());

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
        assert(components.front().soft_clip_read_indices.size() == 1);

        int softclip_candidates = 0;
        for (const auto& bp : components.front().breakpoint_candidates) {
            if (bp.clip_len > 0) {
                softclip_candidates += 1;
            }
        }
        assert(softclip_candidates == 2);

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
        // 5 alt / 1 ref (AF 0.83) is a genuinely ambiguous het-vs-hom call at low
        // depth: the beta-binomial genotyper reports it as a low-confidence het
        // that does not pass the GQ floor.
        assert(d.best_gt == "0/1");
        assert(d.allele_fraction > 0.80);
        assert(d.gq < 20);
        assert(!d.pass);
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
        assert(d.gq < 20);
        assert(!d.pass);
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
        // Sample-level overdispersion estimator (L3): tight ~binomial het sites
        // give rho ~ 0; widely dispersed allele balances give a large rho; too
        // few informative sites fall back to the default.
        std::vector<AltDepthObservation> tight;
        std::vector<AltDepthObservation> dispersed;
        for (int i = 0; i < 40; ++i) {
            const int alt_tight = 10 + ((i % 3) - 1);  // 9,10,11 around 0.5
            tight.push_back({alt_tight, 20});
            const int alt_dispersed = (i % 2 == 0) ? 4 : 16;  // 0.2 / 0.8
            dispersed.push_back({alt_dispersed, 20});
        }
        const double rho_tight = estimate_alt_depth_overdispersion(tight, 0.02);
        const double rho_dispersed = estimate_alt_depth_overdispersion(dispersed, 0.02);
        assert(rho_tight >= 0.0 && rho_tight < 0.05);
        assert(rho_dispersed > 0.2);
        assert(rho_dispersed > rho_tight);

        std::vector<AltDepthObservation> too_few = {{10, 20}, {9, 20}, {11, 20}};
        assert(estimate_alt_depth_overdispersion(too_few, 0.037) == 0.037);
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
        te.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        te.sequence_model_score = 0.40;

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
        set_bam_breakdown(in, 3, 1, 0, 0);
        in.ref_span_reads = 5;
        in.min_depth = 3;
        in.min_gq = 20;
        in.error_rate = 0.02;

        const auto existence = build_event_existence_evidence(in);
        const auto joint = evaluate_joint_hypotheses(existence, seg, te, boundary);
        assert(!joint.emit_te_call);
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
        te.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        te.sequence_model_score = 0.40;

        BoundaryEvidence boundary;
        boundary.geometry_defined = true;
        boundary.canonical_pass = true;
        boundary.evidence_consistent = true;
        boundary.boundary_type = "TSD";
        boundary.boundary_len = 18;
        boundary.score = 1.0;
        boundary.qc = "PASS_BOUNDARY_TSD";

        EventGenotypeInput in;
        in.alt_struct_reads = 12;
        set_bam_breakdown(in, 4, 4, 2, 2);
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
        assert(joint.final_qc == "PASS_TE_CLOSED");
        assert(joint.te_posterior >= 0.85);
        assert(joint.te_vs_artifact_log_odds > 1.0);
        assert(joint.posterior_qc == "PASS_TE_POSTERIOR");
        assert(joint.lfdr_qc == "PASS_TE_LFDR");
        assert(joint.worst_case_lfdr <= 0.10);
        assert(joint.mechanistic_lower_log_bf_te_vs_artifact > 0.0);
        assert(joint.mechanistic_lower_log_bf_te_vs_non_te > 0.0);
        assert(joint.mechanistic_blocks != "NA");
        assert(joint.robust_mechanistic_worst_case_lfdr <= 1.0);
        assert(joint.robust_mechanistic_qc == "PASS_TE_LFDR");
    }

    {
        EventExistenceEvidence existence;
        existence.best_gt = "0/1";
        existence.af = 0.82;
        existence.gq = 70;
        existence.alt_struct_reads = 14;
        set_bam_breakdown(existence, 0, 0, 0, 14);
        existence.ref_span_reads = 4;
        existence.depth = 18;
        existence.score = 2.5;

        EventSegmentationEvidence seg;
        seg.has_consensus = true;
        seg.has_left_flank = true;
        seg.has_right_flank = true;
        seg.has_insert_seq = true;
        seg.pair_valid = true;
        seg.left_align_len = 90;
        seg.right_align_len = 88;
        seg.left_identity = 0.95;
        seg.right_identity = 0.95;
        seg.insert_len = 431;
        seg.score = 1.2;
        seg.qc = "PASS_EVENT_SEGMENTATION";

        TEAlignmentEvidence te;
        te.best_family = "hAT-Ac";
        te.best_subfamily = "hAT-Ac-10";
        te.best_identity = 0.948;
        te.best_query_coverage = 0.97;
        te.cross_family_margin = 0.35;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT";
        te.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        te.sequence_model_score = 0.40;

        BoundaryEvidence boundary;
        boundary.geometry_defined = true;
        boundary.canonical_pass = true;
        boundary.evidence_consistent = true;
        boundary.boundary_type = "TSD";
        boundary.boundary_len = 9;
        boundary.score = 1.0;
        boundary.qc = "PASS_BOUNDARY_TSD";

        const auto lfdr = evaluate_latent_mechanism_lfdr(
            existence,
            seg,
            te,
            boundary);
        assert(lfdr.worst_case_lfdr > 0.10);
        assert(lfdr.lfdr_qc == "TE_LFDR_HIGH");

        const auto joint = evaluate_joint_hypotheses(existence, seg, te, boundary);
        assert(!joint.emit_te_call);
        assert(joint.lfdr_qc == "TE_LFDR_HIGH");
        assert(joint.mechanistic_ref_conflict_signal > 0.0);
        assert(joint.mechanistic_lower_log_bf_te_vs_artifact < 0.0 ||
               joint.lfdr_qc == "TE_LFDR_HIGH");
        assert(joint.robust_mechanistic_qc == "PASS_TE_LFDR");
        assert(joint.robust_mechanistic_worst_case_lfdr <= 0.10);
    }

    {
        EventExistenceEvidence existence;
        existence.best_gt = "0/1";
        existence.af = 0.975;
        existence.gq = 99;
        existence.alt_struct_reads = 40;
        set_bam_breakdown(existence, 9, 8, 12, 11);
        existence.ref_span_reads = 1;
        existence.depth = 41;
        existence.score = 4.0;

        EventSegmentationEvidence seg;
        seg.has_consensus = true;
        seg.has_left_flank = true;
        seg.has_right_flank = true;
        seg.has_insert_seq = true;
        seg.pair_valid = true;
        seg.left_align_len = 91;
        seg.right_align_len = 88;
        seg.left_identity = 0.95;
        seg.right_identity = 0.94;
        seg.insert_len = 431;
        seg.score = 1.25;
        seg.qc = "PASS_EVENT_SEGMENTATION";

        TEAlignmentEvidence te;
        te.best_family = "UNKNOWN";
        te.best_subfamily = "UNKNOWN";
        te.best_identity = 0.655;
        te.best_query_coverage = 0.323;
        te.cross_family_margin = 0.02;
        te.pass = false;
        te.qc_reason = "TE_ALIGNMENT_LOW_COVERAGE";
        te.sequence_model_label = "TE_MODEL_OUTLIER";
        te.sequence_model_score = 0.05;
        te.annotation_confidence = "LOW";

        BoundaryEvidence boundary;
        boundary.geometry_defined = true;
        boundary.canonical_pass = true;
        boundary.evidence_consistent = true;
        boundary.boundary_type = "SMALL_DEL";
        boundary.boundary_len = 10;
        boundary.score = 0.85;
        boundary.qc = "PASS_BOUNDARY_SMALL_DEL";

        const auto joint = evaluate_joint_hypotheses(existence, seg, te, boundary);
        assert(!joint.emit_te_call);
        assert(joint.emit_structural_event_call);
        assert(joint.emit_unknown_te);
        assert(joint.final_qc == "PASS_STRUCTURAL_INSERTION");
        assert(joint.best.kind == FinalHypothesisKind::kInsertionNonTe);
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
        te.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        te.sequence_model_score = 0.40;

        BoundaryEvidence boundary;
        boundary.geometry_defined = true;
        boundary.canonical_pass = true;
        boundary.evidence_consistent = true;
        boundary.boundary_type = "TSD";
        boundary.boundary_len = 11;
        boundary.score = 1.0;
        boundary.qc = "PASS_BOUNDARY_TSD";

        EventGenotypeInput in;
        in.alt_struct_reads = 12;
        set_bam_breakdown(in, 4, 4, 2, 2);
        in.ref_span_reads = 2;
        in.min_depth = 3;
        in.min_gq = 20;
        in.error_rate = 0.02;

        const auto existence = build_event_existence_evidence(in);
        const auto joint = evaluate_joint_hypotheses(existence, seg, te, boundary);
        assert(joint.emit_te_call);
        assert(joint.emit_unknown_te);
        assert(joint.lfdr_qc == "TE_LFDR_HIGH");
    }

    {
        EventSegmentationEvidence seg;
        seg.has_consensus = true;
        seg.has_left_flank = true;
        seg.has_right_flank = true;
        seg.has_insert_seq = true;
        seg.pair_valid = true;
        seg.left_align_len = 82;
        seg.right_align_len = 78;
        seg.left_identity = 0.95;
        seg.right_identity = 0.94;
        seg.insert_len = 600;
        seg.score = 1.15;
        seg.qc = "PASS_EVENT_SEGMENTATION";

        TEAlignmentEvidence te;
        te.best_family = "UNKNOWN";
        te.best_subfamily = "UNKNOWN";
        te.best_identity = 0.555;
        te.best_query_coverage = 0.96;
        te.cross_family_margin = 0.10;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";
        te.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        te.sequence_model_score = 0.40;
        te.annotation_confidence = "LOW";
        te.annotation_residual_fraction = 0.04;
        te.annotation_masked_fraction = 0.0;

        BoundaryEvidence boundary;
        boundary.geometry_defined = true;
        boundary.canonical_pass = true;
        boundary.evidence_consistent = true;
        boundary.boundary_type = "SMALL_DEL";
        boundary.boundary_len = 12;
        boundary.score = 1.0;
        boundary.qc = "PASS_BOUNDARY_SMALL_DEL";

        EventGenotypeInput in;
        in.alt_struct_reads = 16;
        set_bam_breakdown(in, 2, 2, 6, 6);
        in.ref_span_reads = 0;
        in.min_depth = 3;
        in.min_gq = 20;
        in.error_rate = 0.02;

        const auto existence = build_event_existence_evidence(in);
        const auto joint = evaluate_joint_hypotheses(existence, seg, te, boundary);
        assert(joint.emit_te_call);
        assert(joint.emit_unknown_te);
        assert(joint.final_qc == "PASS_TE_CLOSED");
        assert(joint.posterior_qc == "PASS_TE_POSTERIOR");
        assert(joint.lfdr_qc == "TE_LFDR_HIGH");
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
        te.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        te.sequence_model_score = 0.40;

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
        set_bam_breakdown(in, 3, 3, 1, 1);
        in.ref_span_reads = 4;
        in.min_depth = 3;
        in.min_gq = 20;
        in.error_rate = 0.02;

        const auto existence = build_event_existence_evidence(in);
        const auto joint = evaluate_joint_hypotheses(existence, seg, te, boundary);
        assert(joint.emit_te_call);
        assert(joint.emit_unknown_te);
        assert(joint.final_qc == "PASS_TE_CLOSED");
        assert(joint.lfdr_qc == "TE_LFDR_HIGH");
        assert(joint.final_qc.find("RESCUE") == std::string::npos);
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
        te.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        te.sequence_model_score = 0.40;

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
        set_bam_breakdown(in, 3, 3, 0, 0);
        in.ref_span_reads = 9;
        in.min_depth = 3;
        in.min_gq = 20;
        in.error_rate = 0.02;

        const auto existence = build_event_existence_evidence(in);
        const auto joint = evaluate_joint_hypotheses(existence, seg, te, boundary);
        // A confident TE call abstains on ambiguous UNKNOWN identity, but this is
        // a decisive closed structural event (1624 bp insert, TSD boundary, 6
        // precise split/indel reads): it is emitted as a real non-reference
        // structural insertion instead of being dropped as ambiguous.
        assert(!joint.emit_te_call);
        assert(joint.best_explanation == "TE");
        assert(joint.emit_structural_event_call);
        assert(joint.final_qc == "PASS_STRUCTURAL_INSERTION");
        assert(joint.final_qc.find("RESCUE") == std::string::npos);
    }

    {
        EventExistenceEvidence existence;
        existence.best_gt = "0/1";
        existence.af = 0.5;
        existence.gq = 75;
        existence.alt_struct_reads = 3;
        set_bam_breakdown(existence, 3, 0, 0, 0);
        existence.ref_span_reads = 3;
        existence.depth = 6;
        existence.score = 2.75;

        EventSegmentationEvidence seg;
        seg.has_consensus = true;
        seg.has_left_flank = true;
        seg.has_right_flank = false;
        seg.has_insert_seq = true;
        seg.pair_valid = true;
        seg.left_align_len = 64;
        seg.right_align_len = 0;
        seg.left_identity = 0.95;
        seg.right_identity = 0.0;
        seg.insert_len = 613;
        seg.score = 0.6825;
        seg.qc = "PASS_EVENT_SEGMENTATION_ONE_SIDED_LEFT";

        TEAlignmentEvidence te;
        te.best_family = "UNKNOWN";
        te.best_subfamily = "UNKNOWN";
        te.best_identity = 0.52957;
        te.best_query_coverage = 1.0;
        te.cross_family_margin = 0.506731;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";
        te.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        te.sequence_model_score = 0.40;

        BoundaryEvidence boundary;
        boundary.geometry_defined = false;
        boundary.canonical_pass = false;
        boundary.evidence_consistent = false;
        boundary.boundary_type = "REJECT";
        boundary.boundary_len = 0;
        boundary.score = -2.0;
        boundary.qc = "REJECT_BOUNDARY_INVALID_REF_SEGMENTS";

        const auto joint = evaluate_joint_hypotheses(existence, seg, te, boundary);
        assert(!joint.emit_te_call);
        assert(!joint.emit_unknown_te);
        assert(joint.emit_evidence_te_call);
        assert(joint.final_qc == "TE_AMBIGUOUS");
        assert(joint.te_posterior < 0.65);
        assert(joint.posterior_qc == "TE_POSTERIOR_LOW");
        assert(joint.final_qc.find("RESCUE") == std::string::npos);
    }

    {
        EventExistenceEvidence existence;
        existence.best_gt = "0/1";
        existence.af = 1.0;
        existence.gq = 89;
        existence.alt_struct_reads = 17;
        existence.ref_span_reads = 0;
        existence.depth = 17;
        existence.score = 3.0;

        EventSegmentationEvidence seg;
        seg.has_consensus = true;
        seg.has_left_flank = true;
        seg.has_right_flank = true;
        seg.has_insert_seq = true;
        seg.pair_valid = true;
        seg.left_align_len = 80;
        seg.right_align_len = 52;
        seg.left_identity = 0.95;
        seg.right_identity = 0.92;
        seg.insert_len = 236;
        seg.score = 1.04;
        seg.qc = "PASS_EVENT_SEGMENTATION";

        TEAlignmentEvidence te;
        te.best_family = "NA";
        te.best_subfamily = "NA";
        te.best_identity = 0.0;
        te.best_query_coverage = 0.0;
        te.cross_family_margin = 0.0;
        te.pass = false;
        te.qc_reason = "NO_TE_ALIGNMENT_MATCH";

        BoundaryEvidence boundary;
        boundary.geometry_defined = true;
        boundary.canonical_pass = true;
        boundary.evidence_consistent = true;
        boundary.boundary_type = "SMALL_DEL";
        boundary.boundary_len = 28;
        boundary.score = 1.0;
        boundary.qc = "PASS_BOUNDARY_SMALL_DEL";

        const auto joint = evaluate_joint_hypotheses(existence, seg, te, boundary);
        assert(!joint.emit_te_call);
    }

    {
        EventExistenceEvidence existence;
        existence.best_gt = "0/1";
        existence.af = 1.0;
        existence.gq = 89;
        existence.alt_struct_reads = 18;
        set_bam_breakdown(existence, 1, 0, 8, 9);
        existence.ref_span_reads = 0;
        existence.depth = 18;
        existence.score = 3.0;

        EventSegmentationEvidence seg;
        seg.has_consensus = true;
        seg.has_left_flank = false;
        seg.has_right_flank = false;
        seg.has_insert_seq = true;
        seg.pair_valid = true;
        seg.left_align_len = 0;
        seg.right_align_len = 0;
        seg.left_identity = 0.0;
        seg.right_identity = 0.0;
        seg.insert_len = 1719;
        seg.score = -3.0;
        seg.qc = "PASS_EVENT_SEGMENTATION_UNPLACED_INSERT";

        TEAlignmentEvidence te;
        te.best_family = "UNKNOWN";
        te.best_subfamily = "UNKNOWN";
        te.best_identity = 0.539192;
        te.best_query_coverage = 0.991274;
        te.cross_family_margin = 0.52518;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";
        te.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        te.sequence_model_score = 0.40;

        BoundaryEvidence boundary;
        boundary.geometry_defined = false;
        boundary.canonical_pass = false;
        boundary.evidence_consistent = false;
        boundary.boundary_type = "REJECT";
        boundary.boundary_len = 0;
        boundary.score = -2.0;
        boundary.qc = "REJECT_BOUNDARY_INVALID_REF_SEGMENTS";

        const auto joint = evaluate_joint_hypotheses(existence, seg, te, boundary);
        assert(joint.emit_te_call);
        assert(joint.emit_unknown_te);
        assert(joint.emit_evidence_te_call);
        assert(joint.final_qc == "PASS_TE_CLOSED" ||
               joint.final_qc == "PASS_TE_IMPRECISE");
        assert(joint.posterior_qc == "TE_POSTERIOR_LOW");
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
        EventSegmentation seg;
        seg.pass = true;
        seg.qc_reason = "PASS_EVENT_SEGMENTATION_ONE_SIDED_LEFT";
        seg.left_flank_align_len = 82;
        seg.right_flank_align_len = 0;
        seg.left_flank_identity = 0.96;
        seg.right_flank_identity = 0.0;
        seg.insert_seq = std::string(240, 'C');

        const auto evidence = analyze_event_segmentation_for_test(true, seg);
        assert(evidence.has_insert_seq);
        assert(evidence.has_left_flank);
        assert(!evidence.has_right_flank);
        assert(evidence.pair_valid);
        assert(evidence.score > 0.0);
    }

    {
        EventExistenceEvidence existence;
        existence.best_gt = "0/1";
        existence.af = 0.375;
        existence.gq = 99;
        existence.alt_struct_reads = 6;
        existence.ref_span_reads = 10;
        existence.depth = 16;
        existence.score = 3.0;

        EventSegmentation seg_input;
        seg_input.pass = true;
        seg_input.qc_reason = "PASS_EVENT_SEGMENTATION_ONE_SIDED_LEFT";
        seg_input.left_flank_align_len = 50;
        seg_input.right_flank_align_len = 0;
        seg_input.left_flank_identity = 0.96;
        seg_input.right_flank_identity = 0.0;
        seg_input.insert_seq = std::string(125, 'A');
        const auto seg = analyze_event_segmentation_for_test(true, seg_input);

        TEAlignmentEvidence te;
        te.best_family = "UNKNOWN";
        te.best_subfamily = "UNKNOWN";
        te.best_identity = 0.62406;
        te.best_query_coverage = 0.96;
        te.cross_family_margin = 0.0880977;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";

        BoundaryEvidence boundary;
        boundary.geometry_defined = false;
        boundary.canonical_pass = false;
        boundary.evidence_consistent = false;
        boundary.boundary_type = "REJECT";
        boundary.boundary_len = 0;
        boundary.score = -2.0;
        boundary.qc = "REJECT_BOUNDARY_INVALID_REF_SEGMENTS";

        const auto joint = evaluate_joint_hypotheses(existence, seg, te, boundary);
        // Ref-heavy / one-sided UNKNOWN: a TE hypothesis may now rank top as a
        // diagnostic (ranking is no longer gated by read/insert/seg ladders),
        // but the calibrated worst-case lFDR gate must abstain from a final call.
        assert(!joint.emit_te_call);
    }

    {
        EventExistenceEvidence existence;
        existence.best_gt = "0/1";
        existence.af = 0.2;
        existence.gq = 53;
        existence.alt_struct_reads = 2;
        existence.ref_span_reads = 8;
        existence.depth = 10;
        existence.score = 1.65;

        EventSegmentation seg_input;
        seg_input.pass = true;
        seg_input.qc_reason = "PASS_EVENT_SEGMENTATION_ONE_SIDED_LEFT";
        seg_input.left_flank_align_len = 80;
        seg_input.right_flank_align_len = 0;
        seg_input.left_flank_identity = 0.96;
        seg_input.right_flank_identity = 0.0;
        seg_input.insert_seq = std::string(1986, 'C');
        const auto seg = analyze_event_segmentation_for_test(true, seg_input);

        TEAlignmentEvidence te;
        te.best_family = "L2";
        te.best_subfamily = "NA";
        te.best_identity = 0.778575;
        te.best_query_coverage = 0.997986;
        te.cross_family_margin = 0.561628;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_FAMILY_ONLY";

        BoundaryEvidence boundary;
        boundary.geometry_defined = false;
        boundary.canonical_pass = false;
        boundary.evidence_consistent = false;
        boundary.boundary_type = "REJECT";
        boundary.boundary_len = 0;
        boundary.score = -2.0;
        boundary.qc = "REJECT_BOUNDARY_INVALID_REF_SEGMENTS";

        const auto joint = evaluate_joint_hypotheses(existence, seg, te, boundary);
        // Ref-heavy / one-sided UNKNOWN: a TE hypothesis may now rank top as a
        // diagnostic (ranking is no longer gated by read/insert/seg ladders),
        // but the calibrated worst-case lFDR gate must abstain from a final call.
        assert(!joint.emit_te_call);
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

    {
        EventExistenceEvidence existence;
        existence.best_gt = "1/1";
        existence.af = 1.0;
        existence.gq = 5;
        existence.alt_struct_reads = 2;
        existence.ref_span_reads = 0;
        existence.depth = 2;
        existence.score = -0.75;

        EventSegmentationEvidence seg;
        seg.has_consensus = true;
        seg.has_left_flank = true;
        seg.has_right_flank = true;
        seg.has_insert_seq = true;
        seg.pair_valid = true;
        seg.left_align_len = 80;
        seg.right_align_len = 80;
        seg.left_identity = 0.95;
        seg.right_identity = 0.95;
        seg.insert_len = 285;
        seg.score = 0.971429;
        seg.qc = "PASS_EVENT_SEGMENTATION";

        TEAlignmentEvidence te;
        te.best_family = "UNKNOWN";
        te.best_subfamily = "UNKNOWN";
        te.best_identity = 0.578313;
        te.best_query_coverage = 0.989474;
        te.cross_family_margin = 0.526612;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";

        BoundaryEvidence boundary;
        boundary.geometry_defined = true;
        boundary.canonical_pass = true;
        boundary.evidence_consistent = true;
        boundary.boundary_type = "BLUNT";
        boundary.boundary_len = 0;
        boundary.score = 1.0;
        boundary.qc = "PASS_BOUNDARY_BLUNT";

        const auto joint = evaluate_joint_hypotheses(existence, seg, te, boundary);
        assert(!joint.emit_te_call);
        assert(!joint.emit_unknown_te);
        assert(joint.emit_evidence_te_call);
        assert(joint.final_qc == "TE_AMBIGUOUS");
        assert(joint.final_qc.find("RESCUE") == std::string::npos);
    }

    {
        EventExistenceEvidence existence;
        existence.best_gt = "1/1";
        existence.af = 1.0;
        existence.gq = 36;
        existence.alt_struct_reads = 13;
        existence.ref_span_reads = 0;
        existence.depth = 13;
        existence.score = 0.8;

        EventSegmentation seg_input;
        seg_input.pass = true;
        seg_input.qc_reason = "PASS_EVENT_SEGMENTATION_ONE_SIDED_RIGHT";
        seg_input.left_flank_align_len = 0;
        seg_input.right_flank_align_len = 50;
        seg_input.left_flank_identity = 0.0;
        seg_input.right_flank_identity = 0.92;
        seg_input.insert_seq = std::string(283, 'T');

        const auto seg = analyze_event_segmentation_for_test(true, seg_input);
        if (seg.score >= 0.0) {
            return 102;
        }

        TEAlignmentEvidence te;
        te.best_family = "UNKNOWN";
        te.best_subfamily = "UNKNOWN";
        te.best_identity = 0.563863;
        te.best_query_coverage = 0.971731;
        te.cross_family_margin = 0.0115714;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";

        BoundaryEvidence boundary;
        boundary.geometry_defined = false;
        boundary.canonical_pass = false;
        boundary.evidence_consistent = false;
        boundary.boundary_type = "REJECT";
        boundary.boundary_len = 0;
        boundary.score = -2.0;
        boundary.qc = "REJECT_BOUNDARY_INVALID_REF_SEGMENTS";

        const auto joint = evaluate_joint_hypotheses(existence, seg, te, boundary);
        assert(!joint.emit_te_call);
        assert(!joint.emit_unknown_te);
        assert(joint.emit_evidence_te_call);
        assert(joint.final_qc == "TE_AMBIGUOUS");
        assert(joint.posterior_qc == "TE_POSTERIOR_LOW");
        assert(joint.final_qc.find("RESCUE") == std::string::npos);
    }

    {
        EventExistenceEvidence existence;
        existence.best_gt = "0/1";
        existence.af = 0.076923;
        existence.gq = 99;
        existence.alt_struct_reads = 1;
        set_bam_breakdown(existence, 1, 0, 0, 0);
        existence.ref_span_reads = 12;
        existence.depth = 13;
        existence.score = 3.0;

        EventSegmentationEvidence seg;
        seg.has_consensus = true;
        seg.has_left_flank = true;
        seg.has_right_flank = true;
        seg.has_insert_seq = true;
        seg.pair_valid = true;
        seg.left_align_len = 75;
        seg.right_align_len = 75;
        seg.left_identity = 0.95;
        seg.right_identity = 0.95;
        seg.insert_len = 125;
        seg.score = 1.225;
        seg.qc = "PASS_EVENT_SEGMENTATION";

        TEAlignmentEvidence te;
        te.best_family = "UNKNOWN";
        te.best_subfamily = "UNKNOWN";
        te.best_identity = 0.65;
        te.best_query_coverage = 0.99;
        te.cross_family_margin = 0.05;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";
        te.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        te.sequence_model_score = 0.40;
        te.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        te.sequence_model_score = 0.40;

        BoundaryEvidence boundary;
        boundary.geometry_defined = true;
        boundary.canonical_pass = true;
        boundary.evidence_consistent = true;
        boundary.boundary_type = "SMALL_DEL";
        boundary.boundary_len = 20;
        boundary.score = 1.0;
        boundary.qc = "PASS_BOUNDARY_SMALL_DEL";

        const auto joint = evaluate_joint_hypotheses(existence, seg, te, boundary);
        if (joint.emit_te_call) {
            return 104;
        }
    }

    {
        EventExistenceEvidence existence;
        existence.best_gt = "0/1";
        existence.af = 0.666667;
        existence.gq = 24;
        existence.alt_struct_reads = 2;
        existence.ref_span_reads = 1;
        existence.depth = 3;
        existence.score = 0.2;

        EventSegmentationEvidence seg;
        seg.has_consensus = true;
        seg.has_left_flank = true;
        seg.has_right_flank = true;
        seg.has_insert_seq = true;
        seg.pair_valid = true;
        seg.left_align_len = 77;
        seg.right_align_len = 63;
        seg.left_identity = 0.95;
        seg.right_identity = 0.94;
        seg.insert_len = 479;
        seg.score = 1.13;
        seg.qc = "PASS_EVENT_SEGMENTATION";

        TEAlignmentEvidence te;
        te.best_family = "Dong-R4";
        te.best_subfamily = "Dong-R4-1";
        te.best_identity = 0.717391;
        te.best_query_coverage = 0.223382;
        te.cross_family_margin = 0.133112;
        te.pass = false;
        te.qc_reason = "TE_ALIGNMENT_LOW_QUERY_COVERAGE";

        BoundaryEvidence boundary;
        boundary.geometry_defined = true;
        boundary.canonical_pass = true;
        boundary.evidence_consistent = true;
        boundary.boundary_type = "SMALL_DEL";
        boundary.boundary_len = 20;
        boundary.score = 1.0;
        boundary.qc = "PASS_BOUNDARY_SMALL_DEL";

        const auto joint = evaluate_joint_hypotheses(existence, seg, te, boundary);
        if (joint.emit_te_call) {
            return 105;
        }
    }

    {
        EventExistenceEvidence existence;
        existence.best_gt = "0/1";
        existence.af = 0.75;
        existence.gq = 99;
        existence.alt_struct_reads = 18;
        set_bam_breakdown(existence, 6, 4, 4, 4);
        existence.ref_span_reads = 6;
        existence.depth = 24;
        existence.score = 3.0;

        EventSegmentationEvidence seg;
        seg.has_consensus = true;
        seg.has_left_flank = true;
        seg.has_right_flank = true;
        seg.has_insert_seq = true;
        seg.pair_valid = true;
        seg.left_align_len = 80;
        seg.right_align_len = 80;
        seg.left_identity = 0.95;
        seg.right_identity = 0.95;
        seg.insert_len = 655;
        seg.score = 0.702857;
        seg.qc = "PASS_EVENT_SEGMENTATION";

        TEAlignmentEvidence te;
        te.best_family = "UNKNOWN";
        te.best_subfamily = "UNKNOWN";
        te.best_identity = 0.62;
        te.best_query_coverage = 0.98;
        te.cross_family_margin = 0.05;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";
        te.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        te.sequence_model_score = 0.40;

        BoundaryEvidence boundary;
        boundary.geometry_defined = true;
        boundary.canonical_pass = true;
        boundary.evidence_consistent = true;
        boundary.boundary_type = "BLUNT";
        boundary.boundary_len = 0;
        boundary.score = 1.0;
        boundary.qc = "PASS_BOUNDARY_BLUNT";

        const auto joint = evaluate_joint_hypotheses(existence, seg, te, boundary);
        if (!joint.emit_te_call ||
            !joint.emit_unknown_te ||
            joint.lfdr_qc != "TE_LFDR_HIGH") {
            return 106;
        }
    }

    {
        EventExistenceEvidence existence;
        existence.best_gt = "0/1";
        existence.af = 0.75;
        existence.gq = 99;
        existence.alt_struct_reads = 18;
        set_bam_breakdown(existence, 6, 4, 4, 4);
        existence.ref_span_reads = 6;
        existence.depth = 24;
        existence.score = 3.0;

        EventSegmentationEvidence seg;
        seg.has_consensus = true;
        seg.has_left_flank = true;
        seg.has_right_flank = true;
        seg.has_insert_seq = true;
        seg.pair_valid = true;
        seg.left_align_len = 80;
        seg.right_align_len = 80;
        seg.left_identity = 0.95;
        seg.right_identity = 0.95;
        seg.insert_len = 655;
        seg.score = 1.5;
        seg.qc = "PASS_EVENT_SEGMENTATION";

        TEAlignmentEvidence te;
        te.best_family = "UNKNOWN";
        te.best_subfamily = "UNKNOWN";
        te.best_identity = 0.559896;
        te.best_query_coverage = 0.99542;
        te.cross_family_margin = 0.00564167;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";
        te.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        te.sequence_model_score = 0.40;

        BoundaryEvidence boundary;
        boundary.geometry_defined = true;
        boundary.canonical_pass = true;
        boundary.evidence_consistent = true;
        boundary.boundary_type = "BLUNT";
        boundary.boundary_len = 0;
        boundary.score = 1.0;
        boundary.qc = "PASS_BOUNDARY_BLUNT";

        const auto joint = evaluate_joint_hypotheses(existence, seg, te, boundary);
        if (!joint.emit_te_call ||
            !joint.emit_unknown_te ||
            joint.lfdr_qc != "TE_LFDR_HIGH") {
            return 108;
        }
    }

    {
        EventExistenceEvidence existence;
        existence.best_gt = "0/1";
        existence.af = 0.357143;
        existence.gq = 99;
        existence.alt_struct_reads = 5;
        existence.ref_span_reads = 9;
        existence.depth = 14;
        existence.score = 3.0;

        EventSegmentationEvidence seg;
        seg.has_consensus = true;
        seg.has_left_flank = false;
        seg.has_right_flank = true;
        seg.has_insert_seq = true;
        seg.pair_valid = true;
        seg.left_align_len = 0;
        seg.right_align_len = 80;
        seg.left_identity = 0.0;
        seg.right_identity = 0.95;
        seg.insert_len = 180;
        seg.score = 1.0;
        seg.qc = "PASS_EVENT_SEGMENTATION_ONE_SIDED_RIGHT";

        TEAlignmentEvidence te;
        te.best_family = "UNKNOWN";
        te.best_subfamily = "UNKNOWN";
        te.best_identity = 0.60;
        te.best_query_coverage = 0.98;
        te.cross_family_margin = 0.05;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";

        BoundaryEvidence boundary;
        boundary.geometry_defined = false;
        boundary.canonical_pass = false;
        boundary.evidence_consistent = false;
        boundary.boundary_type = "REJECT";
        boundary.boundary_len = 0;
        boundary.score = -2.0;
        boundary.qc = "REJECT_BOUNDARY_INVALID_REF_SEGMENTS";

        const auto joint = evaluate_joint_hypotheses(existence, seg, te, boundary);
        if (joint.emit_te_call) {
            return 107;
        }
    }

    {
        EventExistenceEvidence existence;
        existence.best_gt = "0/1";
        existence.af = 0.222222;
        existence.gq = 0;
        existence.alt_struct_reads = 2;
        existence.ref_span_reads = 7;
        existence.depth = 9;
        existence.score = -1.0;

        EventSegmentationEvidence seg;
        seg.has_consensus = true;
        seg.has_left_flank = true;
        seg.has_right_flank = false;
        seg.has_insert_seq = true;
        seg.pair_valid = true;
        seg.left_align_len = 50;
        seg.right_align_len = 0;
        seg.left_identity = 0.95;
        seg.right_identity = 0.0;
        seg.insert_len = 381;
        seg.score = 0.52;
        seg.qc = "PASS_EVENT_SEGMENTATION_ONE_SIDED_LEFT";

        TEAlignmentEvidence te;
        te.best_family = "UNKNOWN";
        te.best_subfamily = "UNKNOWN";
        te.best_identity = 0.550538;
        te.best_query_coverage = 0.918635;
        te.cross_family_margin = 0.089011;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";

        BoundaryEvidence boundary;
        boundary.geometry_defined = false;
        boundary.canonical_pass = false;
        boundary.evidence_consistent = false;
        boundary.boundary_type = "REJECT";
        boundary.boundary_len = 0;
        boundary.score = -2.0;
        boundary.qc = "REJECT_BOUNDARY_INVALID_REF_SEGMENTS";

        const auto joint = evaluate_joint_hypotheses(existence, seg, te, boundary);
        if (joint.emit_te_call) {
            return 109;
        }
    }

    {
        EventExistenceEvidence existence;
        existence.best_gt = "0/1";
        existence.af = 0.777778;
        existence.gq = 36;
        existence.alt_struct_reads = 7;
        existence.ref_span_reads = 2;
        existence.depth = 9;
        existence.score = 0.8;

        EventSegmentationEvidence seg;
        seg.has_consensus = true;
        seg.has_left_flank = true;
        seg.has_right_flank = true;
        seg.has_insert_seq = true;
        seg.pair_valid = true;
        seg.left_align_len = 64;
        seg.right_align_len = 52;
        seg.left_identity = 0.95;
        seg.right_identity = 0.95;
        seg.insert_len = 327;
        seg.score = 0.655385;
        seg.qc = "PASS_EVENT_SEGMENTATION";

        TEAlignmentEvidence te;
        te.best_family = "Pao";
        te.best_subfamily = "Pao-7";
        te.best_identity = 0.548303;
        te.best_query_coverage = 0.896024;
        te.cross_family_margin = 0.00467437;
        te.pass = false;
        te.qc_reason = "TE_ALIGNMENT_LOW_IDENTITY";

        BoundaryEvidence boundary;
        boundary.geometry_defined = true;
        boundary.canonical_pass = true;
        boundary.evidence_consistent = true;
        boundary.boundary_type = "SMALL_DEL";
        boundary.boundary_len = 16;
        boundary.score = 1.0;
        boundary.qc = "PASS_BOUNDARY_SMALL_DEL";

        const auto joint = evaluate_joint_hypotheses(existence, seg, te, boundary);
        if (joint.emit_te_call) {
            return 110;
        }
    }

    {
        EventExistenceEvidence existence;
        existence.best_gt = "0/1";
        existence.af = 0.2;
        existence.gq = 0;
        existence.alt_struct_reads = 2;
        existence.ref_span_reads = 8;
        existence.depth = 10;
        existence.score = -1.0;

        EventSegmentationEvidence seg;
        seg.has_consensus = true;
        seg.has_left_flank = true;
        seg.has_right_flank = true;
        seg.has_insert_seq = true;
        seg.pair_valid = true;
        seg.left_align_len = 80;
        seg.right_align_len = 80;
        seg.left_identity = 0.95;
        seg.right_identity = 0.95;
        seg.insert_len = 497;
        seg.score = 0.832195;
        seg.qc = "PASS_EVENT_SEGMENTATION";

        TEAlignmentEvidence te;
        te.best_family = "UNKNOWN";
        te.best_subfamily = "UNKNOWN";
        te.best_identity = 0.541195;
        te.best_query_coverage = 0.981891;
        te.cross_family_margin = 0.499202;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";

        BoundaryEvidence boundary;
        boundary.geometry_defined = true;
        boundary.canonical_pass = true;
        boundary.evidence_consistent = true;
        boundary.boundary_type = "SMALL_DEL";
        boundary.boundary_len = 15;
        boundary.score = 1.0;
        boundary.qc = "PASS_BOUNDARY_SMALL_DEL";

        const auto joint = evaluate_joint_hypotheses(existence, seg, te, boundary);
        if (joint.emit_te_call) {
            return 111;
        }
    }

    {
        EventExistenceEvidence existence;
        existence.best_gt = "0/1";
        existence.af = 0.4;
        existence.gq = 99;
        existence.alt_struct_reads = 4;
        existence.ref_span_reads = 6;
        existence.depth = 10;
        existence.score = 3.0;

        EventSegmentationEvidence seg;
        seg.has_consensus = true;
        seg.has_left_flank = true;
        seg.has_right_flank = false;
        seg.has_insert_seq = true;
        seg.pair_valid = true;
        seg.left_align_len = 81;
        seg.right_align_len = 0;
        seg.left_identity = 0.95;
        seg.right_identity = 0.0;
        seg.insert_len = 360;
        seg.score = 1.12;
        seg.qc = "PASS_EVENT_SEGMENTATION_ONE_SIDED_LEFT";

        TEAlignmentEvidence te;
        te.best_family = "Rex-Babar";
        te.best_subfamily = "Rex-Babar-15";
        te.best_identity = 0.994444;
        te.best_query_coverage = 1.0;
        te.cross_family_margin = 0.952778;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT";

        BoundaryEvidence boundary;
        boundary.geometry_defined = false;
        boundary.canonical_pass = false;
        boundary.evidence_consistent = false;
        boundary.boundary_type = "REJECT";
        boundary.boundary_len = 0;
        boundary.score = -2.0;
        boundary.qc = "REJECT_BOUNDARY_INVALID_REF_SEGMENTS";

        const auto joint = evaluate_joint_hypotheses(existence, seg, te, boundary);
        if (joint.emit_te_call) {
            return 112;
        }
    }

    {
        EventExistenceEvidence existence;
        existence.best_gt = "0/1";
        existence.af = 0.533333;
        existence.gq = 99;
        existence.alt_struct_reads = 8;
        existence.ref_span_reads = 7;
        existence.depth = 15;
        existence.score = 3.0;

        EventSegmentationEvidence seg;
        seg.has_consensus = true;
        seg.has_left_flank = true;
        seg.has_right_flank = false;
        seg.has_insert_seq = true;
        seg.pair_valid = true;
        seg.left_align_len = 97;
        seg.right_align_len = 0;
        seg.left_identity = 0.95;
        seg.right_identity = 0.0;
        seg.insert_len = 4727;
        seg.score = 1.02763;
        seg.qc = "PASS_EVENT_SEGMENTATION_ONE_SIDED_LEFT";

        TEAlignmentEvidence te;
        te.best_family = "UNKNOWN";
        te.best_subfamily = "UNKNOWN";
        te.best_identity = 0.604518;
        te.best_query_coverage = 0.495875;
        te.cross_family_margin = 0.296169;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";

        BoundaryEvidence boundary;
        boundary.geometry_defined = false;
        boundary.canonical_pass = false;
        boundary.evidence_consistent = false;
        boundary.boundary_type = "REJECT";
        boundary.boundary_len = 0;
        boundary.score = -2.0;
        boundary.qc = "REJECT_BOUNDARY_INVALID_REF_SEGMENTS";

        const auto joint = evaluate_joint_hypotheses(existence, seg, te, boundary);
        if (joint.emit_te_call) {
            return 113;
        }
    }

    {
        EventExistenceEvidence existence;
        existence.best_gt = "0/1";
        existence.af = 0.227273;
        existence.gq = 0;
        existence.alt_struct_reads = 5;
        existence.ref_span_reads = 17;
        existence.depth = 22;
        existence.score = 0.0;

        EventSegmentationEvidence seg;
        seg.has_consensus = true;
        seg.has_left_flank = false;
        seg.has_right_flank = true;
        seg.has_insert_seq = true;
        seg.pair_valid = true;
        seg.left_align_len = 0;
        seg.right_align_len = 50;
        seg.left_identity = 0.0;
        seg.right_identity = 0.92;
        seg.insert_len = 881;
        seg.score = -0.0659155;
        seg.qc = "PASS_EVENT_SEGMENTATION_ONE_SIDED_RIGHT";

        TEAlignmentEvidence te;
        te.best_family = "UNKNOWN";
        te.best_subfamily = "UNKNOWN";
        te.best_identity = 0.530418;
        te.best_query_coverage = 0.994325;
        te.cross_family_margin = 0.305307;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";

        BoundaryEvidence boundary;
        boundary.geometry_defined = false;
        boundary.canonical_pass = false;
        boundary.evidence_consistent = false;
        boundary.boundary_type = "REJECT";
        boundary.boundary_len = 0;
        boundary.score = -2.0;
        boundary.qc = "REJECT_BOUNDARY_INVALID_REF_SEGMENTS";

        const auto joint = evaluate_joint_hypotheses(existence, seg, te, boundary);
        if (joint.emit_te_call) {
            return 114;
        }
    }

    {
        const auto existence = make_igv_existence(10, 0, 1, 3, 6, 2, 99, 3.0);
        const auto segmentation = make_igv_complete_segmentation(209, 1.44, 80, 72);
        const auto boundary = make_igv_boundary("SMALL_DEL", 8);
        const auto te = make_igv_te(
            "UNKNOWN",
            "UNKNOWN",
            0.615702,
            1.0,
            0.543932,
            true,
            "PASS_INSERT_TE_ALIGNMENT_UNKNOWN",
            "TE_MODEL_OUTLIER",
            -0.45);

        assert_igv_recall_emit(existence, segmentation, te, boundary, true);
    }

    {
        const auto existence = make_igv_existence(18, 0, 14, 1, 3, 10, 99, 3.0);
        const auto segmentation = make_igv_complete_segmentation(381, 1.60);
        const auto boundary = make_igv_boundary("BLUNT", 0);
        const auto te = make_igv_te(
            "UNKNOWN",
            "UNKNOWN",
            0.561881,
            0.952756,
            0.0199867,
            true,
            "PASS_INSERT_TE_ALIGNMENT_UNKNOWN",
            "TE_MODEL_EDGE",
            0.0);

        assert_igv_recall_emit(existence, segmentation, te, boundary, true);
    }

    {
        const auto existence = make_igv_existence(3, 0, 1, 1, 1, 7, 15, -0.25);
        const auto segmentation = make_igv_complete_segmentation(4859, 1.08, 85, 60);
        const auto boundary = make_igv_boundary("SMALL_DEL", 20);
        const auto te = make_igv_te(
            "UNKNOWN",
            "UNKNOWN",
            0.517575,
            0.415723,
            0.167919,
            false,
            "TE_ALIGNMENT_LOW_IDENTITY",
            "TE_MODEL_EDGE",
            0.0);

        assert_igv_recall_emit(existence, segmentation, te, boundary, true);
    }

    {
        const auto existence = make_igv_existence(9, 0, 3, 1, 5, 4, 99, 3.0);
        const auto segmentation = make_igv_complete_segmentation(108, 1.60);
        const auto boundary = make_igv_boundary("BLUNT", 0);
        const auto te = make_igv_te(
            "UNKNOWN",
            "UNKNOWN",
            0.598425,
            0.935185,
            0.41149,
            true,
            "PASS_INSERT_TE_ALIGNMENT_UNKNOWN",
            "TE_MODEL_OUTLIER",
            -0.45);

        assert_igv_recall_emit(existence, segmentation, te, boundary, true);
    }

    {
        const auto existence = make_igv_existence(1, 0, 1, 0, 0, 3, 31, 0.55);
        const auto segmentation = make_igv_complete_segmentation(491, 1.43, 78, 80);
        const auto boundary = make_igv_boundary("SMALL_DEL", 2);
        const auto te = make_igv_te(
            "Rex-Babar",
            "Rex-Babar-2",
            0.985537,
            0.985743,
            0.928717,
            true,
            "PASS_INSERT_TE_ALIGNMENT",
            "TE_MODEL_IN_DISTRIBUTION",
            0.40);

        assert_igv_recall_emit(existence, segmentation, te, boundary, false);
    }

    {
        const auto existence = make_igv_existence(6, 0, 1, 1, 4, 14, 86, 3.0);
        const auto segmentation = make_igv_complete_segmentation(244, 1.30, 73, 72);
        const auto boundary = make_igv_boundary("SMALL_DEL", 16);
        auto te = make_igv_te(
            "UNKNOWN",
            "UNKNOWN",
            0.902913,
            0.422131,
            0.0271013,
            true,
            "PASS_INSERT_TE_ALIGNMENT_UNKNOWN",
            "TE_MODEL_OUTLIER",
            -0.45);
        te.annotation_confidence = "LOW";
        te.annotation_residual_fraction = 0.106557;
        te.annotation_masked_fraction = 0.860656;

        assert_igv_recall_emit(existence, segmentation, te, boundary, true);
    }

    {
        const auto existence = make_igv_existence(6, 0, 4, 2, 0, 3, 86, 3.0);
        const auto segmentation = make_igv_complete_segmentation(297, 1.60);
        const auto boundary = make_igv_boundary("BLUNT", 0);
        const auto te = make_igv_te(
            "Unknown",
            "NA",
            0.764516,
            0.983165,
            0.292972,
            true,
            "PASS_INSERT_TE_ALIGNMENT_FAMILY_ONLY",
            "TE_MODEL_IN_DISTRIBUTION",
            0.40);

        assert_igv_recall_emit(existence, segmentation, te, boundary, false);
    }

    {
        const auto existence = make_igv_existence(9, 0, 3, 4, 2, 1, 83, 3.0);
        const auto segmentation =
            make_igv_complete_segmentation(408, 1.43179, 86, 78);
        const auto boundary = make_igv_boundary("TSD", 6);
        auto te = make_igv_te(
            "ERV",
            "NA",
            0.86692,
            0.612745,
            0.496887,
            true,
            "PASS_INSERT_TE_ALIGNMENT_FAMILY_ONLY",
            "TE_MODEL_IN_DISTRIBUTION",
            0.40);
        te.annotation_confidence = "MEDIUM";
        te.annotation_residual_fraction = 0.387255;
        te.annotation_masked_fraction = 0.0;

        assert_igv_recall_emit(existence, segmentation, te, boundary, false);
    }

    {
        const auto existence = make_igv_existence(90, 0, 0, 90, 42, 0, 99, 3.0);
        EventSegmentationEvidence segmentation;
        segmentation.has_consensus = true;
        segmentation.has_left_flank = false;
        segmentation.has_right_flank = false;
        segmentation.has_insert_seq = false;
        segmentation.pair_valid = false;
        segmentation.insert_len = 0;
        segmentation.score = -3.0;
        segmentation.qc = "NO_LEFT_FLANK_MATCH";

        TEAlignmentEvidence te;
        te.pass = false;
        te.qc_reason = "NO_TE_ALIGNMENT";
        te.sequence_model_label = "TE_MODEL_UNAVAILABLE";

        BoundaryEvidence boundary;
        boundary.geometry_defined = false;
        boundary.qc = "REJECT_BOUNDARY_UNSET";
        boundary.score = -3.0;

        const auto joint = evaluate_joint_hypotheses(
            existence,
            segmentation,
            te,
            boundary);
        assert(!joint.emit_te_call);
    }

    {
        const auto existence = make_igv_existence(8, 0, 0, 5, 3, 4, 20, 0.0);
        EventSegmentationEvidence segmentation;
        segmentation.has_consensus = true;
        segmentation.has_left_flank = false;
        segmentation.has_right_flank = false;
        segmentation.has_insert_seq = false;
        segmentation.pair_valid = false;
        segmentation.insert_len = 0;
        segmentation.score = -3.0;
        segmentation.qc = "NO_LEFT_FLANK_MATCH";

        TEAlignmentEvidence te;
        te.pass = false;
        te.qc_reason = "NO_TE_ALIGNMENT";
        te.sequence_model_label = "TE_MODEL_UNAVAILABLE";

        BoundaryEvidence boundary;
        boundary.geometry_defined = false;
        boundary.qc = "REJECT_BOUNDARY_UNSET";
        boundary.score = -3.0;

        const auto joint = evaluate_joint_hypotheses(
            existence,
            segmentation,
            te,
            boundary);
        assert(!joint.emit_te_call);
    }

    {
        EventExistenceEvidence existence;
        existence.best_gt = "0/1";
        existence.af = 4.0 / 35.0;
        existence.gq = 45;
        existence.alt_struct_reads = 4;
        set_bam_breakdown(existence, 0, 1, 1, 2);
        existence.ref_span_reads = 31;
        existence.depth = 35;
        existence.score = 1.2;

        EventSegmentationEvidence seg;
        seg.has_consensus = true;
        seg.has_left_flank = true;
        seg.has_right_flank = true;
        seg.has_insert_seq = true;
        seg.pair_valid = true;
        seg.left_align_len = 82;
        seg.right_align_len = 78;
        seg.left_identity = 0.95;
        seg.right_identity = 0.94;
        seg.insert_len = 289;
        seg.score = 1.05;
        seg.qc = "PASS_EVENT_SEGMENTATION";

        TEAlignmentEvidence te;
        te.best_family = "UNKNOWN";
        te.best_subfamily = "UNKNOWN";
        te.best_identity = 0.853448;
        te.best_query_coverage = 0.822695;
        te.cross_family_margin = 0.253448;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";
        te.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        te.sequence_model_score = 0.40;

        BoundaryEvidence boundary;
        boundary.geometry_defined = true;
        boundary.canonical_pass = true;
        boundary.evidence_consistent = true;
        boundary.boundary_type = "SMALL_DEL";
        boundary.boundary_len = 16;
        boundary.score = 1.0;
        boundary.qc = "PASS_BOUNDARY_SMALL_DEL";

        const auto joint = evaluate_joint_hypotheses(existence, seg, te, boundary);
        assert(!joint.emit_te_call);
        assert(joint.emit_structural_event_call);
        assert(joint.emit_unknown_te);
        assert(joint.emit_evidence_te_call);
        assert(joint.final_qc == "PASS_STRUCTURAL_INSERTION");
    }

    return 0;
}
