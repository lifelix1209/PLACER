#include "decision_policy.h"
#include "pipeline.h"

#include <cassert>
#include <string>
#include <vector>

int main() {
    using namespace placer;

    {
        InsertionEvidence e;
        e.tier = 1;
        e.support_reads = 1;
        InsertionAcceptanceParams p;
        p.emit_low_confidence_calls = false;
        const auto d = evaluate_insertion_acceptance(e, p);
        assert(d.pass);
        assert(!d.low_confidence);
    }

    {
        InsertionEvidence e;
        e.tier = 3;
        e.support_reads = 10;
        InsertionAcceptanceParams p;
        p.emit_low_confidence_calls = false;
        const auto d = evaluate_insertion_acceptance(e, p);
        assert(!d.pass);
    }

    {
        InsertionEvidence e;
        e.tier = 2;
        e.support_reads = 2;
        InsertionAcceptanceParams p;
        p.emit_low_confidence_calls = true;
        p.low_conf_min_support_reads = 3;
        p.low_conf_max_tier = 3;
        const auto d = evaluate_insertion_acceptance(e, p);
        assert(d.pass);
        assert(!d.low_confidence);
    }

    {
        InsertionEvidence e;
        e.tier = 3;
        e.support_reads = 5;
        InsertionAcceptanceParams p;
        p.emit_low_confidence_calls = true;
        p.low_conf_min_support_reads = 4;
        p.low_conf_max_tier = 3;
        const auto d = evaluate_insertion_acceptance(e, p);
        assert(d.pass);
        assert(d.low_confidence);
    }

    {
        InsertionEvidence e;
        e.tier = 3;
        e.support_reads = 2;
        InsertionAcceptanceParams p;
        p.emit_low_confidence_calls = true;
        p.low_conf_min_support_reads = 4;
        p.low_conf_max_tier = 3;
        const auto d = evaluate_insertion_acceptance(e, p);
        assert(!d.pass);
    }

    {
        TeOpenSetInput in;
        in.te_gate_pass = true;
        in.te_gate_uncertain_path = false;
        in.has_proxy_signal = false;
        in.proxy_certain = false;
        in.te_confidence_prob = 0.95;
        TeOpenSetParams p;
        const auto d = classify_te_open_set(in, p);
        assert(d.status == "TE_CERTAIN");
        assert(!d.promote_top1_name);
        assert(d.qc_proxy_tag.empty());
    }

    {
        TeOpenSetInput in;
        in.te_gate_pass = true;
        in.te_gate_uncertain_path = true;
        in.has_proxy_signal = true;
        in.proxy_certain = true;
        in.te_confidence_prob = 0.10;
        TeOpenSetParams p;
        const auto d = classify_te_open_set(in, p);
        assert(d.status == "TE_CERTAIN");
        assert(d.promote_top1_name);
        assert(d.qc_proxy_tag == "TE_PROXY_CERTAIN");
    }

    {
        TeOpenSetInput in;
        in.te_gate_pass = true;
        in.te_gate_uncertain_path = true;
        in.has_proxy_signal = true;
        in.proxy_certain = false;
        in.te_confidence_prob = 0.50;
        TeOpenSetParams p;
        p.conf_certain_min = 0.85;
        p.conf_uncertain_min = 0.35;
        const auto d = classify_te_open_set(in, p);
        assert(d.status == "TE_UNCERTAIN");
        assert(!d.promote_top1_name);
        assert(d.qc_proxy_tag == "TE_PROXY_WEAK");
    }

    {
        TeOpenSetInput in;
        in.te_gate_pass = true;
        in.te_gate_uncertain_path = true;
        in.has_proxy_signal = false;
        in.proxy_certain = false;
        in.te_confidence_prob = 0.90;
        TeOpenSetParams p;
        const auto d = classify_te_open_set(in, p);
        assert(d.status == "NON_TE");
        assert(d.qc_proxy_tag == "TE_PROXY_NONE");
    }

    {
        TeOpenSetInput in;
        in.te_gate_pass = false;
        in.te_gate_uncertain_path = false;
        in.has_proxy_signal = true;
        in.proxy_certain = false;
        in.te_confidence_prob = 0.20;
        TeOpenSetParams p;
        const auto d = classify_te_open_set(in, p);
        assert(d.add_te_gate_fail_tag);
    }

    {
        BreakpointConsistencyInput in;
        in.split_count = 2;
        in.clip_count = 6;
        in.indel_count = 0;
        in.all_breakpoint_mad = 72.0;
        in.split_breakpoint_mad = 6.0;
        in.split_like_breakpoint_mad = 6.0;
        in.clip_breakpoint_mad = 150.0;
        in.indel_breakpoint_mad = 0.0;
        in.split_clip_core_delta = 180.0;

        BreakpointConsistencyParams p;
        p.breakpoint_mad_max = 80.0;
        p.clip_breakpoint_mad_max = 120.0;
        p.split_clip_core_delta_max = 150.0;
        p.tight_split_like_breakpoint_mad_max = 40.0;
        p.min_split_like_count = 1;
        p.min_clip_count = 2;

        const auto d = evaluate_breakpoint_consistency(in, p);
        assert(d.pass_breakpoint_mad);
        assert(!d.pass_split_clip_consistency);
        assert(d.severe_inconsistency);
        assert(d.placeability_penalty >= 1.0);
    }

    {
        BreakpointConsistencyInput in;
        in.split_count = 0;
        in.clip_count = 4;
        in.indel_count = 3;
        in.all_breakpoint_mad = 44.0;
        in.split_breakpoint_mad = 0.0;
        in.split_like_breakpoint_mad = 8.0;
        in.clip_breakpoint_mad = 24.0;
        in.indel_breakpoint_mad = 7.0;
        in.split_clip_core_delta = 165.0;

        BreakpointConsistencyParams p;
        p.breakpoint_mad_max = 80.0;
        p.clip_breakpoint_mad_max = 120.0;
        p.split_clip_core_delta_max = 150.0;
        p.tight_split_like_breakpoint_mad_max = 40.0;
        p.min_split_like_count = 1;
        p.min_clip_count = 2;

        const auto d = evaluate_breakpoint_consistency(in, p);
        assert(d.pass_breakpoint_mad);
        assert(!d.pass_split_clip_consistency);
        assert(d.severe_inconsistency);
    }

    {
        BreakpointConsistencyInput in;
        in.split_count = 3;
        in.clip_count = 5;
        in.indel_count = 1;
        in.all_breakpoint_mad = 18.0;
        in.split_breakpoint_mad = 7.0;
        in.split_like_breakpoint_mad = 8.0;
        in.clip_breakpoint_mad = 20.0;
        in.indel_breakpoint_mad = 5.0;
        in.split_clip_core_delta = 12.0;

        BreakpointConsistencyParams p;
        p.breakpoint_mad_max = 80.0;
        p.clip_breakpoint_mad_max = 120.0;
        p.split_clip_core_delta_max = 150.0;
        p.tight_split_like_breakpoint_mad_max = 40.0;
        p.min_split_like_count = 1;
        p.min_clip_count = 2;

        const auto d = evaluate_breakpoint_consistency(in, p);
        assert(d.pass_breakpoint_mad);
        assert(d.pass_split_clip_consistency);
        assert(!d.severe_inconsistency);
        assert(d.placeability_penalty > 0.0);
        assert(d.placeability_penalty < 0.2);
    }

    {
        TeOpenSetInput in;
        in.te_gate_pass = true;
        in.force_non_te = true;
        in.has_proxy_signal = true;
        in.proxy_certain = true;
        in.te_confidence_prob = 0.99;
        TeOpenSetParams p;
        const auto d = classify_te_open_set(in, p);
        assert(d.status == "NON_TE");
    }

    {
        FinalCall call;
        call.te_qc = "PASS_CLASSIC|ANCHOR_FAIL_THETA_UNCERTAIN";
        call.tsd_type = "NONE";
        call.bp_source_counts = "split:1,clip:6,indel:0";
        call.te_top1_name = "L1:L1HS";
        call.te_top2_name = "L1:L1PA2";
        call.te_posterior_top2 = 0.31;
        call.te_posterior_margin = 0.08;
        call.te_confidence_prob = 0.93;

        PipelineConfig config;
        config.te_certain_posterior_margin_min = 0.20;
        config.te_confidence_anchor_fail_penalty = 0.15;
        config.te_confidence_no_tsd_penalty = 0.10;
        config.te_confidence_low_margin_penalty = 0.12;
        config.te_force_non_te_on_combined_weakness = true;
        config.te_force_non_te_on_anchor_weak_tsd = true;

        const auto d = apply_te_evidence_gates(call, config);
        assert(d.force_te_uncertain);
        assert(d.force_non_te);
        assert(call.te_confidence_prob <= 0.175);
        assert(call.te_qc.find("TE_FORCE_UNCERTAIN_ANCHOR") != std::string::npos);
        assert(call.te_qc.find("TE_FORCE_UNCERTAIN_NO_TSD") != std::string::npos);
        assert(call.te_qc.find("TE_FORCE_UNCERTAIN_LOW_MARGIN") != std::string::npos);
        assert(call.te_qc.find("TE_FORCE_NON_TE_ANCHOR_WEAK_TSD") != std::string::npos);
        assert(call.te_qc.find("TE_FORCE_NON_TE_COMBINED_WEAKNESS") != std::string::npos);
    }

    {
        FinalCall call;
        call.te_qc = "PASS_CLASSIC|ANCHOR_PASS";
        call.tsd_type = "NONE";
        call.bp_source_counts = "split:4,clip:2,indel:1";
        call.te_top1_name = "ALU:Ya5";
        call.te_top2_name = "ALU:Yb8";
        call.te_posterior_top2 = 0.12;
        call.te_posterior_margin = 0.30;
        call.te_confidence_prob = 0.88;

        PipelineConfig config;
        config.te_certain_posterior_margin_min = 0.20;

        const auto d = apply_te_evidence_gates(call, config);
        assert(!d.force_te_uncertain);
        assert(!d.force_non_te);
        assert(call.te_confidence_prob == 0.88);
        assert(call.te_qc.find("TE_FORCE_UNCERTAIN_") == std::string::npos);
    }

    {
        FinalCall call;
        call.te_qc = "PASS_CLASSIC|ANCHOR_PASS";
        call.tsd_type = "DEL";
        call.bp_source_counts = "split:5,clip:1,indel:1";
        call.te_top1_name = "Gypsy:Gypsy-95";
        call.te_top2_name = "Gypsy:Gypsy-40";
        call.te_posterior_top2 = 0.26;
        call.te_posterior_margin = 0.28;
        call.te_confidence_prob = 0.90;

        PipelineConfig config;
        config.te_certain_posterior_margin_min = 0.20;
        config.te_same_family_ambiguity_margin_max = 0.30;
        config.te_confidence_low_margin_penalty = 0.16;

        const auto d = apply_te_evidence_gates(call, config);
        assert(d.force_te_uncertain);
        assert(!d.force_non_te);
        assert(call.te_confidence_prob <= 0.60);
        assert(call.te_qc.find("TE_FORCE_UNCERTAIN_FAMILY_AMBIGUOUS") != std::string::npos);
    }

    {
        FinalCall call;
        call.te_qc = "PASS_CLASSIC|ANCHOR_NO_CORE_RESULT";
        call.tsd_type = "DEL";
        call.bp_source_counts = "split:2,clip:1,indel:1";
        call.te_top1_name = "ALU:Ya5";
        call.te_top2_name = "ALU:Yb8";
        call.te_posterior_top2 = 0.12;
        call.te_posterior_margin = 0.35;
        call.te_confidence_prob = 0.92;

        PipelineConfig config;
        config.te_confidence_prob_certain_min = 0.85;
        config.te_confidence_prob_uncertain_min = 0.35;
        config.te_confidence_anchor_fail_penalty = 0.15;
        config.te_force_non_te_on_anchor_weak_tsd = true;

        const auto d = apply_te_evidence_gates(call, config);
        assert(d.force_te_uncertain);
        assert(!d.force_non_te);
        assert(call.te_confidence_prob <= 0.60);
        assert(call.te_qc.find("TE_FORCE_UNCERTAIN_ANCHOR") != std::string::npos);
        assert(call.te_qc.find("TE_FORCE_NON_TE_ANCHOR_WEAK_TSD") == std::string::npos);
    }

    {
        FinalCall call;
        call.te_name = "L1:L1HS";
        call.te_median_identity = 0.01;
        call.te_qc = "PASS_SPLIT_INDEL_RESCUE|ANCHOR_PASS";
        call.tsd_type = "DEL";
        call.bp_source_counts = "split:6,clip:0,indel:4";
        call.te_confidence_prob = 0.91;

        PipelineConfig config;
        config.te_rescue_median_identity_min = 0.20;
        config.te_confidence_prob_uncertain_min = 0.35;

        const auto d = apply_te_evidence_gates(call, config);
        assert(!d.force_te_uncertain);
        assert(d.force_non_te);
        assert(call.te_qc.find("TE_FORCE_NON_TE_LOW_IDENTITY") != std::string::npos);
        assert(call.te_confidence_prob <= 0.175);
    }

    {
        FinalCall call;
        call.te_qc = "PASS_ASM_RESCUE|ANCHOR_PASS";
        call.tsd_type = "UNCERTAIN";
        call.bp_source_counts = "split:0,clip:5,indel:0";
        call.te_confidence_prob = 0.81;

        PipelineConfig config;
        config.te_confidence_prob_certain_min = 0.85;
        config.te_confidence_prob_uncertain_min = 0.35;
        config.te_confidence_no_tsd_penalty = 0.10;

        const auto d = apply_te_evidence_gates(call, config);
        assert(d.force_te_uncertain);
        assert(!d.force_non_te);
        assert(call.te_confidence_prob <= 0.60);
        assert(call.te_qc.find("TE_FORCE_UNCERTAIN_TSD_UNCERTAIN") != std::string::npos);
    }

    {
        ClusterTECall te_call;
        te_call.te_name = "Gypsy:Gypsy-95";
        te_call.fragment_count = 6;
        te_call.vote_fraction = 0.62;
        te_call.median_identity = 0.71;

        AssemblyCall assembly;
        assembly.qc_pass = true;
        assembly.identity_est = 0.92;
        assembly.consensus_len = 420;

        ComponentCall component;
        component.soft_clip_read_indices = {0, 1, 2, 3, 4, 5, 6, 7, 8};
        component.split_sa_read_indices = {9, 10, 11};
        component.breakpoint_candidates.push_back({"1", 1000, false, 0, 0, 0, "r0", 0, kCandidateSoftClip});
        component.breakpoint_candidates.push_back({"1", 1002, false, 0, 0, 0, "r1", 9, kCandidateSplitSaSupplementary});

        PipelineConfig config;
        config.te_min_fragments_for_vote = 2;
        config.te_vote_fraction_min = 0.40;
        config.te_median_identity_min = 0.30;
        config.te_mixed_min_non_softclip_reads = 4;
        config.te_mixed_min_non_softclip_frac = 0.35;
        config.te_mixed_clip_dominant_ratio = 2.0;

        const auto d = evaluate_post_assembly_te_decision(
            te_call,
            assembly,
            component,
            config);
        assert(d.pass);
        assert(d.force_te_uncertain);
        assert(d.qc == "PASS_INSERTION_TE_UNCERTAIN_CLIP_DOMINANT");
    }

    {
        ClusterTECall te_call;
        te_call.te_name = "Gypsy:Gypsy-95";
        te_call.fragment_count = 6;
        te_call.vote_fraction = 0.62;
        te_call.median_identity = 0.71;

        AssemblyCall assembly;
        assembly.qc_pass = true;
        assembly.identity_est = 0.92;
        assembly.consensus_len = 420;

        ComponentCall component;
        component.soft_clip_read_indices = {0, 1, 2, 3, 4, 5};
        component.split_sa_read_indices = {6, 7, 8, 9};
        component.breakpoint_candidates.push_back({"1", 1000, false, 0, 0, 0, "r0", 0, kCandidateSoftClip});
        component.breakpoint_candidates.push_back({"1", 1002, false, 0, 0, 0, "r1", 6, kCandidateSplitSaSupplementary});

        PipelineConfig config;
        config.te_min_fragments_for_vote = 2;
        config.te_vote_fraction_min = 0.40;
        config.te_median_identity_min = 0.30;
        config.te_mixed_min_non_softclip_reads = 4;
        config.te_mixed_min_non_softclip_frac = 0.35;
        config.te_mixed_clip_dominant_ratio = 2.0;

        const auto d = evaluate_post_assembly_te_decision(
            te_call,
            assembly,
            component,
            config);
        assert(d.pass);
        assert(!d.force_te_uncertain);
        assert(d.qc == "PASS_CLASSIC");
    }

    {
        TeOpenSetInput in;
        in.te_gate_pass = true;
        in.force_te_uncertain = true;
        in.has_proxy_signal = true;
        in.proxy_certain = true;
        in.te_confidence_prob = 0.99;
        TeOpenSetParams p;
        const auto d = classify_te_open_set(in, p);
        assert(d.status == "TE_UNCERTAIN");
    }

    {
        auto make_call = [](
                             int32_t pos,
                             const std::string& te_name,
                             int32_t support_reads,
                             double te_confidence_prob,
                             double te_posterior_margin,
                             const std::string& te_status) {
            FinalCall call;
            call.chrom = "1";
            call.tid = 0;
            call.pos = pos;
            call.te_name = te_name;
            call.support_reads = support_reads;
            call.te_confidence_prob = te_confidence_prob;
            call.te_posterior_margin = te_posterior_margin;
            call.te_status = te_status;
            call.confidence = "HIGH";
            return call;
        };

        PipelineResult result;
        result.final_calls.push_back(make_call(1000, "ALU:Ya5", 8, 0.82, 0.20, "TE_UNCERTAIN"));
        result.final_calls.push_back(make_call(1040, "L1:L1HS", 12, 0.91, 0.45, "TE_CERTAIN"));
        finalize_final_calls(result);
        assert(result.final_calls.size() == 1);
        assert(result.final_calls.front().pos == 1040);
        assert(result.final_calls.front().te_name == "L1:L1HS");
        assert(result.final_te_certain == 1);
        assert(result.final_te_uncertain == 0);
    }

    {
        auto make_call = [](
                             int32_t pos,
                             int32_t support_reads,
                             double te_confidence_prob) {
            FinalCall call;
            call.chrom = "1";
            call.tid = 0;
            call.pos = pos;
            call.te_name = "ALU:Ya5";
            call.support_reads = support_reads;
            call.te_confidence_prob = te_confidence_prob;
            call.te_status = "TE_CERTAIN";
            call.confidence = "HIGH";
            return call;
        };

        PipelineResult result;
        result.final_calls.push_back(make_call(1000, 5, 0.70));
        result.final_calls.push_back(make_call(1080, 9, 0.92));
        result.final_calls.push_back(make_call(1185, 7, 0.88));
        finalize_final_calls(result);
        assert(result.final_calls.size() == 2);
        assert(result.final_calls[0].pos == 1080);
        assert(result.final_calls[1].pos == 1185);
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

    return 0;
}
