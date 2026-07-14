#ifdef NDEBUG
#undef NDEBUG
#endif

#include "pipeline.h"

#include <cassert>
#include <string>
#include <vector>

namespace {

placer::FinalCall make_call(
    int32_t pos,
    int32_t alt_reads,
    double identity,
    double coverage,
    double margin,
    int32_t ref_span_reads) {
    placer::FinalCall call;
    call.chrom = "chr1";
    call.tid = 0;
    call.pos = pos;
    call.bp_left = pos;
    call.bp_right = pos;
    call.family = "L1";
    call.subfamily = "L1HS";
    call.family_committed = true;
    call.final_qc = "PASS_TE_CLOSED";
    call.alt_struct_reads = alt_reads;
    call.support_reads = alt_reads;
    call.ref_span_reads = ref_span_reads;
    call.max_raw_cigar_insert_len = 100;
    call.best_te_identity = identity;
    call.best_te_query_coverage = coverage;
    call.cross_family_margin = margin;
    return call;
}

placer::EvidenceLedgerRow make_null_control(
    int32_t pos,
    int32_t alt_reads,
    double identity,
    double coverage,
    double margin,
    int32_t ref_span_reads) {
    placer::EvidenceLedgerRow row;
    row.chrom = "chr1";
    row.tid = 0;
    row.pos = pos;
    row.bp_left = pos;
    row.bp_right = pos;
    row.final_qc = "REFERENCE_OR_ARTIFACT";
    row.alt_struct_reads = alt_reads;
    row.ref_span_reads = ref_span_reads;
    row.best_te_identity = identity;
    row.best_te_query_coverage = coverage;
    row.cross_family_margin = margin;
    row.mechanistic_blocks = "event;sequence";
    return row;
}

placer::EvidenceLedgerRow make_promotable_ledger_row(
    int32_t pos,
    int32_t alt_reads,
    double identity,
    double coverage,
    double margin,
    int32_t ref_span_reads) {
    placer::EvidenceLedgerRow row;
    row.chrom = "chr1";
    row.tid = 0;
    row.pos = pos;
    row.bp_left = pos;
    row.bp_right = pos;
    row.coverage_left = pos;
    row.coverage_right = pos;
    row.family = "L1";
    row.subfamily = "L1HS";
    row.family_alignment_resolved = true;
    row.final_qc = "TE_AMBIGUOUS";
    row.posterior_qc = "PASS_TE_POSTERIOR";
    row.lfdr_qc = "PASS_TE_LFDR";
    row.candidate_retention_reason = "EXPENSIVE_STAGE_CHALLENGER";
    row.alt_struct_reads = alt_reads;
    row.ref_span_reads = ref_span_reads;
    row.best_te_identity = identity;
    row.best_te_query_coverage = coverage;
    row.cross_family_margin = margin;
    row.te_posterior = 0.97;
    row.lfdr = 0.03;
    row.worst_case_lfdr = 0.04;
    row.mechanistic_lower_log_bf_te_vs_artifact = 5.0;
    row.mechanistic_lower_log_bf_te_vs_non_te = 4.0;
    row.mechanistic_ref_conflict_signal = 0.0;
    row.mechanistic_blocks = "event;sequence;boundary";
    row.robust_mechanistic_lfdr = 0.04;
    row.robust_mechanistic_worst_case_lfdr = 0.04;
    row.robust_mechanistic_qc = "PASS_TE_LFDR";
    return row;
}

std::vector<std::string> make_qnames(const std::string& prefix, int32_t count) {
    std::vector<std::string> out;
    out.reserve(static_cast<size_t>(std::max(0, count)));
    for (int32_t i = 0; i < count; ++i) {
        out.push_back(prefix + std::to_string(i));
    }
    return out;
}

placer::FinalCallFilterConfig te_calibrated_filter_config() {
    placer::FinalCallFilterConfig filter_config;
    filter_config.report_mode = placer::FinalReportMode::TeCalibrated;
    return filter_config;
}

void finalization_applies_sample_local_conformal_selector() {
    placer::PipelineResult result;
    result.final_calls.push_back(make_call(1000, 20, 0.98, 0.95, 0.20, 0));
    result.final_calls.push_back(make_call(2000, 3, 0.60, 0.50, 0.01, 6));
    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            3000 + i,
            i % 5,
            0.60,
            0.50,
            0.01,
            5));
    }

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().pos == 1000);
    assert(result.final_calls.front().conformal_qc == "PASS_CONFORMAL_FDR");
    assert(result.final_calls.front().conformal_null_p > 0.0);
    assert(result.final_calls.front().conformal_null_count == 40);
}

void finalization_preserves_legacy_dedup_when_no_ledger_is_available() {
    placer::PipelineResult result;
    result.final_calls.push_back(make_call(1000, 6, 0.90, 0.90, 0.10, 0));
    result.final_calls.push_back(make_call(1010, 10, 0.95, 0.95, 0.20, 0));

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().pos == 1010);
    assert(result.final_calls.front().conformal_qc == "CONFORMAL_NOT_EVALUATED");
}

void finalization_commits_sequence_certified_family_after_selection() {
    placer::PipelineResult result;
    auto call = make_call(1000, 10, 0.97, 0.94, 0.24, 0);
    call.family = "UNKNOWN";
    call.subfamily = "UNKNOWN";
    call.family_committed = false;
    call.te_name = "UNKNOWN";
    call.final_qc = "PASS_TE_IMPRECISE";
    call.sequence_family_candidate = "Gypsy";
    call.sequence_subfamily_candidate = "Gypsy-12";
    call.sequence_family_commit_eligible = true;
    result.final_calls.push_back(call);

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().family == "Gypsy");
    assert(result.final_calls.front().subfamily == "Gypsy-12");
    assert(result.final_calls.front().te_name == "Gypsy-12");
    assert(result.final_calls.front().family_committed);
    assert(result.final_calls.front().final_qc.find("SEQUENCE_FAMILY_COMMITTED") !=
           std::string::npos);
}

void finalization_keeps_pareto_weaker_distant_final_candidate() {
    placer::PipelineResult result;
    result.final_calls.push_back(make_call(1000, 20, 0.98, 0.95, 0.20, 0));
    result.final_calls.push_back(make_call(5000, 6, 0.80, 0.70, 0.10, 0));
    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            i % 3,
            0.50,
            0.50,
            0.01,
            5));
    }

    placer::finalize_final_calls(result);

    bool saw_first = false;
    bool saw_second = false;
    for (const auto& call : result.final_calls) {
        saw_first = saw_first || call.pos == 1000;
        saw_second = saw_second || call.pos == 5000;
        assert(call.conformal_qc == "PASS_CONFORMAL_FDR");
    }
    assert(result.final_calls.size() == 2);
    assert(saw_first);
    assert(saw_second);
}

void finalization_does_not_remove_pareto_weaker_distant_event() {
    placer::PipelineResult result;
    result.final_calls.push_back(make_call(1000, 20, 0.98, 0.95, 0.20, 0));
    result.final_calls.push_back(make_call(5000, 18, 0.96, 0.93, 0.18, 0));
    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            i % 5,
            0.50,
            0.50,
            0.01,
            5));
    }

    placer::finalize_final_calls(result);

    bool saw_first = false;
    bool saw_second = false;
    for (const auto& call : result.final_calls) {
        saw_first = saw_first || call.pos == 1000;
        saw_second = saw_second || call.pos == 5000;
    }
    assert(result.final_calls.size() == 2);
    assert(saw_first);
    assert(saw_second);
}

void finalization_promotes_event_cluster_from_strong_ledger_evidence() {
    placer::PipelineResult result;
    result.evidence_ledger.push_back(
        make_promotable_ledger_row(1000, 18, 0.96, 0.92, 0.18, 0));
    result.evidence_ledger.push_back(
        make_promotable_ledger_row(1015, 12, 0.91, 0.88, 0.12, 1));
    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            i % 3,
            0.50,
            0.50,
            0.01,
            5));
    }

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().pos == 1000);
    assert(result.final_calls.front().final_qc.find("EVENT_CLUSTER_PROMOTED") !=
           std::string::npos);
    assert(result.final_calls.front().conformal_qc == "PASS_CONFORMAL_FDR");
}

void finalization_recovers_pass_te_ledger_row_missing_from_final_calls() {
    placer::PipelineResult result;
    auto row = make_promotable_ledger_row(1000, 18, 0.96, 0.92, 0.18, 0);
    row.final_qc = "PASS_TE_CLOSED";
    row.left_flank_align_len = 76;
    row.right_flank_align_len = 82;
    row.insert_seq = "ACGTACGTACGT";
    result.evidence_ledger.push_back(row);
    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            i % 3,
            0.50,
            0.50,
            0.01,
            5));
    }

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().pos == 1000);
    assert(result.final_calls.front().final_qc.find("EVENT_CLUSTER_PROMOTED") !=
           std::string::npos);
    assert(result.final_calls.front().conformal_qc == "PASS_CONFORMAL_FDR");
    assert(result.final_calls.front().left_flank_align_len == 76);
    assert(result.final_calls.front().right_flank_align_len == 82);
    assert(result.final_calls.front().insert_seq == "ACGTACGTACGT");
}

void finalization_keeps_distinct_promoted_events_in_same_coverage_bin() {
    placer::PipelineResult result;
    auto first = make_promotable_ledger_row(1000, 18, 0.96, 0.92, 0.18, 0);
    first.coverage_left = 0;
    first.coverage_right = 9999;
    auto second = make_promotable_ledger_row(1800, 18, 0.96, 0.92, 0.18, 0);
    second.coverage_left = 0;
    second.coverage_right = 9999;
    result.evidence_ledger.push_back(first);
    result.evidence_ledger.push_back(second);
    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            i % 3,
            0.50,
            0.50,
            0.01,
            5));
    }

    placer::finalize_final_calls(result);

    bool saw_first = false;
    bool saw_second = false;
    for (const auto& call : result.final_calls) {
        saw_first = saw_first || call.pos == 1000;
        saw_second = saw_second || call.pos == 1800;
        assert(call.final_qc.find("EVENT_CLUSTER_PROMOTED") != std::string::npos);
    }
    assert(result.final_calls.size() == 2);
    assert(saw_first);
    assert(saw_second);
}

void finalization_does_not_promote_reference_conflict_without_statistical_support() {
    placer::PipelineResult result;
    auto conflict = make_promotable_ledger_row(1000, 18, 0.96, 0.92, 0.18, 12);
    conflict.mechanistic_ref_conflict_signal = 0.80;
    conflict.final_qc = "REFERENCE_OR_ARTIFACT";
    conflict.posterior_qc = "TE_POSTERIOR_LOW";
    conflict.lfdr_qc = "TE_LFDR_HIGH";
    conflict.te_posterior = 0.05;
    conflict.worst_case_lfdr = 0.95;
    conflict.robust_mechanistic_qc = "TE_LFDR_HIGH";
    conflict.robust_mechanistic_worst_case_lfdr = 0.95;
    conflict.mechanistic_lower_log_bf_te_vs_artifact = -3.0;
    conflict.mechanistic_lower_log_bf_te_vs_non_te = -2.0;
    conflict.mechanistic_ambiguity_width = 0.50;
    result.evidence_ledger.push_back(conflict);
    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            i % 3,
            0.50,
            0.50,
            0.01,
            5));
    }

    placer::finalize_final_calls(result);

    assert(result.final_calls.empty());
}

void finalization_prefers_stronger_promoted_cluster_over_weak_final_call() {
    placer::PipelineResult result;
    result.final_calls.push_back(make_call(1000, 2, 0.56, 0.52, 0.02, 8));
    auto strong = make_promotable_ledger_row(1010, 22, 0.97, 0.95, 0.22, 0);
    strong.final_qc = "PASS_TE_CLOSED";
    strong.posterior_qc = "PASS_TE_POSTERIOR";
    strong.lfdr_qc = "PASS_TE_LFDR";
    strong.robust_mechanistic_qc = "PASS_TE_LFDR";
    result.evidence_ledger.push_back(strong);
    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            i % 3,
            0.50,
            0.50,
            0.01,
            5));
    }

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().pos == 1010);
    assert(result.final_calls.front().final_qc.find("EVENT_CLUSTER_PROMOTED") !=
           std::string::npos);
}

void finalization_uses_event_lfdr_when_conformal_null_is_too_small() {
    placer::PipelineResult result;
    result.evidence_ledger.push_back(
        make_promotable_ledger_row(1000, 18, 0.96, 0.92, 0.18, 0));
    for (int i = 0; i < 3; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            i % 3,
            0.50,
            0.50,
            0.01,
            5));
    }

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().pos == 1000);
    assert(result.final_calls.front().conformal_qc == "PASS_EVENT_LFDR_FALLBACK");
}

void finalization_uses_event_ebh_when_null_is_too_small_but_bf_is_strong() {
    placer::PipelineResult result;
    auto strong = make_promotable_ledger_row(1000, 18, 0.80, 0.75, 0.12, 0);
    strong.final_qc = "PASS_TE_CLOSED";
    strong.posterior_qc = "TE_POSTERIOR_LOW";
    strong.lfdr_qc = "TE_LFDR_HIGH";
    strong.te_posterior = 0.55;
    strong.worst_case_lfdr = 0.95;
    strong.robust_mechanistic_qc = "TE_LFDR_HIGH";
    strong.robust_mechanistic_worst_case_lfdr = 0.35;
    strong.mechanistic_lower_log_bf_te_vs_artifact = 5.5;
    strong.mechanistic_lower_log_bf_te_vs_non_te = 4.2;
    strong.mechanistic_ambiguity_width = 0.30;
    result.evidence_ledger.push_back(strong);

    auto weak = make_promotable_ledger_row(5000, 18, 0.80, 0.75, 0.12, 0);
    weak.final_qc = "PASS_TE_CLOSED";
    weak.posterior_qc = "TE_POSTERIOR_LOW";
    weak.lfdr_qc = "TE_LFDR_HIGH";
    weak.te_posterior = 0.55;
    weak.worst_case_lfdr = 0.95;
    weak.robust_mechanistic_qc = "TE_LFDR_HIGH";
    weak.robust_mechanistic_worst_case_lfdr = 0.35;
    weak.mechanistic_lower_log_bf_te_vs_artifact = 1.0;
    weak.mechanistic_lower_log_bf_te_vs_non_te = 0.8;
    weak.mechanistic_ambiguity_width = 0.30;
    result.evidence_ledger.push_back(weak);

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().pos == 1000);
    assert(result.final_calls.front().conformal_qc == "PASS_EVENT_EBH");
}

void finalization_promotes_ambiguous_ledger_row_for_event_ebh_certificate() {
    placer::PipelineResult result;
    auto strong = make_promotable_ledger_row(1000, 22, 0.52, 0.99, 0.51, 0);
    strong.final_qc = "TE_AMBIGUOUS";
    strong.posterior_qc = "TE_POSTERIOR_LOW";
    strong.lfdr_qc = "TE_LFDR_HIGH";
    strong.te_posterior = 0.39;
    strong.lfdr = 0.84;
    strong.worst_case_lfdr = 0.98;
    strong.robust_mechanistic_qc = "TE_LFDR_HIGH";
    strong.robust_mechanistic_worst_case_lfdr = 0.52;
    strong.mechanistic_lower_log_bf_te_vs_artifact = 5.2;
    strong.mechanistic_lower_log_bf_te_vs_non_te = 3.5;
    strong.mechanistic_ambiguity_width = 0.30;
    result.evidence_ledger.push_back(strong);

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().pos == 1000);
    assert(result.final_calls.front().final_qc.find("EVENT_CLUSTER_PROMOTED") !=
           std::string::npos);
    assert(result.final_calls.front().conformal_qc == "PASS_EVENT_EBH");
}

void finalization_keeps_spatially_distinct_ebh_event_when_bfdr_event_exists() {
    placer::PipelineResult result;
    auto far_bfdr = make_promotable_ledger_row(1000, 4, 0.54, 0.99, 0.003, 0);
    far_bfdr.final_qc = "PASS_NONTE_INSERTION";
    far_bfdr.posterior_qc = "PASS_TE_POSTERIOR";
    far_bfdr.lfdr_qc = "TE_LFDR_HIGH";
    far_bfdr.te_posterior = 0.82;
    far_bfdr.lfdr = 0.94;
    far_bfdr.worst_case_lfdr = 0.99;
    far_bfdr.robust_mechanistic_qc = "TE_LFDR_HIGH";
    far_bfdr.robust_mechanistic_worst_case_lfdr = 0.87;
    far_bfdr.mechanistic_lower_log_bf_te_vs_artifact = 1.9;
    far_bfdr.mechanistic_lower_log_bf_te_vs_non_te = 1.7;
    far_bfdr.mechanistic_ambiguity_width = 0.39;
    result.evidence_ledger.push_back(far_bfdr);

    auto true_ebh = make_promotable_ledger_row(7000, 21, 0.52, 0.99, 0.51, 0);
    true_ebh.final_qc = "TE_AMBIGUOUS";
    true_ebh.posterior_qc = "TE_POSTERIOR_LOW";
    true_ebh.lfdr_qc = "TE_LFDR_HIGH";
    true_ebh.te_posterior = 0.39;
    true_ebh.lfdr = 0.83;
    true_ebh.worst_case_lfdr = 0.98;
    true_ebh.robust_mechanistic_qc = "TE_LFDR_HIGH";
    true_ebh.robust_mechanistic_worst_case_lfdr = 0.51;
    true_ebh.mechanistic_lower_log_bf_te_vs_artifact = 5.17;
    true_ebh.mechanistic_lower_log_bf_te_vs_non_te = 3.47;
    true_ebh.mechanistic_ambiguity_width = 0.30;
    result.evidence_ledger.push_back(true_ebh);

    placer::finalize_final_calls(result);

    bool saw_far = false;
    bool saw_true = false;
    for (const auto& call : result.final_calls) {
        if (call.pos == 1000) {
            saw_far = true;
        }
        if (call.pos == 7000) {
            saw_true = true;
            assert(call.conformal_qc == "PASS_EVENT_EBH");
        }
    }
    assert(saw_far);
    assert(saw_true);
}

void finalization_uses_event_bayesian_fdr_for_calibrated_lfdr_with_reference_span() {
    placer::PipelineResult result;
    auto strong = make_promotable_ledger_row(1000, 12, 0.86, 0.82, 0.12, 11);
    strong.mechanistic_ref_conflict_signal = 0.55;
    strong.final_qc = "TE_AMBIGUOUS";
    strong.posterior_qc = "PASS_TE_POSTERIOR";
    strong.lfdr_qc = "PASS_TE_LFDR";
    strong.te_posterior = 0.99;
    strong.worst_case_lfdr = 0.01;
    strong.robust_mechanistic_qc = "TE_LFDR_HIGH";
    strong.robust_mechanistic_worst_case_lfdr = 0.75;
    strong.mechanistic_lower_log_bf_te_vs_artifact = 1.0;
    strong.mechanistic_lower_log_bf_te_vs_non_te = 0.8;
    strong.mechanistic_ambiguity_width = 0.30;
    result.evidence_ledger.push_back(strong);

    auto weak = make_promotable_ledger_row(5000, 12, 0.86, 0.82, 0.12, 11);
    weak.mechanistic_ref_conflict_signal = 0.55;
    weak.final_qc = "TE_AMBIGUOUS";
    weak.posterior_qc = "TE_POSTERIOR_LOW";
    weak.lfdr_qc = "TE_LFDR_HIGH";
    weak.te_posterior = 0.20;
    weak.worst_case_lfdr = 0.45;
    weak.robust_mechanistic_qc = "TE_LFDR_HIGH";
    weak.robust_mechanistic_worst_case_lfdr = 0.75;
    weak.mechanistic_lower_log_bf_te_vs_artifact = 1.0;
    weak.mechanistic_lower_log_bf_te_vs_non_te = 0.8;
    weak.mechanistic_ambiguity_width = 0.30;
    result.evidence_ledger.push_back(weak);

    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            i % 3,
            0.50,
            0.50,
            0.01,
            5));
    }

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().pos == 1000);
    assert(result.final_calls.front().conformal_qc == "PASS_CONFORMAL_FDR" ||
           result.final_calls.front().conformal_qc == "PASS_EVENT_BFDR");
}

void finalization_uses_event_bayesian_fdr_for_pass_posterior_qc() {
    placer::PipelineResult result;
    auto strong = make_promotable_ledger_row(1000, 3, 0.80, 0.75, 0.12, 1);
    strong.final_qc = "PASS_TE_CLOSED";
    strong.posterior_qc = "PASS_TE_POSTERIOR";
    strong.lfdr_qc = "TE_LFDR_HIGH";
    strong.te_posterior = 0.82;
    strong.worst_case_lfdr = 0.95;
    strong.robust_mechanistic_qc = "TE_LFDR_HIGH";
    strong.robust_mechanistic_worst_case_lfdr = 0.62;
    strong.mechanistic_lower_log_bf_te_vs_artifact = 4.0;
    strong.mechanistic_lower_log_bf_te_vs_non_te = 3.0;
    strong.mechanistic_ambiguity_width = 0.30;
    result.evidence_ledger.push_back(strong);

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().pos == 1000);
    assert(result.final_calls.front().conformal_qc == "PASS_EVENT_BFDR" ||
           result.final_calls.front().conformal_qc == "PASS_EVENT_EBH");
}

void finalization_promotes_calibrated_near_complete_te_structure_despite_robust_uncertainty() {
    placer::PipelineResult result;
    auto row = make_promotable_ledger_row(1000, 4, 0.943038, 0.997845, 0.419673, 18);
    row.final_qc = "TE_AMBIGUOUS";
    row.posterior_qc = "PASS_TE_POSTERIOR";
    row.lfdr_qc = "PASS_TE_LFDR";
    row.te_posterior = 0.993225;
    row.non_te_posterior = 0.00370143;
    row.artifact_posterior = 0.00307329;
    row.lfdr = 0.03;
    row.worst_case_lfdr = 0.04;
    row.robust_mechanistic_qc = "TE_LFDR_HIGH";
    row.robust_mechanistic_worst_case_lfdr = 0.295719;
    row.raw_cigar_insert_reads = 15;
    row.max_raw_cigar_insert_len = 456;
    row.event_consensus_len = 608;
    row.te_structure_path =
        "0-463:TE_CORE:L1/UNKNOWN,463-464:UNEXPLAINED:HIGH_COMPLEXITY_RESIDUAL";
    row.te_structure_log_evidence = 4.07764;
    row.nonte_structure_log_evidence = -2.29372;
    row.artifact_structure_log_evidence = -1.75523;
    row.te_structure_path_confidence = 0.997079;
    row.te_core_coverage = 0.997845;
    row.unexplained_high_complexity_bp = 1;
    row.polyA_posterior = 0.0;
    row.transduction_posterior = 0.284391;
    row.mechanistic_lower_log_bf_te_vs_artifact = -1.2;
    row.mechanistic_lower_log_bf_te_vs_non_te = -0.4;
    row.mechanistic_ref_conflict_signal = 0.70;
    row.mechanistic_ambiguity_width = 0.30;
    row.mechanistic_blocks = "event;independent;sequence;structure;boundary;ref_conflict";
    result.evidence_ledger.push_back(row);

    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            i % 3,
            0.50,
            0.50,
            0.01,
            5));
    }

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().pos == 1000);
    assert(result.final_calls.front().final_qc.find("PASS_TE_IMPRECISE") == 0);
    assert(result.final_calls.front().conformal_qc == "PASS_EVENT_BFDR" ||
           result.final_calls.front().conformal_qc == "PASS_CONFORMAL_FDR");
}

void finalization_keeps_structural_insertion_with_event_existence_certificate() {
    placer::PipelineResult result;
    placer::FinalCall call;
    call.chrom = "chr1";
    call.tid = 0;
    call.pos = 1000;
    call.bp_left = 996;
    call.bp_right = 1005;
    call.family = "UNKNOWN";
    call.subfamily = "UNKNOWN";
    call.te_name = "UNKNOWN";
    call.final_qc = "PASS_STRUCTURAL_INSERTION";
    call.alt_struct_reads = 40;
    call.support_reads = 40;
    call.ref_span_reads = 1;
    call.best_te_identity = 0.655;
    call.best_te_query_coverage = 0.323;
    call.cross_family_margin = 0.02;
    call.non_te_posterior = 0.91;
    call.artifact_posterior = 0.04;
    call.posterior_qc = "TE_POSTERIOR_LOW";
    call.lfdr_qc = "TE_LFDR_HIGH";
    call.worst_case_lfdr = 0.96;
    call.mechanistic_lower_log_bf_te_vs_artifact = 4.8;
    call.mechanistic_lower_log_bf_te_vs_non_te = -2.0;
    call.mechanistic_ref_conflict_signal = 0.05;
    call.mechanistic_ambiguity_width = 0.30;
    call.mechanistic_blocks = "event;independent;boundary";
    call.robust_mechanistic_qc = "TE_LFDR_HIGH";
    call.robust_mechanistic_worst_case_lfdr = 0.80;
    result.final_calls.push_back(call);

    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            45,
            0.70,
            0.40,
            0.03,
            0));
    }

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().pos == 1000);
    assert(result.final_calls.front().final_qc == "PASS_STRUCTURAL_INSERTION");
    assert(result.final_calls.front().conformal_qc == "PASS_EVENT_EXISTENCE_BFDR");
}

void te_calibrated_report_mode_removes_structural_event_existence_call() {
    placer::PipelineResult result;
    placer::FinalCall call;
    call.chrom = "chr1";
    call.tid = 0;
    call.pos = 1000;
    call.bp_left = 996;
    call.bp_right = 1005;
    call.family = "UNKNOWN";
    call.subfamily = "UNKNOWN";
    call.te_name = "UNKNOWN";
    call.final_qc = "PASS_STRUCTURAL_INSERTION";
    call.alt_struct_reads = 40;
    call.support_reads = 40;
    call.ref_span_reads = 1;
    call.non_te_posterior = 0.91;
    call.artifact_posterior = 0.04;
    call.posterior_qc = "TE_POSTERIOR_LOW";
    call.lfdr_qc = "TE_LFDR_HIGH";
    call.mechanistic_lower_log_bf_te_vs_artifact = 4.8;
    call.mechanistic_lower_log_bf_te_vs_non_te = -2.0;
    call.mechanistic_ref_conflict_signal = 0.05;
    call.mechanistic_ambiguity_width = 0.30;
    call.mechanistic_blocks = "event;independent;boundary";
    result.final_calls.push_back(call);

    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            45,
            0.70,
            0.40,
            0.03,
            0));
    }

    placer::finalize_final_calls(result, 0.10, te_calibrated_filter_config());

    assert(result.final_calls.empty());
}

void te_calibrated_report_mode_keeps_calibrated_closed_te_call() {
    placer::PipelineResult result;
    auto call = make_call(1000, 24, 0.98, 0.95, 0.20, 0);
    call.final_qc = "PASS_TE_CLOSED";
    call.posterior_qc = "PASS_TE_POSTERIOR";
    call.lfdr_qc = "PASS_TE_LFDR";
    result.final_calls.push_back(call);

    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            i % 3,
            0.50,
            0.50,
            0.01,
            5));
    }

    placer::finalize_final_calls(result, 0.10, te_calibrated_filter_config());

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().pos == 1000);
    assert(result.final_calls.front().final_qc.find("PASS_TE") == 0);
    assert(result.final_calls.front().posterior_qc == "PASS_TE_POSTERIOR");
    assert(result.final_calls.front().lfdr_qc == "PASS_TE_LFDR");
    assert(result.final_calls.front().conformal_qc == "PASS_CONFORMAL_FDR");
}

void te_calibrated_report_mode_keeps_calibrated_promoted_te_call() {
    placer::PipelineResult result;
    auto row = make_promotable_ledger_row(1000, 18, 0.96, 0.92, 0.18, 0);
    row.final_qc = "PASS_TE_CLOSED";
    row.posterior_qc = "PASS_TE_POSTERIOR";
    row.lfdr_qc = "PASS_TE_LFDR";
    row.robust_mechanistic_qc = "PASS_TE_LFDR";
    result.evidence_ledger.push_back(row);

    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            i % 3,
            0.50,
            0.50,
            0.01,
            5));
    }

    placer::finalize_final_calls(result, 0.10, te_calibrated_filter_config());

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().pos == 1000);
    assert(result.final_calls.front().final_qc.find("PASS_TE_IMPRECISE") == 0);
    assert(result.final_calls.front().final_qc.find("EVENT_CLUSTER_PROMOTED") !=
           std::string::npos);
    assert(result.final_calls.front().posterior_qc == "PASS_TE_POSTERIOR");
    assert(result.final_calls.front().lfdr_qc == "PASS_TE_LFDR");
}

void te_calibrated_report_mode_removes_te_call_without_lfdr_certificate() {
    placer::PipelineResult result;
    auto call = make_call(1000, 24, 0.98, 0.95, 0.20, 0);
    call.final_qc = "PASS_TE_CLOSED";
    call.posterior_qc = "PASS_TE_POSTERIOR";
    call.lfdr_qc = "TE_LFDR_HIGH";
    result.final_calls.push_back(call);

    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            i % 3,
            0.50,
            0.50,
            0.01,
            5));
    }

    placer::finalize_final_calls(result, 0.10, te_calibrated_filter_config());

    assert(result.final_calls.empty());
}

void finalization_promotes_structural_ledger_event_without_te_calibration() {
    placer::PipelineResult result;
    auto structural = make_promotable_ledger_row(1000, 26, 0.52, 0.32, 0.02, 0);
    structural.family = "UNKNOWN";
    structural.subfamily = "UNKNOWN";
    structural.final_qc = "REFERENCE_OR_ARTIFACT";
    structural.posterior_qc = "TE_POSTERIOR_LOW";
    structural.lfdr_qc = "TE_LFDR_HIGH";
    structural.te_posterior = 0.08;
    structural.worst_case_lfdr = 0.96;
    structural.robust_mechanistic_qc = "TE_LFDR_HIGH";
    structural.robust_mechanistic_worst_case_lfdr = 0.80;
    structural.mechanistic_lower_log_bf_te_vs_artifact = 4.8;
    structural.mechanistic_lower_log_bf_te_vs_non_te = -2.0;
    structural.mechanistic_ambiguity_width = 0.30;
    structural.mechanistic_ref_conflict_signal = 0.0;
    structural.mechanistic_blocks = "event;independent;boundary";
    result.evidence_ledger.push_back(structural);

    auto distant_te_like = make_promotable_ledger_row(7000, 24, 0.85, 0.80, 0.20, 0);
    distant_te_like.final_qc = "TE_AMBIGUOUS";
    distant_te_like.posterior_qc = "TE_POSTERIOR_LOW";
    distant_te_like.lfdr_qc = "TE_LFDR_HIGH";
    distant_te_like.te_posterior = 0.45;
    distant_te_like.worst_case_lfdr = 0.90;
    distant_te_like.robust_mechanistic_qc = "TE_LFDR_HIGH";
    distant_te_like.robust_mechanistic_worst_case_lfdr = 0.70;
    distant_te_like.mechanistic_lower_log_bf_te_vs_artifact = 1.2;
    distant_te_like.mechanistic_lower_log_bf_te_vs_non_te = 0.8;
    distant_te_like.mechanistic_ambiguity_width = 0.30;
    result.evidence_ledger.push_back(distant_te_like);

    placer::finalize_final_calls(result);

    bool saw_structural = false;
    for (const auto& call : result.final_calls) {
        if (call.pos == 1000) {
            saw_structural = true;
            assert(call.final_qc.find("PASS_STRUCTURAL_INSERTION") == 0);
            assert(call.final_qc.find("EVENT_CLUSTER_PROMOTED") != std::string::npos);
            assert(call.family == "UNKNOWN");
            assert(call.subfamily == "UNKNOWN");
            assert(call.conformal_qc == "PASS_EVENT_EXISTENCE_BFDR");
        }
    }
    assert(saw_structural);
}

void finalization_promotes_balanced_ambiguous_ledger_event() {
    placer::PipelineResult result;
    auto balanced = make_promotable_ledger_row(1000, 24, 0.564655, 0.965458, 0.0175953, 13);
    balanced.family = "UNKNOWN";
    balanced.subfamily = "UNKNOWN";
    balanced.final_qc = "TE_AMBIGUOUS";
    balanced.posterior_qc = "TE_POSTERIOR_LOW";
    balanced.lfdr_qc = "TE_LFDR_HIGH";
    balanced.te_posterior = 0.180482;
    balanced.lfdr = 0.807638;
    balanced.worst_case_lfdr = 0.988411;
    balanced.mechanistic_lower_log_bf_te_vs_artifact = 1.02316;
    balanced.mechanistic_lower_log_bf_te_vs_non_te = 1.13198;
    balanced.mechanistic_ref_conflict_signal = 0.481853;
    balanced.mechanistic_ambiguity_width = 0.364489;
    balanced.robust_mechanistic_qc = "TE_LFDR_HIGH";
    balanced.robust_mechanistic_worst_case_lfdr = 0.928256;
    balanced.mechanistic_blocks = "event;independent;sequence;boundary";
    result.evidence_ledger.push_back(balanced);

    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            45,
            0.70,
            0.40,
            0.03,
            0));
    }

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().pos == 1000);
    assert(result.final_calls.front().final_qc.find("PASS_STRUCTURAL_INSERTION") == 0);
    assert(result.final_calls.front().final_qc.find("EVENT_CLUSTER_PROMOTED") !=
           std::string::npos);
    assert(result.final_calls.front().conformal_qc == "PASS_EVENT_EXISTENCE_BFDR");
}

void finalization_promotes_low_allele_fraction_event_with_error_model_certificate() {
    placer::PipelineResult result;
    auto low_af = make_promotable_ledger_row(1000, 7, 0.673684, 0.915441, 0.02, 35);
    low_af.family = "UNKNOWN";
    low_af.subfamily = "UNKNOWN";
    low_af.final_qc = "TE_AMBIGUOUS";
    low_af.posterior_qc = "TE_POSTERIOR_LOW";
    low_af.lfdr_qc = "TE_LFDR_HIGH";
    low_af.te_posterior = 0.12;
    low_af.lfdr = 0.90;
    low_af.worst_case_lfdr = 0.99;
    low_af.mechanistic_lower_log_bf_te_vs_artifact = -5.10;
    low_af.mechanistic_lower_log_bf_te_vs_non_te = -1.74;
    low_af.mechanistic_ref_conflict_signal = 0.83;
    low_af.mechanistic_ambiguity_width = 0.42;
    low_af.robust_mechanistic_qc = "TE_LFDR_HIGH";
    low_af.robust_mechanistic_worst_case_lfdr = 0.999;
    low_af.mechanistic_blocks = "event;independent;sequence;boundary;ref_conflict";
    low_af.support_qnames = make_qnames("low_af_support_", 7);
    result.evidence_ledger.push_back(low_af);

    auto corroborating = low_af;
    corroborating.pos = 1012;
    corroborating.bp_left = 1012;
    corroborating.bp_right = 1012;
    corroborating.coverage_left = 1012;
    corroborating.coverage_right = 1012;
    corroborating.alt_struct_reads = 7;
    corroborating.ref_span_reads = 36;
    corroborating.support_qnames = {
        "low_af_support_0",
        "low_af_support_1",
        "low_af_support_2",
        "low_af_support_3",
        "low_af_support_4",
        "low_af_support_5",
        "low_af_support_6",
    };
    result.evidence_ledger.push_back(corroborating);

    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            8,
            0.70,
            0.95,
            0.03,
            8));
    }

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().pos == 1000);
    assert(result.final_calls.front().family == "UNKNOWN");
    assert(result.final_calls.front().subfamily == "UNKNOWN");
    assert(result.final_calls.front().final_qc.find("PASS_STRUCTURAL_INSERTION") == 0);
    assert(result.final_calls.front().conformal_qc == "PASS_EVENT_EBH");
}

void finalization_rejects_isolated_low_allele_fraction_event_without_stable_community() {
    placer::PipelineResult result;
    auto low_af = make_promotable_ledger_row(1000, 7, 0.673684, 0.915441, 0.02, 35);
    low_af.family = "UNKNOWN";
    low_af.subfamily = "UNKNOWN";
    low_af.final_qc = "TE_AMBIGUOUS";
    low_af.posterior_qc = "TE_POSTERIOR_LOW";
    low_af.lfdr_qc = "TE_LFDR_HIGH";
    low_af.te_posterior = 0.12;
    low_af.lfdr = 0.90;
    low_af.worst_case_lfdr = 0.99;
    low_af.mechanistic_lower_log_bf_te_vs_artifact = -5.10;
    low_af.mechanistic_lower_log_bf_te_vs_non_te = -1.74;
    low_af.mechanistic_ref_conflict_signal = 0.83;
    low_af.mechanistic_ambiguity_width = 0.42;
    low_af.robust_mechanistic_qc = "TE_LFDR_HIGH";
    low_af.robust_mechanistic_worst_case_lfdr = 0.999;
    low_af.mechanistic_blocks = "event;independent;sequence;boundary;ref_conflict";
    low_af.support_qnames = make_qnames("isolated_low_af_support_", 7);
    result.evidence_ledger.push_back(low_af);

    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            8,
            0.70,
            0.95,
            0.03,
            8));
    }

    placer::finalize_final_calls(result);

    assert(result.final_calls.empty());
}

void finalization_promotes_low_af_event_with_self_stable_interval() {
    placer::PipelineResult result;
    auto interval = make_promotable_ledger_row(1057, 7, 0.540323, 0.973030, 0.228243, 26);
    interval.bp_left = 1000;
    interval.bp_right = 1114;
    interval.coverage_left = 1000;
    interval.coverage_right = 1114;
    interval.family = "UNKNOWN";
    interval.subfamily = "UNKNOWN";
    interval.final_qc = "TE_AMBIGUOUS";
    interval.posterior_qc = "TE_POSTERIOR_LOW";
    interval.lfdr_qc = "TE_LFDR_HIGH";
    interval.te_posterior = 0.0064;
    interval.lfdr = 0.95;
    interval.worst_case_lfdr = 0.997;
    interval.mechanistic_lower_log_bf_te_vs_artifact = -5.2;
    interval.mechanistic_lower_log_bf_te_vs_non_te = -1.8;
    interval.mechanistic_ref_conflict_signal = 0.79;
    interval.mechanistic_ambiguity_width = 0.40;
    interval.robust_mechanistic_qc = "TE_LFDR_HIGH";
    interval.robust_mechanistic_worst_case_lfdr = 0.997;
    interval.mechanistic_blocks = "event;independent;sequence;boundary;ref_conflict";
    interval.support_qnames = make_qnames("self_stable_interval_", 7);
    result.evidence_ledger.push_back(interval);

    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            8,
            0.70,
            0.95,
            0.03,
            8));
    }

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().family == "UNKNOWN");
    assert(result.final_calls.front().subfamily == "UNKNOWN");
    assert(result.final_calls.front().final_qc.find("PASS_STRUCTURAL_INSERTION") == 0);
    assert(result.final_calls.front().final_qc.find("EVENT_INTERVAL_STABLE") !=
           std::string::npos);
    assert(result.final_calls.front().conformal_qc == "PASS_EVENT_EBH");
}

void finalization_uses_low_af_evidence_not_fixed_read_floor_for_stable_interval() {
    placer::PipelineResult result;
    auto interval = make_promotable_ledger_row(1057, 4, 0.619718, 0.990566, 0.547834, 9);
    interval.bp_left = 1000;
    interval.bp_right = 1114;
    interval.coverage_left = 1000;
    interval.coverage_right = 1114;
    interval.family = "UNKNOWN";
    interval.subfamily = "UNKNOWN";
    interval.final_qc = "TE_AMBIGUOUS";
    interval.posterior_qc = "TE_POSTERIOR_LOW";
    interval.lfdr_qc = "TE_LFDR_HIGH";
    interval.te_posterior = 0.14;
    interval.lfdr = 0.90;
    interval.worst_case_lfdr = 0.99;
    interval.mechanistic_lower_log_bf_te_vs_artifact = -5.2;
    interval.mechanistic_lower_log_bf_te_vs_non_te = -1.8;
    interval.mechanistic_ref_conflict_signal = 0.69;
    interval.mechanistic_ambiguity_width = 0.40;
    interval.robust_mechanistic_qc = "TE_LFDR_HIGH";
    interval.robust_mechanistic_worst_case_lfdr = 0.99;
    interval.mechanistic_blocks = "event;independent;sequence;boundary;ref_conflict";
    interval.support_qnames = make_qnames("statistical_interval_", 4);
    result.evidence_ledger.push_back(interval);

    auto insufficient = interval;
    insufficient.pos = 3000;
    insufficient.bp_left = 2960;
    insufficient.bp_right = 3060;
    insufficient.coverage_left = 2960;
    insufficient.coverage_right = 3060;
    insufficient.alt_struct_reads = 3;
    insufficient.ref_span_reads = 9;
    insufficient.support_qnames = make_qnames("insufficient_interval_", 3);
    result.evidence_ledger.push_back(insufficient);

    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            8,
            0.70,
            0.95,
            0.03,
            8));
    }

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().pos == 1057);
    assert(result.final_calls.front().family == "UNKNOWN");
    assert(result.final_calls.front().subfamily == "UNKNOWN");
    assert(result.final_calls.front().final_qc.find("EVENT_INTERVAL_STABLE") !=
           std::string::npos);
    assert(result.final_calls.front().conformal_qc == "PASS_EVENT_EBH");
}

void finalization_does_not_collapse_nonoverlapping_stable_interval_into_adjacent_te_call() {
    placer::PipelineResult result;

    auto interval = make_promotable_ledger_row(1057, 4, 0.619718, 0.990566, 0.547834, 9);
    interval.bp_left = 1000;
    interval.bp_right = 1114;
    interval.coverage_left = 1000;
    interval.coverage_right = 1114;
    interval.family = "UNKNOWN";
    interval.subfamily = "UNKNOWN";
    interval.final_qc = "TE_AMBIGUOUS";
    interval.posterior_qc = "TE_POSTERIOR_LOW";
    interval.lfdr_qc = "TE_LFDR_HIGH";
    interval.te_posterior = 0.141246;
    interval.lfdr = 0.474954;
    interval.worst_case_lfdr = 0.993724;
    interval.mechanistic_lower_log_bf_te_vs_artifact = -3.08068;
    interval.mechanistic_lower_log_bf_te_vs_non_te = -0.633782;
    interval.mechanistic_ref_conflict_signal = 0.692308;
    interval.mechanistic_ambiguity_width = 0.453553;
    interval.robust_mechanistic_qc = "TE_LFDR_HIGH";
    interval.robust_mechanistic_worst_case_lfdr = 0.998834;
    interval.mechanistic_blocks = "event;independent;sequence;boundary;ref_conflict";
    interval.support_qnames = make_qnames("stable_interval_left_", 4);
    result.evidence_ledger.push_back(interval);

    auto adjacent_te = make_promotable_ledger_row(1218, 5, 0.963080, 0.655785, 0.628290, 9);
    adjacent_te.bp_left = 1162;
    adjacent_te.bp_right = 1275;
    adjacent_te.coverage_left = 1162;
    adjacent_te.coverage_right = 1275;
    adjacent_te.family = "Rex-Babar";
    adjacent_te.subfamily = "NA";
    adjacent_te.final_qc = "TE_AMBIGUOUS";
    adjacent_te.posterior_qc = "TE_POSTERIOR_LOW";
    adjacent_te.lfdr_qc = "PASS_TE_LFDR";
    adjacent_te.te_posterior = 0.357850;
    adjacent_te.lfdr = 0.000203884;
    adjacent_te.worst_case_lfdr = 0.0827198;
    adjacent_te.mechanistic_lower_log_bf_te_vs_artifact = -3.64502;
    adjacent_te.mechanistic_lower_log_bf_te_vs_non_te = -1.06423;
    adjacent_te.mechanistic_ref_conflict_signal = 0.642857;
    adjacent_te.mechanistic_ambiguity_width = 0.470716;
    adjacent_te.robust_mechanistic_qc = "TE_LFDR_HIGH";
    adjacent_te.robust_mechanistic_worst_case_lfdr = 0.999348;
    adjacent_te.mechanistic_blocks = "event;independent;sequence;boundary;ref_conflict";
    adjacent_te.support_qnames = make_qnames("adjacent_te_right_", 5);
    result.evidence_ledger.push_back(adjacent_te);

    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            8,
            0.70,
            0.95,
            0.03,
            8));
    }

    placer::finalize_final_calls(result);

    bool saw_stable_interval = false;
    bool saw_adjacent_te = false;
    for (const auto& call : result.final_calls) {
        if (call.pos == 1057) {
            saw_stable_interval = true;
            assert(call.family == "UNKNOWN");
            assert(call.subfamily == "UNKNOWN");
            assert(call.final_qc.find("PASS_STRUCTURAL_INSERTION") == 0);
            assert(call.final_qc.find("EVENT_INTERVAL_STABLE") != std::string::npos);
            assert(call.conformal_qc == "PASS_EVENT_EBH");
        }
        if (call.pos == 1218) {
            saw_adjacent_te = true;
        }
    }
    assert(saw_stable_interval);
    assert(saw_adjacent_te);
}

void finalization_promotes_low_af_event_with_subthreshold_stable_neighbor() {
    placer::PipelineResult result;
    auto focal = make_promotable_ledger_row(1000, 7, 0.673684, 0.915441, 0.02, 35);
    focal.family = "UNKNOWN";
    focal.subfamily = "UNKNOWN";
    focal.final_qc = "TE_AMBIGUOUS";
    focal.posterior_qc = "TE_POSTERIOR_LOW";
    focal.lfdr_qc = "TE_LFDR_HIGH";
    focal.te_posterior = 0.12;
    focal.lfdr = 0.90;
    focal.worst_case_lfdr = 0.99;
    focal.mechanistic_lower_log_bf_te_vs_artifact = -5.10;
    focal.mechanistic_lower_log_bf_te_vs_non_te = -1.74;
    focal.mechanistic_ref_conflict_signal = 0.83;
    focal.mechanistic_ambiguity_width = 0.42;
    focal.robust_mechanistic_qc = "TE_LFDR_HIGH";
    focal.robust_mechanistic_worst_case_lfdr = 0.999;
    focal.mechanistic_blocks = "event;independent;sequence;boundary;ref_conflict";
    focal.support_qnames = make_qnames("stable_subthreshold_support_", 7);
    result.evidence_ledger.push_back(focal);

    auto subthreshold_neighbor = focal;
    subthreshold_neighbor.pos = 1012;
    subthreshold_neighbor.bp_left = 1004;
    subthreshold_neighbor.bp_right = 1012;
    subthreshold_neighbor.coverage_left = 1004;
    subthreshold_neighbor.coverage_right = 1012;
    subthreshold_neighbor.alt_struct_reads = 4;
    subthreshold_neighbor.ref_span_reads = 36;
    subthreshold_neighbor.best_te_identity = 0.72;
    subthreshold_neighbor.best_te_query_coverage = 0.86;
    subthreshold_neighbor.cross_family_margin = 0.004;
    subthreshold_neighbor.support_qnames = {
        "stable_subthreshold_support_0",
        "stable_subthreshold_support_1",
        "stable_subthreshold_support_2",
        "stable_subthreshold_support_3",
    };
    result.evidence_ledger.push_back(subthreshold_neighbor);

    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            8,
            0.70,
            0.95,
            0.03,
            8));
    }

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().final_qc.find("EVENT_COMMUNITY_STABLE") !=
           std::string::npos);
    assert(result.final_calls.front().conformal_qc == "PASS_EVENT_EBH");
}

void finalization_prefers_stable_low_af_community_over_uncertified_direct_call() {
    placer::PipelineResult result;
    placer::FinalCall direct;
    direct.chrom = "chr1";
    direct.tid = 0;
    direct.pos = 1000;
    direct.bp_left = 1000;
    direct.bp_right = 1000;
    direct.family = "UNKNOWN";
    direct.subfamily = "UNKNOWN";
    direct.te_name = "UNKNOWN";
    direct.final_qc = "PASS_STRUCTURAL_INSERTION";
    direct.alt_struct_reads = 7;
    direct.support_reads = 7;
    direct.ref_span_reads = 35;
    direct.best_te_identity = 0.673684;
    direct.best_te_query_coverage = 0.915441;
    direct.cross_family_margin = 0.02;
    direct.mechanistic_lower_log_bf_te_vs_artifact = -5.10;
    direct.mechanistic_lower_log_bf_te_vs_non_te = -1.74;
    direct.mechanistic_ref_conflict_signal = 0.83;
    direct.mechanistic_ambiguity_width = 0.42;
    direct.mechanistic_blocks = "event;independent;sequence;boundary;ref_conflict";
    result.final_calls.push_back(direct);

    auto low_af = make_promotable_ledger_row(1000, 7, 0.673684, 0.915441, 0.02, 35);
    low_af.family = "UNKNOWN";
    low_af.subfamily = "UNKNOWN";
    low_af.final_qc = "TE_AMBIGUOUS";
    low_af.posterior_qc = "TE_POSTERIOR_LOW";
    low_af.lfdr_qc = "TE_LFDR_HIGH";
    low_af.te_posterior = 0.12;
    low_af.lfdr = 0.90;
    low_af.worst_case_lfdr = 0.99;
    low_af.mechanistic_lower_log_bf_te_vs_artifact = -5.10;
    low_af.mechanistic_lower_log_bf_te_vs_non_te = -1.74;
    low_af.mechanistic_ref_conflict_signal = 0.83;
    low_af.mechanistic_ambiguity_width = 0.42;
    low_af.robust_mechanistic_qc = "TE_LFDR_HIGH";
    low_af.robust_mechanistic_worst_case_lfdr = 0.999;
    low_af.mechanistic_blocks = "event;independent;sequence;boundary;ref_conflict";
    low_af.support_qnames = make_qnames("community_support_", 7);
    result.evidence_ledger.push_back(low_af);

    auto corroborating = low_af;
    corroborating.pos = 1012;
    corroborating.bp_left = 1012;
    corroborating.bp_right = 1012;
    corroborating.coverage_left = 1012;
    corroborating.coverage_right = 1012;
    result.evidence_ledger.push_back(corroborating);

    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            8,
            0.70,
            0.95,
            0.03,
            8));
    }

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().final_qc.find("EVENT_COMMUNITY_STABLE") !=
           std::string::npos);
    assert(result.final_calls.front().conformal_qc == "PASS_EVENT_EBH");
}

void finalization_aggregates_connected_weak_event_community() {
    placer::PipelineResult result;

    auto weak = make_promotable_ledger_row(1000, 3, 0.0, 0.0, 0.0, 24);
    weak.family = "UNKNOWN";
    weak.subfamily = "UNKNOWN";
    weak.final_qc = "REFERENCE_OR_ARTIFACT";
    weak.posterior_qc = "TE_POSTERIOR_LOW";
    weak.lfdr_qc = "TE_LFDR_HIGH";
    weak.te_posterior = 0.08;
    weak.lfdr = 0.92;
    weak.worst_case_lfdr = 0.99;
    weak.mechanistic_lower_log_bf_te_vs_artifact = -4.4;
    weak.mechanistic_lower_log_bf_te_vs_non_te = -2.1;
    weak.mechanistic_ref_conflict_signal = 0.86;
    weak.mechanistic_ambiguity_width = 0.40;
    weak.robust_mechanistic_qc = "TE_LFDR_HIGH";
    weak.robust_mechanistic_worst_case_lfdr = 0.99;
    weak.mechanistic_blocks = "event;independent;boundary;ref_conflict";
    weak.support_qnames = {"weak_read_0", "weak_read_1", "weak_read_2"};
    result.evidence_ledger.push_back(weak);

    auto bridge = weak;
    bridge.pos = 1010;
    bridge.bp_left = 1010;
    bridge.bp_right = 1010;
    bridge.coverage_left = 1010;
    bridge.coverage_right = 1010;
    bridge.alt_struct_reads = 4;
    bridge.ref_span_reads = 25;
    bridge.support_qnames = {
        "weak_read_1",
        "weak_read_2",
        "weak_read_3",
        "weak_read_4",
    };
    result.evidence_ledger.push_back(bridge);

    auto tail = weak;
    tail.pos = 1020;
    tail.bp_left = 1020;
    tail.bp_right = 1020;
    tail.coverage_left = 1020;
    tail.coverage_right = 1020;
    tail.alt_struct_reads = 4;
    tail.ref_span_reads = 26;
    tail.support_qnames = {
        "weak_read_3",
        "weak_read_4",
        "weak_read_5",
        "weak_read_6",
    };
    result.evidence_ledger.push_back(tail);

    auto nearby_disconnected = weak;
    nearby_disconnected.pos = 1030;
    nearby_disconnected.bp_left = 1030;
    nearby_disconnected.bp_right = 1030;
    nearby_disconnected.coverage_left = 1030;
    nearby_disconnected.coverage_right = 1030;
    nearby_disconnected.alt_struct_reads = 3;
    nearby_disconnected.ref_span_reads = 25;
    nearby_disconnected.support_qnames = {
        "nearby_noise_0",
        "nearby_noise_1",
        "nearby_noise_2",
    };
    result.evidence_ledger.push_back(nearby_disconnected);

    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            8,
            0.70,
            0.95,
            0.03,
            8));
    }

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().family == "UNKNOWN");
    assert(result.final_calls.front().subfamily == "UNKNOWN");
    assert(result.final_calls.front().final_qc.find("PASS_STRUCTURAL_INSERTION") == 0);
    assert(result.final_calls.front().final_qc.find("EVENT_COMMUNITY_AGGREGATED") !=
           std::string::npos);
    assert(result.final_calls.front().final_qc.find("EVENT_COMMUNITY_STABLE") !=
           std::string::npos);
    assert(result.final_calls.front().alt_struct_reads == 7);
    assert(result.final_calls.front().support_reads == 7);
    assert(result.final_calls.front().ref_span_reads == 25);
    assert(result.final_calls.front().bp_left == 950);
    assert(result.final_calls.front().bp_right == 1070);
    assert(result.final_calls.front().conformal_qc == "PASS_EVENT_EBH");
}

void finalization_competes_overlapping_promoted_event_envelopes() {
    placer::PipelineResult result;

    placer::FinalCall community;
    community.chrom = "chr1";
    community.tid = 0;
    community.pos = 1120;
    community.bp_left = 1000;
    community.bp_right = 1200;
    community.family = "UNKNOWN";
    community.subfamily = "UNKNOWN";
    community.te_name = "UNKNOWN";
    community.final_qc =
        "PASS_STRUCTURAL_INSERTION|EVENT_CLUSTER_PROMOTED|"
        "EVENT_COMMUNITY_STABLE|EVENT_COMMUNITY_AGGREGATED";
    community.best_explanation = "EVENT_COMMUNITY";
    community.explanation_path = "READ_OVERLAP_COMMUNITY";
    community.support_reads = 4;
    community.alt_struct_reads = 4;
    community.ref_span_reads = 10;
    community.event_consensus_len = 156;
    community.non_te_posterior = 0.76;
    community.artifact_posterior = 0.23;
    community.mechanistic_blocks = "event;community;read_overlap;boundary;ref_conflict";
    result.final_calls.push_back(community);

    placer::FinalCall imprecise_te;
    imprecise_te.chrom = "chr1";
    imprecise_te.tid = 0;
    imprecise_te.pos = 1450;
    imprecise_te.bp_left = 1450;
    imprecise_te.bp_right = 1450;
    imprecise_te.family = "L2";
    imprecise_te.subfamily = "L2";
    imprecise_te.te_name = "L2";
    imprecise_te.final_qc = "PASS_TE_IMPRECISE|EVENT_CLUSTER_PROMOTED";
    imprecise_te.best_explanation = "EVENT_CLUSTER";
    imprecise_te.support_reads = 3;
    imprecise_te.alt_struct_reads = 3;
    imprecise_te.ref_span_reads = 6;
    imprecise_te.event_consensus_len = 900;
    imprecise_te.best_te_identity = 0.77;
    imprecise_te.best_te_query_coverage = 0.98;
    imprecise_te.cross_family_margin = 0.74;
    imprecise_te.te_posterior = 0.84;
    imprecise_te.artifact_posterior = 0.12;
    imprecise_te.mechanistic_blocks = "event;independent;sequence;boundary;ref_conflict";
    result.final_calls.push_back(imprecise_te);

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().pos == 1450);
    assert(result.final_calls.front().final_qc.find("PASS_TE_IMPRECISE") == 0);
}

void finalization_keeps_precise_te_origin_over_broad_structural_envelope() {
    placer::PipelineResult result;

    placer::FinalCall precise_te;
    precise_te.chrom = "chr1";
    precise_te.tid = 0;
    precise_te.pos = 1000;
    precise_te.bp_left = 1000;
    precise_te.bp_right = 1000;
    precise_te.family = "L1";
    precise_te.subfamily = "L1";
    precise_te.te_name = "L1";
    precise_te.final_qc = "PASS_TE_IMPRECISE|EVENT_CLUSTER_PROMOTED";
    precise_te.best_explanation = "EVENT_CLUSTER";
    precise_te.support_reads = 4;
    precise_te.alt_struct_reads = 4;
    precise_te.ref_span_reads = 18;
    precise_te.event_consensus_len = 608;
    precise_te.best_te_identity = 0.943038;
    precise_te.best_te_query_coverage = 0.997845;
    precise_te.cross_family_margin = 0.419673;
    precise_te.te_posterior = 0.993225;
    precise_te.non_te_posterior = 0.00370143;
    precise_te.artifact_posterior = 0.00307329;
    precise_te.posterior_qc = "PASS_TE_POSTERIOR";
    precise_te.lfdr_qc = "PASS_TE_LFDR";
    precise_te.worst_case_lfdr = 0.00796805;
    precise_te.te_structure_path =
        "0-463:TE_CORE:L1/UNKNOWN,463-464:UNEXPLAINED:HIGH_COMPLEXITY_RESIDUAL";
    precise_te.te_structure_log_evidence = 4.07764;
    precise_te.nonte_structure_log_evidence = -2.29372;
    precise_te.artifact_structure_log_evidence = -1.75523;
    precise_te.te_structure_path_confidence = 0.997079;
    precise_te.mechanistic_blocks =
        "event;independent;sequence;structure;boundary;ref_conflict";
    result.final_calls.push_back(precise_te);

    placer::FinalCall broad_structural;
    broad_structural.chrom = "chr1";
    broad_structural.tid = 0;
    broad_structural.pos = 2186;
    broad_structural.bp_left = 2186;
    broad_structural.bp_right = 2186;
    broad_structural.family = "UNKNOWN";
    broad_structural.subfamily = "UNKNOWN";
    broad_structural.te_name = "UNKNOWN";
    broad_structural.final_qc = "PASS_STRUCTURAL_INSERTION|EVENT_CLUSTER_PROMOTED";
    broad_structural.best_explanation = "EVENT_CLUSTER";
    broad_structural.support_reads = 13;
    broad_structural.alt_struct_reads = 13;
    broad_structural.ref_span_reads = 4;
    broad_structural.event_consensus_len = 5549;
    broad_structural.best_te_identity = 0.984375;
    broad_structural.best_te_query_coverage = 0.80825;
    broad_structural.cross_family_margin = 0.792648;
    broad_structural.te_posterior = 0.996436;
    broad_structural.non_te_posterior = 0.00333243;
    broad_structural.artifact_posterior = 0.0002312;
    broad_structural.posterior_qc = "PASS_TE_POSTERIOR";
    broad_structural.lfdr_qc = "PASS_TE_LFDR";
    broad_structural.worst_case_lfdr = 0.000159606;
    broad_structural.te_structure_path =
        "0-4350:TE_CORE:Gypsy/Gypsy-34,4350-5381:TRANSDUCTION:HIGH_COMPLEXITY_TAIL";
    broad_structural.te_structure_log_evidence = 3.36813;
    broad_structural.nonte_structure_log_evidence = -1.66855;
    broad_structural.artifact_structure_log_evidence = -1.60612;
    broad_structural.te_structure_path_confidence = 0.993134;
    broad_structural.mechanistic_blocks =
        "event;independent;sequence;structure;boundary;ref_conflict";
    result.final_calls.push_back(broad_structural);

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().pos == 1000);
    assert(result.final_calls.front().final_qc.find("PASS_TE_IMPRECISE") == 0);
}

void finalization_keeps_overlapping_structural_promoted_envelopes_distinct() {
    placer::PipelineResult result;

    placer::FinalCall left;
    left.chrom = "chr1";
    left.tid = 0;
    left.pos = 1000;
    left.bp_left = 950;
    left.bp_right = 1050;
    left.family = "UNKNOWN";
    left.subfamily = "UNKNOWN";
    left.te_name = "UNKNOWN";
    left.final_qc = "PASS_STRUCTURAL_INSERTION|EVENT_CLUSTER_PROMOTED|EVENT_COMMUNITY_STABLE";
    left.support_reads = 10;
    left.alt_struct_reads = 10;
    left.ref_span_reads = 30;
    left.event_consensus_len = 900;
    left.te_posterior = 0.15;
    left.non_te_posterior = 0.70;
    left.artifact_posterior = 0.15;
    left.mechanistic_blocks = "event;community;boundary;ref_conflict";
    result.final_calls.push_back(left);

    placer::FinalCall right = left;
    right.pos = 1450;
    right.bp_left = 1400;
    right.bp_right = 1500;
    right.support_reads = 9;
    right.alt_struct_reads = 9;
    right.ref_span_reads = 25;
    right.event_consensus_len = 900;
    right.te_posterior = 0.20;
    right.non_te_posterior = 0.72;
    right.artifact_posterior = 0.08;
    result.final_calls.push_back(right);

    placer::finalize_final_calls(result);

    bool saw_left = false;
    bool saw_right = false;
    for (const auto& call : result.final_calls) {
        saw_left = saw_left || call.pos == 1000;
        saw_right = saw_right || call.pos == 1450;
    }
    assert(result.final_calls.size() == 2);
    assert(saw_left);
    assert(saw_right);
}

void finalization_collapses_shared_support_structural_fragments() {
    placer::PipelineResult result;

    auto left = make_promotable_ledger_row(1000, 12, 0.58, 0.52, 0.05, 3);
    left.bp_left = 950;
    left.bp_right = 1050;
    left.coverage_left = 950;
    left.coverage_right = 1050;
    left.family = "UNKNOWN";
    left.subfamily = "UNKNOWN";
    left.final_qc = "TE_AMBIGUOUS";
    left.posterior_qc = "TE_POSTERIOR_LOW";
    left.lfdr_qc = "TE_LFDR_HIGH";
    left.te_posterior = 0.05;
    left.non_te_posterior = 0.84;
    left.artifact_posterior = 0.11;
    left.lfdr = 0.90;
    left.worst_case_lfdr = 0.99;
    left.event_consensus_len = 900;
    left.mechanistic_lower_log_bf_te_vs_artifact = 4.5;
    left.mechanistic_lower_log_bf_te_vs_non_te = 1.5;
    left.mechanistic_ref_conflict_signal = 0.18;
    left.mechanistic_ambiguity_width = 0.10;
    left.robust_mechanistic_qc = "TE_LFDR_HIGH";
    left.robust_mechanistic_worst_case_lfdr = 0.99;
    left.mechanistic_blocks = "event;independent;sequence;boundary";
    left.support_qnames = make_qnames("shared_fragment_support_", 12);
    result.evidence_ledger.push_back(left);

    auto right = left;
    right.pos = 1450;
    right.bp_left = 1400;
    right.bp_right = 1500;
    right.coverage_left = 1400;
    right.coverage_right = 1500;
    right.alt_struct_reads = 11;
    right.ref_span_reads = 3;
    right.non_te_posterior = 0.86;
    right.artifact_posterior = 0.09;
    right.mechanistic_lower_log_bf_te_vs_artifact = 4.2;
    result.evidence_ledger.push_back(right);

    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            i % 3,
            0.50,
            0.50,
            0.01,
            5));
    }

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().final_qc.find("PASS_STRUCTURAL_INSERTION") == 0);
    assert(result.final_calls.front().final_qc.find("EVENT_FRAGMENT_COLLAPSED") !=
           std::string::npos);
}

void finalization_does_not_collapse_disjoint_fragments_with_weak_support_overlap() {
    placer::PipelineResult result;

    placer::FinalCall left;
    left.chrom = "chr1";
    left.tid = 0;
    left.pos = 1000;
    left.bp_left = 950;
    left.bp_right = 1050;
    left.family = "UNKNOWN";
    left.subfamily = "UNKNOWN";
    left.te_name = "UNKNOWN";
    left.final_qc = "PASS_STRUCTURAL_INSERTION|EVENT_CLUSTER_PROMOTED";
    left.support_reads = 10;
    left.alt_struct_reads = 10;
    left.ref_span_reads = 4;
    left.event_consensus_len = 900;
    left.non_te_posterior = 0.80;
    left.artifact_posterior = 0.10;
    left.mechanistic_blocks = "event;boundary;read_overlap";
    left.support_qnames = {
        "shared_0", "shared_1", "left_0", "left_1", "left_2",
        "left_3", "left_4", "left_5", "left_6", "left_7"};
    result.final_calls.push_back(left);

    placer::FinalCall right = left;
    right.pos = 1450;
    right.bp_left = 1400;
    right.bp_right = 1500;
    right.support_qnames = {
        "shared_0", "shared_1", "right_0", "right_1", "right_2",
        "right_3", "right_4", "right_5", "right_6", "right_7"};
    result.final_calls.push_back(right);

    placer::finalize_final_calls(result);

    bool saw_left = false;
    bool saw_right = false;
    for (const auto& call : result.final_calls) {
        saw_left = saw_left || call.pos == 1000;
        saw_right = saw_right || call.pos == 1450;
    }
    assert(result.final_calls.size() == 2);
    assert(saw_left);
    assert(saw_right);
}

void finalization_collapses_enveloped_fragments_with_shared_event_context() {
    placer::PipelineResult result;

    placer::FinalCall left;
    left.chrom = "chr1";
    left.tid = 0;
    left.pos = 1000;
    left.bp_left = 940;
    left.bp_right = 1060;
    left.window_start = 900;
    left.window_end = 1600;
    left.family = "UNKNOWN";
    left.subfamily = "UNKNOWN";
    left.te_name = "UNKNOWN";
    left.final_qc = "PASS_STRUCTURAL_INSERTION|EVENT_CLUSTER_PROMOTED";
    left.support_reads = 8;
    left.alt_struct_reads = 8;
    left.ref_span_reads = 4;
    left.max_raw_cigar_insert_len = 120;
    left.event_consensus_len = 900;
    left.non_te_posterior = 0.80;
    left.artifact_posterior = 0.10;
    left.mechanistic_blocks = "event;boundary;read_overlap";
    left.support_qnames = {
        "left_0", "left_1", "left_2", "left_3", "left_4",
        "left_5", "left_6", "left_7"};
    result.final_calls.push_back(left);

    placer::FinalCall right = left;
    right.pos = 1400;
    right.bp_left = 1340;
    right.bp_right = 1460;
    right.window_start = 920;
    right.window_end = 1620;
    right.support_qnames = {
        "right_0", "right_1", "right_2", "right_3", "right_4",
        "right_5", "right_6", "right_7"};
    result.final_calls.push_back(right);

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().final_qc.find("PASS_STRUCTURAL_INSERTION") == 0);
    assert(result.final_calls.front().final_qc.find("EVENT_FRAGMENT_COLLAPSED") !=
           std::string::npos);
    assert(result.final_calls.front().bp_left <= 940);
    assert(result.final_calls.front().bp_right >= 1460);
}

void finalization_derives_promoted_structural_context_from_event_envelope() {
    placer::PipelineResult result;

    auto left = make_promotable_ledger_row(1000, 8, 0.55, 0.95, 0.20, 2);
    left.family = "UNKNOWN";
    left.subfamily = "UNKNOWN";
    left.final_qc = "TE_AMBIGUOUS";
    left.bp_left = 940;
    left.bp_right = 1060;
    left.coverage_left = 940;
    left.coverage_right = 1060;
    left.owner_context_left = 900;
    left.owner_context_right = 1600;
    left.event_consensus_len = 900;
    left.support_qnames = make_qnames("left_support_", 8);
    result.evidence_ledger.push_back(left);

    auto right = left;
    right.pos = 1400;
    right.bp_left = 1340;
    right.bp_right = 1460;
    right.coverage_left = 1340;
    right.coverage_right = 1460;
    right.owner_context_left = 920;
    right.owner_context_right = 1620;
    right.support_qnames = make_qnames("right_support_", 8);
    result.evidence_ledger.push_back(right);

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().final_qc.find("PASS_STRUCTURAL_INSERTION") == 0);
    assert(result.final_calls.front().final_qc.find("EVENT_FRAGMENT_COLLAPSED") !=
           std::string::npos);
    assert(result.final_calls.front().bp_left <= 550);
    assert(result.final_calls.front().bp_right >= 1850);
}

void finalization_collapses_structural_candidates_from_same_owner_context() {
    placer::PipelineResult result;

    placer::FinalCall left;
    left.chrom = "chr1";
    left.tid = 0;
    left.pos = 1000;
    left.bp_left = 960;
    left.bp_right = 1040;
    left.window_start = 500;
    left.window_end = 2500;
    left.family = "UNKNOWN";
    left.subfamily = "UNKNOWN";
    left.te_name = "UNKNOWN";
    left.final_qc = "PASS_STRUCTURAL_INSERTION|EVENT_CLUSTER_PROMOTED";
    left.support_reads = 9;
    left.alt_struct_reads = 9;
    left.ref_span_reads = 4;
    left.max_raw_cigar_insert_len = 120;
    left.event_consensus_len = 120;
    left.non_te_posterior = 0.76;
    left.artifact_posterior = 0.11;
    left.mechanistic_blocks = "event;boundary;read_overlap";
    left.support_qnames = make_qnames("left_owner_context_", 9);
    result.final_calls.push_back(left);

    placer::FinalCall right = left;
    right.pos = 1900;
    right.bp_left = 1860;
    right.bp_right = 1940;
    right.support_reads = 8;
    right.alt_struct_reads = 8;
    right.support_qnames = make_qnames("right_owner_context_", 8);
    result.final_calls.push_back(right);

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().final_qc.find("PASS_STRUCTURAL_INSERTION") == 0);
    assert(result.final_calls.front().final_qc.find("EVENT_FRAGMENT_COLLAPSED") !=
           std::string::npos);
    assert(result.final_calls.front().bp_left <= 940);
    assert(result.final_calls.front().bp_right >= 1960);
}

void finalization_collapses_mixed_te_and_structural_same_owner_context() {
    placer::PipelineResult result;

    placer::FinalCall structural;
    structural.chrom = "chr1";
    structural.tid = 0;
    structural.pos = 1000;
    structural.bp_left = 940;
    structural.bp_right = 1060;
    structural.window_start = 500;
    structural.window_end = 2500;
    structural.family = "UNKNOWN";
    structural.subfamily = "UNKNOWN";
    structural.te_name = "UNKNOWN";
    structural.final_qc = "PASS_STRUCTURAL_INSERTION|EVENT_CLUSTER_PROMOTED";
    structural.support_reads = 12;
    structural.alt_struct_reads = 12;
    structural.ref_span_reads = 8;
    structural.max_raw_cigar_insert_len = 160;
    structural.event_consensus_len = 600;
    structural.non_te_posterior = 0.72;
    structural.artifact_posterior = 0.18;
    structural.mechanistic_blocks = "event;boundary;read_overlap";
    result.final_calls.push_back(structural);

    placer::FinalCall te;
    te.chrom = "chr1";
    te.tid = 0;
    te.pos = 1900;
    te.bp_left = 1880;
    te.bp_right = 1920;
    te.window_start = 520;
    te.window_end = 2520;
    te.family = "L2";
    te.subfamily = "L2";
    te.te_name = "L2";
    te.final_qc = "PASS_TE_IMPRECISE|EVENT_CLUSTER_PROMOTED";
    te.support_reads = 10;
    te.alt_struct_reads = 10;
    te.ref_span_reads = 2;
    te.max_raw_cigar_insert_len = 160;
    te.event_consensus_len = 600;
    te.best_te_identity = 0.92;
    te.best_te_query_coverage = 0.94;
    te.cross_family_margin = 0.30;
    te.te_posterior = 0.91;
    te.artifact_posterior = 0.04;
    te.posterior_qc = "PASS_TE_POSTERIOR";
    te.lfdr_qc = "PASS_TE_LFDR";
    te.worst_case_lfdr = 0.02;
    te.te_structure_path = "0-550:TE_CORE:L2/L2";
    te.te_structure_path_confidence = 0.98;
    te.mechanistic_blocks = "event;sequence;structure;boundary;read_overlap";
    result.final_calls.push_back(te);

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().final_qc.find("PASS_TE_IMPRECISE") == 0);
    assert(result.final_calls.front().final_qc.find("EVENT_FRAGMENT_COLLAPSED") !=
           std::string::npos);
    assert(result.final_calls.front().bp_left <= 700);
    assert(result.final_calls.front().bp_right >= 2200);
}

void finalization_collapses_overlapping_te_representatives_to_precise_call() {
    placer::PipelineResult result;

    auto precise = make_call(1000, 24, 0.93, 0.97, 0.20, 0);
    precise.bp_left = 1000;
    precise.bp_right = 1000;
    precise.insert_len = 240;
    precise.event_consensus_len = 408;
    precise.final_qc =
        "PASS_TE_CLOSED|PASS_INSERT_TE_ALIGNMENT_UNKNOWN|"
        "PASS_BOUNDARY_BLUNT|PASS_EVENT_CONSENSUS|PASS_EVENT_SEGMENTATION";
    precise.robust_mechanistic_qc = "PASS_TE_LFDR";
    precise.robust_mechanistic_worst_case_lfdr = 0.012;
    precise.posterior_qc = "PASS_TE_POSTERIOR";
    precise.lfdr_qc = "TE_LFDR_HIGH";
    precise.te_posterior = 0.74;
    precise.worst_case_lfdr = 0.98;
    precise.mechanistic_lower_log_bf_te_vs_artifact = 9.3;
    precise.mechanistic_lower_log_bf_te_vs_non_te = 7.7;
    precise.mechanistic_blocks = "event;independent;sequence;structure;boundary;ref_conflict";
    result.final_calls.push_back(precise);

    auto imprecise = precise;
    imprecise.pos = 1100;
    imprecise.bp_left = 1000;
    imprecise.bp_right = 1200;
    imprecise.insert_len = 0;
    imprecise.support_reads = 25;
    imprecise.alt_struct_reads = 25;
    imprecise.final_qc = "PASS_TE_IMPRECISE|EVENT_CLUSTER_PROMOTED";
    imprecise.best_explanation = "EVENT_CLUSTER";
    result.final_calls.push_back(imprecise);

    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            i % 3,
            0.50,
            0.50,
            0.01,
            5));
    }

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().pos == 1000);
    assert(result.final_calls.front().final_qc.find("PASS_TE_CLOSED") == 0);
}

void finalization_collapsed_fragment_component_preserves_union_interval() {
    placer::PipelineResult result;

    placer::FinalCall left;
    left.chrom = "chr3";
    left.tid = 2;
    left.pos = 78518954;
    left.bp_left = 78518854;
    left.bp_right = 78519054;
    left.family = "UNKNOWN";
    left.subfamily = "UNKNOWN";
    left.te_name = "UNKNOWN";
    left.final_qc = "PASS_STRUCTURAL_INSERTION|EVENT_CLUSTER_PROMOTED";
    left.alt_struct_reads = 3;
    left.support_reads = 3;
    left.raw_cigar_insert_reads = 15;
    left.max_raw_cigar_insert_len = 5;
    left.ref_span_reads = 30;
    left.event_consensus_len = 1838;
    left.mechanistic_blocks = "event;boundary;read_overlap";
    left.support_qnames = {"shared_0", "shared_1", "left_0"};
    result.final_calls.push_back(left);

    placer::FinalCall right = left;
    right.pos = 78519316;
    right.bp_left = 78519316;
    right.bp_right = 78519316;
    right.final_qc = "PASS_STRUCTURAL_INSERTION|EVENT_CLUSTER_PROMOTED";
    right.support_qnames = {"shared_0", "shared_1", "left_0"};
    result.final_calls.push_back(right);

    placer::FinalCallFilterConfig filter_config;
    filter_config.min_raw_cigar_insert_len_bp = 0;
    placer::finalize_final_calls(result, 0.10, filter_config);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().bp_left <= 78518854);
    assert(result.final_calls.front().bp_right >= 78519316);
}

void finalization_keeps_adjacent_nonoverlapping_te_events_distinct() {
    placer::PipelineResult result;

    auto left = make_call(1000, 24, 0.98, 0.96, 0.25, 0);
    left.bp_left = 990;
    left.bp_right = 1010;
    left.event_consensus_len = 120;
    left.final_qc = "PASS_TE_CLOSED|PASS_INSERT_TE_ALIGNMENT|PASS_BOUNDARY_BLUNT";
    left.robust_mechanistic_qc = "PASS_TE_LFDR";
    left.robust_mechanistic_worst_case_lfdr = 0.01;
    left.mechanistic_blocks = "event;independent;sequence;structure;boundary;ref_conflict";
    result.final_calls.push_back(left);

    auto right = left;
    right.pos = 1230;
    right.bp_left = 1220;
    right.bp_right = 1240;
    right.support_reads = 23;
    right.alt_struct_reads = 23;
    result.final_calls.push_back(right);

    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            i % 3,
            0.50,
            0.50,
            0.01,
            5));
    }

    placer::finalize_final_calls(result);

    bool saw_left = false;
    bool saw_right = false;
    for (const auto& call : result.final_calls) {
        saw_left = saw_left || call.pos == 1000;
        saw_right = saw_right || call.pos == 1230;
    }
    assert(result.final_calls.size() == 2);
    assert(saw_left);
    assert(saw_right);
}

void finalization_filters_calls_without_long_raw_cigar_insertion() {
    placer::PipelineResult result;

    auto no_long_cigar_insert = make_call(1000, 12, 0.96, 0.94, 0.20, 0);
    no_long_cigar_insert.insert_len = 800;
    no_long_cigar_insert.max_raw_cigar_insert_len = 50;
    result.final_calls.push_back(no_long_cigar_insert);

    auto long_cigar_insert = make_call(2000, 12, 0.96, 0.94, 0.20, 0);
    long_cigar_insert.insert_len = 0;
    long_cigar_insert.max_raw_cigar_insert_len = 51;
    result.final_calls.push_back(long_cigar_insert);

    placer::FinalCallFilterConfig filter_config;
    filter_config.min_raw_cigar_insert_len_bp = 50;
    placer::finalize_final_calls(result, 0.10, filter_config);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().pos == 2000);
}

void finalization_filters_imprecise_call_without_long_raw_cigar_insertion() {
    placer::PipelineResult result;

    placer::FinalCall imprecise_te;
    imprecise_te.chrom = "chr1";
    imprecise_te.tid = 0;
    imprecise_te.pos = 3000;
    imprecise_te.bp_left = 3000;
    imprecise_te.bp_right = 3000;
    imprecise_te.family = "L2";
    imprecise_te.subfamily = "L2";
    imprecise_te.te_name = "L2";
    imprecise_te.final_qc = "PASS_TE_IMPRECISE|EVENT_CLUSTER_PROMOTED";
    imprecise_te.insert_len = 0;
    imprecise_te.event_consensus_len = 900;
    imprecise_te.max_raw_cigar_insert_len = 0;
    imprecise_te.support_reads = 4;
    imprecise_te.alt_struct_reads = 4;
    imprecise_te.ref_span_reads = 6;
    result.final_calls.push_back(imprecise_te);

    placer::FinalCallFilterConfig filter_config;
    filter_config.min_raw_cigar_insert_len_bp = 50;
    placer::finalize_final_calls(result, 0.10, filter_config);

    assert(result.final_calls.empty());
}

void finalization_allows_disabling_raw_cigar_insert_filter() {
    placer::PipelineResult result;

    auto short_insert = make_call(1000, 12, 0.96, 0.94, 0.20, 0);
    short_insert.insert_len = 50;
    short_insert.max_raw_cigar_insert_len = 0;
    result.final_calls.push_back(short_insert);

    placer::FinalCallFilterConfig filter_config;
    filter_config.min_raw_cigar_insert_len_bp = 0;
    placer::finalize_final_calls(result, 0.10, filter_config);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().pos == 1000);
}

void finalization_uses_configured_raw_cigar_insert_filter() {
    placer::PipelineResult result;

    auto medium_insert = make_call(1000, 12, 0.96, 0.94, 0.20, 0);
    medium_insert.max_raw_cigar_insert_len = 80;
    result.final_calls.push_back(medium_insert);

    auto long_insert = make_call(2000, 12, 0.96, 0.94, 0.20, 0);
    long_insert.max_raw_cigar_insert_len = 120;
    result.final_calls.push_back(long_insert);

    placer::FinalCallFilterConfig filter_config;
    filter_config.min_raw_cigar_insert_len_bp = 100;
    placer::finalize_final_calls(result, 0.10, filter_config);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().pos == 2000);
}

void finalization_keeps_assembled_long_insertion_without_single_long_raw_cigar_insert() {
    placer::PipelineResult result;

    placer::FinalCall assembled;
    assembled.chrom = "chr1";
    assembled.tid = 0;
    assembled.pos = 11995657;
    assembled.bp_left = 11995657;
    assembled.bp_right = 11995657;
    assembled.family = "UNKNOWN";
    assembled.subfamily = "UNKNOWN";
    assembled.te_name = "UNKNOWN";
    assembled.final_qc = "PASS_STRUCTURAL_INSERTION";
    assembled.alt_struct_reads = 4;
    assembled.support_reads = 4;
    assembled.ref_span_reads = 33;
    assembled.raw_cigar_insert_reads = 10;
    assembled.max_raw_cigar_insert_len = 4;
    assembled.event_consensus_len = 255;
    assembled.best_te_identity = 0.572034;
    assembled.best_te_query_coverage = 0.960199;
    assembled.cross_family_margin = 0.48459;
    assembled.te_structure_path = "0-245:TE_CORE:UNKNOWN/UNKNOWN";
    assembled.te_structure_log_evidence = 3.00705;
    assembled.nonte_structure_log_evidence = -0.25;
    assembled.artifact_structure_log_evidence = 0.95;
    assembled.te_structure_path_confidence = 0.988436;
    assembled.te_posterior = 0.02;
    assembled.non_te_posterior = 0.033;
    assembled.artifact_posterior = 0.947;
    assembled.mechanistic_ref_conflict_signal = 0.40;
    assembled.mechanistic_blocks =
        "event;independent;sequence;structure;boundary;ref_conflict";
    result.final_calls.push_back(assembled);

    placer::FinalCallFilterConfig filter_config;
    filter_config.min_raw_cigar_insert_len_bp = 50;
    placer::finalize_final_calls(result, 0.10, filter_config);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().pos == 11995657);
}

void finalization_filters_assembled_artifact_fragment_without_structure_support() {
    placer::PipelineResult result;

    placer::FinalCall artifact_fragment;
    artifact_fragment.chrom = "chr1";
    artifact_fragment.tid = 0;
    artifact_fragment.pos = 78520929;
    artifact_fragment.bp_left = 78520862;
    artifact_fragment.bp_right = 78520997;
    artifact_fragment.family = "UNKNOWN";
    artifact_fragment.subfamily = "UNKNOWN";
    artifact_fragment.te_name = "UNKNOWN";
    artifact_fragment.final_qc = "PASS_STRUCTURAL_INSERTION|EVENT_CLUSTER_PROMOTED";
    artifact_fragment.alt_struct_reads = 3;
    artifact_fragment.support_reads = 3;
    artifact_fragment.ref_span_reads = 23;
    artifact_fragment.raw_cigar_insert_reads = 12;
    artifact_fragment.max_raw_cigar_insert_len = 4;
    artifact_fragment.event_consensus_len = 170;
    artifact_fragment.te_structure_path = "NA";
    artifact_fragment.te_structure_log_evidence = -3.02504;
    artifact_fragment.nonte_structure_log_evidence = -0.25;
    artifact_fragment.artifact_structure_log_evidence = -0.10;
    artifact_fragment.te_posterior = 0.0000214942;
    artifact_fragment.non_te_posterior = 0.00677296;
    artifact_fragment.artifact_posterior = 0.993206;
    artifact_fragment.mechanistic_ref_conflict_signal = 0.884615;
    artifact_fragment.mechanistic_blocks =
        "event;independent;sequence;structure;boundary;ref_conflict";
    result.final_calls.push_back(artifact_fragment);

    placer::FinalCallFilterConfig filter_config;
    filter_config.min_raw_cigar_insert_len_bp = 50;
    placer::finalize_final_calls(result, 0.10, filter_config);

    assert(result.final_calls.empty());
}

void finalization_promotes_ref_unopposed_bilateral_partial_anchor_as_unknown_structural_event() {
    placer::PipelineResult result;

    auto bilateral = make_promotable_ledger_row(6035565, 2, 0.0, 0.0, 0.0, 0);
    bilateral.bp_left = 6035531;
    bilateral.bp_right = 6035600;
    bilateral.coverage_left = 6035531;
    bilateral.coverage_right = 6035600;
    bilateral.family = "NA";
    bilateral.subfamily = "NA";
    bilateral.final_qc = "NO_CALL_INCOMPLETE";
    bilateral.posterior_qc = "TE_POSTERIOR_LOW";
    bilateral.lfdr_qc = "TE_LFDR_HIGH";
    bilateral.te_posterior = 0.000761763;
    bilateral.non_te_posterior = 0.00;
    bilateral.artifact_posterior = 0.99;
    bilateral.lfdr = 1.0;
    bilateral.worst_case_lfdr = 1.0;
    bilateral.mechanistic_lower_log_bf_te_vs_artifact = -4.03431;
    bilateral.mechanistic_lower_log_bf_te_vs_non_te = -3.35965;
    bilateral.mechanistic_ref_conflict_signal = 0.0;
    bilateral.mechanistic_ambiguity_width = 0.491361;
    bilateral.mechanistic_blocks = "event;independent;boundary;ref_conflict";
    bilateral.robust_mechanistic_qc = "TE_LFDR_HIGH";
    bilateral.robust_mechanistic_worst_case_lfdr = 1.0;
    bilateral.alt_split_reads = 0;
    bilateral.alt_indel_reads = 0;
    bilateral.alt_left_clip_reads = 1;
    bilateral.alt_right_clip_reads = 1;
    bilateral.partial_context_input_reads = 2;
    bilateral.input_event_reads = 2;
    bilateral.left_anchor_input_reads = 1;
    bilateral.right_anchor_input_reads = 1;
    bilateral.event_consensus_len = 185;
    bilateral.support_qnames = {"partial_left_anchor", "partial_right_anchor"};
    result.evidence_ledger.push_back(bilateral);

    auto one_sided = bilateral;
    one_sided.pos = 6037000;
    one_sided.bp_left = 6036970;
    one_sided.bp_right = 6037030;
    one_sided.coverage_left = 6036970;
    one_sided.coverage_right = 6037030;
    one_sided.alt_right_clip_reads = 0;
    one_sided.right_anchor_input_reads = 0;
    one_sided.support_qnames = {"one_sided_0", "one_sided_1"};
    result.evidence_ledger.push_back(one_sided);

    auto ref_opposed = bilateral;
    ref_opposed.pos = 6039000;
    ref_opposed.bp_left = 6038970;
    ref_opposed.bp_right = 6039030;
    ref_opposed.coverage_left = 6038970;
    ref_opposed.coverage_right = 6039030;
    ref_opposed.ref_span_reads = 1;
    ref_opposed.support_qnames = {"ref_opposed_0", "ref_opposed_1"};
    result.evidence_ledger.push_back(ref_opposed);

    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            4,
            0.50,
            0.50,
            0.01,
            5));
    }

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    const auto& call = result.final_calls.front();
    assert(call.pos == 6035565);
    assert(call.bp_left <= 6035365);
    assert(call.bp_right >= 6035361);
    assert(call.family == "UNKNOWN");
    assert(call.subfamily == "UNKNOWN");
    assert(call.te_name == "UNKNOWN");
    assert(call.artifact_posterior == 0.99);
    assert(call.final_qc.find("PASS_STRUCTURAL_INSERTION") == 0);
    assert(call.final_qc.find("EVENT_BILATERAL_PARTIAL_ANCHOR") !=
           std::string::npos);
    assert(call.conformal_qc == "PASS_EVENT_EBH");
}

void finalization_rejects_low_allele_fraction_event_without_error_model_support() {
    placer::PipelineResult result;
    auto low_support = make_promotable_ledger_row(1000, 2, 0.58, 0.97, 0.09, 30);
    low_support.family = "UNKNOWN";
    low_support.subfamily = "UNKNOWN";
    low_support.final_qc = "TE_AMBIGUOUS";
    low_support.posterior_qc = "TE_POSTERIOR_LOW";
    low_support.lfdr_qc = "TE_LFDR_HIGH";
    low_support.te_posterior = 0.12;
    low_support.lfdr = 0.90;
    low_support.worst_case_lfdr = 0.99;
    low_support.mechanistic_lower_log_bf_te_vs_artifact = -5.70;
    low_support.mechanistic_lower_log_bf_te_vs_non_te = -1.80;
    low_support.mechanistic_ref_conflict_signal = 0.94;
    low_support.mechanistic_ambiguity_width = 0.48;
    low_support.robust_mechanistic_qc = "TE_LFDR_HIGH";
    low_support.robust_mechanistic_worst_case_lfdr = 0.999;
    low_support.mechanistic_blocks = "event;independent;sequence;boundary;ref_conflict";
    result.evidence_ledger.push_back(low_support);

    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            i % 3,
            0.50,
            0.50,
            0.01,
            8));
    }

    placer::finalize_final_calls(result);

    assert(result.final_calls.empty());
}

void finalization_keeps_balanced_heterozygous_structural_insertion() {
    placer::PipelineResult result;
    placer::FinalCall call;
    call.chrom = "chr1";
    call.tid = 0;
    call.pos = 1000;
    call.bp_left = 1000;
    call.bp_right = 1000;
    call.family = "UNKNOWN";
    call.subfamily = "UNKNOWN";
    call.te_name = "UNKNOWN";
    call.final_qc = "PASS_STRUCTURAL_INSERTION";
    call.alt_struct_reads = 20;
    call.support_reads = 20;
    call.ref_span_reads = 16;
    call.best_te_identity = 0.0;
    call.best_te_query_coverage = 0.0;
    call.cross_family_margin = 0.0;
    call.non_te_posterior = 0.72;
    call.artifact_posterior = 0.25;
    call.posterior_qc = "TE_POSTERIOR_LOW";
    call.lfdr_qc = "TE_LFDR_HIGH";
    call.mechanistic_lower_log_bf_te_vs_artifact = -0.91;
    call.mechanistic_lower_log_bf_te_vs_non_te = -0.54;
    call.mechanistic_ref_conflict_signal = 0.52;
    call.mechanistic_ambiguity_width = 0.38;
    call.mechanistic_blocks = "event;independent;boundary";
    result.final_calls.push_back(call);

    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            45,
            0.70,
            0.40,
            0.03,
            0));
    }

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().pos == 1000);
    assert(result.final_calls.front().conformal_qc == "PASS_EVENT_EXISTENCE_BFDR");
}

void finalization_rejects_ref_dominated_structural_candidate() {
    placer::PipelineResult result;
    placer::FinalCall call;
    call.chrom = "chr1";
    call.tid = 0;
    call.pos = 1000;
    call.bp_left = 1000;
    call.bp_right = 1000;
    call.family = "UNKNOWN";
    call.subfamily = "UNKNOWN";
    call.te_name = "UNKNOWN";
    call.final_qc = "PASS_STRUCTURAL_INSERTION";
    call.alt_struct_reads = 2;
    call.support_reads = 2;
    call.ref_span_reads = 30;
    call.best_te_identity = 0.58;
    call.best_te_query_coverage = 0.97;
    call.cross_family_margin = 0.09;
    call.non_te_posterior = 0.20;
    call.artifact_posterior = 0.74;
    call.posterior_qc = "TE_POSTERIOR_LOW";
    call.lfdr_qc = "TE_LFDR_HIGH";
    call.mechanistic_lower_log_bf_te_vs_artifact = -5.7;
    call.mechanistic_lower_log_bf_te_vs_non_te = -1.8;
    call.mechanistic_ref_conflict_signal = 0.94;
    call.mechanistic_ambiguity_width = 0.48;
    call.mechanistic_blocks = "event;independent;boundary";
    result.final_calls.push_back(call);

    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            4,
            0.50,
            0.50,
            0.01,
            5));
    }

    placer::finalize_final_calls(result);

    assert(result.final_calls.empty());
}

void finalization_promotes_competing_nonte_long_insertion_event() {
    placer::PipelineResult result;

    auto nonte = make_promotable_ledger_row(1000, 18, 0.58, 0.56, 0.35, 3);
    nonte.family = "UNKNOWN";
    nonte.subfamily = "UNKNOWN";
    nonte.final_qc = "TE_AMBIGUOUS";
    nonte.posterior_qc = "TE_POSTERIOR_LOW";
    nonte.lfdr_qc = "TE_LFDR_HIGH";
    nonte.te_posterior = 0.04;
    nonte.non_te_posterior = 0.74;
    nonte.artifact_posterior = 0.22;
    nonte.lfdr = 0.90;
    nonte.worst_case_lfdr = 0.99;
    nonte.robust_mechanistic_qc = "TE_LFDR_HIGH";
    nonte.robust_mechanistic_worst_case_lfdr = 0.98;
    nonte.mechanistic_lower_log_bf_te_vs_artifact = -0.30;
    nonte.mechanistic_lower_log_bf_te_vs_non_te = -0.10;
    nonte.mechanistic_ref_conflict_signal = 0.18;
    nonte.mechanistic_ambiguity_width = 0.36;
    nonte.raw_cigar_insert_reads = 36;
    nonte.max_raw_cigar_insert_len = 740;
    nonte.event_consensus_len = 1200;
    nonte.te_core_coverage = 0.56;
    nonte.unexplained_high_complexity_bp = 520;
    nonte.mechanistic_blocks =
        "event;independent;sequence;boundary;ref_conflict";
    result.evidence_ledger.push_back(nonte);

    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            i % 3,
            0.50,
            0.50,
            0.01,
            5));
    }

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().pos == 1000);
    assert(result.final_calls.front().family == "UNKNOWN");
    assert(result.final_calls.front().subfamily == "UNKNOWN");
    assert(result.final_calls.front().final_qc.find("PASS_STRUCTURAL_INSERTION") == 0);
    assert(result.final_calls.front().conformal_qc == "PASS_CONFORMAL_FDR" ||
           result.final_calls.front().conformal_qc == "PASS_EVENT_EXISTENCE_BFDR" ||
           result.final_calls.front().conformal_qc == "PASS_EVENT_EBH");
}

void finalization_promotes_robust_structural_insertion_when_conformal_nulls_are_stronger() {
    placer::PipelineResult result;

    auto row = make_promotable_ledger_row(1000, 40, 0.544992, 0.997984, 0.511635, 0);
    row.family = "Gypsy";
    row.subfamily = "Gypsy-26";
    row.final_qc = "PASS_TE_CLOSED";
    row.posterior_qc = "TE_POSTERIOR_LOW";
    row.lfdr_qc = "TE_LFDR_HIGH";
    row.te_posterior = 0.404234;
    row.non_te_posterior = 0.587561;
    row.artifact_posterior = 0.00820492;
    row.lfdr = 0.822969;
    row.worst_case_lfdr = 0.990197;
    row.robust_mechanistic_qc = "PASS_TE_LFDR";
    row.robust_mechanistic_worst_case_lfdr = 0.0415926;
    row.raw_cigar_insert_reads = 27;
    row.max_raw_cigar_insert_len = 469;
    row.event_consensus_len = 626;
    row.te_structure_path =
        "0-495:TE_CORE:UNKNOWN/UNKNOWN,495-496:UNEXPLAINED:HIGH_COMPLEXITY_RESIDUAL";
    row.te_structure_log_evidence = 2.47112;
    row.nonte_structure_log_evidence = -1.95588;
    row.artifact_structure_log_evidence = -1.45683;
    row.te_structure_path_confidence = 0.980696;
    row.te_core_coverage = 0.997984;
    row.unexplained_high_complexity_bp = 1;
    row.ref_span_reads = 25;
    row.mechanistic_lower_log_bf_te_vs_artifact = 0.20;
    row.mechanistic_lower_log_bf_te_vs_non_te = 0.10;
    row.mechanistic_ref_conflict_signal = 0.0;
    row.mechanistic_ambiguity_width = 0.294207;
    row.mechanistic_blocks =
        "event;independent;sequence;structure;boundary;ref_conflict";
    result.evidence_ledger.push_back(row);

    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            41,
            0.60,
            0.99,
            0.52,
            25));
    }

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().pos == 1000);
    assert(result.final_calls.front().family == "Gypsy");
    assert(result.final_calls.front().final_qc.find("PASS_TE_IMPRECISE") == 0);
    assert(result.final_calls.front().conformal_qc == "PASS_CONFORMAL_FDR" ||
           result.final_calls.front().conformal_qc == "PASS_EVENT_BFDR" ||
           result.final_calls.front().conformal_qc == "PASS_EVENT_EXISTENCE_BFDR");
}

void finalization_combines_te_structure_origin_with_posterior_evidence() {
    placer::PipelineResult result;

    auto row = make_promotable_ledger_row(1000, 1, 0.698198, 0.990654, 0.667374, 6);
    row.family = "ERV1";
    row.subfamily = "NA";
    row.final_qc = "TE_AMBIGUOUS";
    row.posterior_qc = "PASS_TE_POSTERIOR";
    row.lfdr_qc = "TE_LFDR_HIGH";
    row.te_posterior = 0.892946;
    row.non_te_posterior = 0.0571762;
    row.artifact_posterior = 0.0498781;
    row.lfdr = 0.940697;
    row.worst_case_lfdr = 0.999674;
    row.robust_mechanistic_qc = "TE_LFDR_HIGH";
    row.robust_mechanistic_worst_case_lfdr = 0.952534;
    row.raw_cigar_insert_reads = 3;
    row.max_raw_cigar_insert_len = 539;
    row.event_consensus_len = 695;
    row.te_structure_path =
        "0-530:TE_CORE:ERV1/UNKNOWN,530-535:UNEXPLAINED:HIGH_COMPLEXITY_RESIDUAL";
    row.te_structure_log_evidence = 2.9971;
    row.nonte_structure_log_evidence = -2.11043;
    row.artifact_structure_log_evidence = -1.44141;
    row.te_structure_path_confidence = 0.988324;
    row.te_core_coverage = 0.990654;
    row.mechanistic_lower_log_bf_te_vs_artifact = 0.612855;
    row.mechanistic_lower_log_bf_te_vs_non_te = 3.90561;
    row.mechanistic_ref_conflict_signal = 0.857143;
    row.mechanistic_ambiguity_width = 0.393099;
    row.mechanistic_blocks =
        "event;independent;sequence;structure;boundary;ref_conflict";
    result.evidence_ledger.push_back(row);

    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            i % 3,
            0.50,
            0.50,
            0.01,
            5));
    }

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().pos == 1000);
    assert(result.final_calls.front().final_qc.find("PASS_TE_IMPRECISE") == 0);
    assert(result.final_calls.front().conformal_qc == "PASS_EVENT_BFDR" ||
           result.final_calls.front().conformal_qc == "PASS_CONFORMAL_FDR");
}

void finalization_promotes_directly_observed_nonte_raw_insertion_event() {
    placer::PipelineResult result;

    auto raw_observed = make_promotable_ledger_row(1000, 1, 0.0, 0.0, 0.0, 22);
    raw_observed.family = "UNKNOWN";
    raw_observed.subfamily = "UNKNOWN";
    raw_observed.final_qc = "REFERENCE_OR_ARTIFACT";
    raw_observed.posterior_qc = "TE_POSTERIOR_LOW";
    raw_observed.lfdr_qc = "TE_LFDR_HIGH";
    raw_observed.te_posterior = 0.004865;
    raw_observed.non_te_posterior = 0.833204;
    raw_observed.artifact_posterior = 0.161931;
    raw_observed.lfdr = 1.0;
    raw_observed.worst_case_lfdr = 1.0;
    raw_observed.robust_mechanistic_qc = "TE_LFDR_HIGH";
    raw_observed.robust_mechanistic_worst_case_lfdr = 1.0;
    raw_observed.raw_cigar_insert_reads = 5;
    raw_observed.max_raw_cigar_insert_len = 118;
    raw_observed.event_consensus_len = 278;
    raw_observed.te_structure_path = "NA";
    raw_observed.te_structure_log_evidence = -3.71014;
    raw_observed.nonte_structure_log_evidence = -0.281865;
    raw_observed.artifact_structure_log_evidence = 1.15;
    raw_observed.te_structure_path_confidence = 0.00768978;
    raw_observed.mechanistic_lower_log_bf_te_vs_artifact = -11.7551;
    raw_observed.mechanistic_lower_log_bf_te_vs_non_te = -7.86282;
    raw_observed.mechanistic_ref_conflict_signal = 0.956522;
    raw_observed.mechanistic_ambiguity_width = 0.533527;
    raw_observed.mechanistic_blocks =
        "event;independent;sequence;structure;boundary;ref_conflict";
    result.evidence_ledger.push_back(raw_observed);

    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            i % 3,
            0.50,
            0.50,
            0.01,
            5));
    }

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().pos == 1000);
    assert(result.final_calls.front().family == "UNKNOWN");
    assert(result.final_calls.front().final_qc.find("PASS_STRUCTURAL_INSERTION") == 0);
    assert(result.final_calls.front().conformal_qc == "PASS_EVENT_EXISTENCE_BFDR" ||
           result.final_calls.front().conformal_qc == "PASS_CONFORMAL_FDR");
}

void finalization_keeps_s16_like_low_af_closed_structural_insertion() {
    placer::PipelineResult result;

    placer::FinalCall call;
    call.chrom = "chr1";
    call.tid = 0;
    call.pos = 43481588;
    call.bp_left = 43481588;
    call.bp_right = 43481588;
    call.family = "UNKNOWN";
    call.subfamily = "UNKNOWN";
    call.te_name = "UNKNOWN";
    call.final_qc = "PASS_STRUCTURAL_INSERTION";
    call.alt_struct_reads = 4;
    call.support_reads = 4;
    call.ref_span_reads = 31;
    call.raw_cigar_insert_reads = 3;
    call.max_raw_cigar_insert_len = 129;
    call.event_consensus_len = 289;
    call.best_te_identity = 0.853448;
    call.best_te_query_coverage = 0.822695;
    call.cross_family_margin = 0.253448;
    call.te_structure_path = "0-141:TE_CORE:ERV1/UNKNOWN";
    call.te_structure_log_evidence = 2.17092;
    call.nonte_structure_log_evidence = -2.22543;
    call.artifact_structure_log_evidence = 0.425871;
    call.te_structure_path_confidence = 0.851327;
    call.te_posterior = 0.372479;
    call.non_te_posterior = 0.217595;
    call.artifact_posterior = 0.409926;
    call.posterior_qc = "TE_POSTERIOR_LOW";
    call.lfdr_qc = "TE_LFDR_HIGH";
    call.worst_case_lfdr = 0.999359;
    call.robust_mechanistic_qc = "TE_LFDR_HIGH";
    call.robust_mechanistic_worst_case_lfdr = 0.999359;
    call.mechanistic_lower_log_bf_te_vs_artifact = -3.70692;
    call.mechanistic_lower_log_bf_te_vs_non_te = -0.26883;
    call.mechanistic_ref_conflict_signal = 0.885714;
    call.mechanistic_ambiguity_width = 0.426672;
    call.mechanistic_blocks =
        "event;independent;sequence;structure;boundary;ref_conflict";
    result.final_calls.push_back(call);

    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            4,
            0.90,
            0.90,
            0.30,
            5));
    }

    placer::finalize_final_calls(result);

    assert(result.final_calls.size() == 1);
    assert(result.final_calls.front().pos == 43481588);
    assert(result.final_calls.front().final_qc == "PASS_STRUCTURAL_INSERTION");
    assert(result.final_calls.front().conformal_qc == "PASS_EVENT_EXISTENCE_BFDR");
}

void finalization_rejects_nonte_event_when_artifact_competitor_wins() {
    placer::PipelineResult result;

    auto artifact_like = make_promotable_ledger_row(1000, 7, 0.58, 0.96, 0.45, 31);
    artifact_like.family = "UNKNOWN";
    artifact_like.subfamily = "UNKNOWN";
    artifact_like.final_qc = "TE_AMBIGUOUS";
    artifact_like.posterior_qc = "TE_POSTERIOR_LOW";
    artifact_like.lfdr_qc = "TE_LFDR_HIGH";
    artifact_like.te_posterior = 0.05;
    artifact_like.non_te_posterior = 0.20;
    artifact_like.artifact_posterior = 0.75;
    artifact_like.lfdr = 0.95;
    artifact_like.worst_case_lfdr = 1.0;
    artifact_like.robust_mechanistic_qc = "TE_LFDR_HIGH";
    artifact_like.robust_mechanistic_worst_case_lfdr = 1.0;
    artifact_like.mechanistic_lower_log_bf_te_vs_artifact = -4.8;
    artifact_like.mechanistic_lower_log_bf_te_vs_non_te = -1.5;
    artifact_like.mechanistic_ref_conflict_signal = 0.84;
    artifact_like.mechanistic_ambiguity_width = 0.45;
    artifact_like.raw_cigar_insert_reads = 15;
    artifact_like.max_raw_cigar_insert_len = 240;
    artifact_like.mechanistic_blocks =
        "event;independent;sequence;boundary;ref_conflict";
    result.evidence_ledger.push_back(artifact_like);

    for (int i = 0; i < 40; ++i) {
        result.evidence_ledger.push_back(make_null_control(
            7000 + i,
            i % 3,
            0.50,
            0.50,
            0.01,
            5));
    }

    placer::finalize_final_calls(result);

    assert(result.final_calls.empty());
}

}  // namespace

int main() {
    finalization_applies_sample_local_conformal_selector();
    finalization_preserves_legacy_dedup_when_no_ledger_is_available();
    finalization_commits_sequence_certified_family_after_selection();
    finalization_keeps_pareto_weaker_distant_final_candidate();
    finalization_does_not_remove_pareto_weaker_distant_event();
    finalization_promotes_event_cluster_from_strong_ledger_evidence();
    finalization_recovers_pass_te_ledger_row_missing_from_final_calls();
    finalization_keeps_distinct_promoted_events_in_same_coverage_bin();
    finalization_does_not_promote_reference_conflict_without_statistical_support();
    finalization_prefers_stronger_promoted_cluster_over_weak_final_call();
    finalization_uses_event_lfdr_when_conformal_null_is_too_small();
    finalization_uses_event_ebh_when_null_is_too_small_but_bf_is_strong();
    finalization_promotes_ambiguous_ledger_row_for_event_ebh_certificate();
    finalization_keeps_spatially_distinct_ebh_event_when_bfdr_event_exists();
    finalization_uses_event_bayesian_fdr_for_calibrated_lfdr_with_reference_span();
    finalization_uses_event_bayesian_fdr_for_pass_posterior_qc();
    finalization_promotes_calibrated_near_complete_te_structure_despite_robust_uncertainty();
    finalization_keeps_structural_insertion_with_event_existence_certificate();
    te_calibrated_report_mode_removes_structural_event_existence_call();
    te_calibrated_report_mode_keeps_calibrated_closed_te_call();
    te_calibrated_report_mode_keeps_calibrated_promoted_te_call();
    te_calibrated_report_mode_removes_te_call_without_lfdr_certificate();
    finalization_promotes_structural_ledger_event_without_te_calibration();
    finalization_promotes_balanced_ambiguous_ledger_event();
    finalization_promotes_low_allele_fraction_event_with_error_model_certificate();
    finalization_rejects_isolated_low_allele_fraction_event_without_stable_community();
    finalization_promotes_low_af_event_with_self_stable_interval();
    finalization_uses_low_af_evidence_not_fixed_read_floor_for_stable_interval();
    finalization_does_not_collapse_nonoverlapping_stable_interval_into_adjacent_te_call();
    finalization_promotes_low_af_event_with_subthreshold_stable_neighbor();
    finalization_prefers_stable_low_af_community_over_uncertified_direct_call();
    finalization_aggregates_connected_weak_event_community();
    finalization_competes_overlapping_promoted_event_envelopes();
    finalization_keeps_precise_te_origin_over_broad_structural_envelope();
    finalization_keeps_overlapping_structural_promoted_envelopes_distinct();
    finalization_collapses_shared_support_structural_fragments();
    finalization_does_not_collapse_disjoint_fragments_with_weak_support_overlap();
    finalization_collapses_enveloped_fragments_with_shared_event_context();
    finalization_derives_promoted_structural_context_from_event_envelope();
    finalization_collapses_structural_candidates_from_same_owner_context();
    finalization_collapses_mixed_te_and_structural_same_owner_context();
    finalization_collapses_overlapping_te_representatives_to_precise_call();
    finalization_collapsed_fragment_component_preserves_union_interval();
    finalization_keeps_adjacent_nonoverlapping_te_events_distinct();
    finalization_filters_calls_without_long_raw_cigar_insertion();
    finalization_filters_imprecise_call_without_long_raw_cigar_insertion();
    finalization_allows_disabling_raw_cigar_insert_filter();
    finalization_uses_configured_raw_cigar_insert_filter();
    finalization_keeps_assembled_long_insertion_without_single_long_raw_cigar_insert();
    finalization_filters_assembled_artifact_fragment_without_structure_support();
    finalization_promotes_ref_unopposed_bilateral_partial_anchor_as_unknown_structural_event();
    finalization_rejects_low_allele_fraction_event_without_error_model_support();
    finalization_keeps_balanced_heterozygous_structural_insertion();
    finalization_rejects_ref_dominated_structural_candidate();
    finalization_promotes_competing_nonte_long_insertion_event();
    finalization_promotes_robust_structural_insertion_when_conformal_nulls_are_stronger();
    finalization_combines_te_structure_origin_with_posterior_evidence();
    finalization_promotes_directly_observed_nonte_raw_insertion_event();
    finalization_keeps_s16_like_low_af_closed_structural_insertion();
    finalization_rejects_nonte_event_when_artifact_competitor_wins();
    return 0;
}
