namespace {

constexpr double kDefaultFinalConformalFdr = 0.10;

bool starts_with_text(const std::string& value, const char* prefix) {
    const std::string_view text(value);
    const std::string_view needle(prefix);
    return text.size() >= needle.size() &&
           text.substr(0, needle.size()) == needle;
}

bool is_final_te_qc_for_conformal(const std::string& final_qc) {
    return starts_with_text(final_qc, "PASS_TE");
}

bool is_final_event_qc_for_conformal(const std::string& final_qc) {
    return is_final_te_qc_for_conformal(final_qc) ||
           starts_with_text(final_qc, "PASS_STRUCTURAL_INSERTION");
}

bool final_call_has_reportable_conformal_certificate(const FinalCall& call) {
    return call.conformal_qc == "PASS_CONFORMAL_FDR" ||
           call.conformal_qc == "PASS_EVENT_LFDR_FALLBACK" ||
           call.conformal_qc == "PASS_EVENT_EBH" ||
           call.conformal_qc == "PASS_EVENT_BFDR" ||
           call.conformal_qc == "PASS_EVENT_EXISTENCE_BFDR";
}

bool final_call_passes_te_calibrated_report_gate(const FinalCall& call) {
    return starts_with_text(call.final_qc, "PASS_TE") &&
           call.posterior_qc == "PASS_TE_POSTERIOR" &&
           call.lfdr_qc == "PASS_TE_LFDR" &&
           final_call_has_reportable_conformal_certificate(call);
}

bool has_conformal_feature_certificate(const EvidenceLedgerRow& row) {
    return !row.mechanistic_blocks.empty() && row.mechanistic_blocks != "NA";
}

bool is_promoted_event_cluster_call(const FinalCall& call) {
    return call.final_qc.find("EVENT_CLUSTER_PROMOTED") != std::string::npos;
}

bool final_call_passes_min_raw_cigar_insert_length(
    const FinalCall& call,
    const FinalCallFilterConfig& filter_config) {
    const int32_t threshold = filter_config.min_raw_cigar_insert_len_bp;
    if (threshold <= 0) {
        return true;
    }
    if (call.max_raw_cigar_insert_len > threshold) {
        return true;
    }
    if (call.alt_struct_reads <= 0 || call.event_consensus_len <= threshold) {
        return false;
    }
    if (call.mechanistic_blocks.empty() || call.mechanistic_blocks == "NA") {
        return false;
    }
    const bool has_long_insertion_structure =
        call.mechanistic_blocks.find("sequence") != std::string::npos ||
        call.mechanistic_blocks.find("structure") != std::string::npos;
    if (!has_long_insertion_structure) {
        return false;
    }
    const double event_structure_log_bf =
        std::max(call.te_structure_log_evidence, call.nonte_structure_log_evidence) -
        call.artifact_structure_log_evidence;
    const double event_posterior =
        std::clamp(call.te_posterior, 0.0, 1.0) +
        std::clamp(call.non_te_posterior, 0.0, 1.0);
    const double artifact_posterior =
        std::clamp(call.artifact_posterior, 0.0, 1.0);
    return event_structure_log_bf > 0.0 || event_posterior > artifact_posterior;
}

void filter_short_final_insertions(
    PipelineResult& result,
    const FinalCallFilterConfig& filter_config) {
    result.final_calls.erase(
        std::remove_if(
            result.final_calls.begin(),
            result.final_calls.end(),
            [&filter_config](const FinalCall& call) {
                return !final_call_passes_min_raw_cigar_insert_length(
                    call,
                    filter_config);
            }),
        result.final_calls.end());
}

std::pair<int32_t, int32_t> final_call_event_envelope(const FinalCall& call) {
    int32_t left = std::min(call.bp_left, call.bp_right);
    int32_t right = std::max(call.bp_left, call.bp_right);
    if (left < 0 || right < 0) {
        left = call.pos;
        right = call.pos;
    }
    if (call.event_consensus_len > 0 && call.pos >= 0) {
        const int32_t uncertainty = call.event_consensus_len / 2;
        left = std::min(left, std::max(0, call.pos - uncertainty));
        right = std::max(right, call.pos + uncertainty);
    }
    if (left < 0 || right < 0) {
        return {-1, -1};
    }
    if (left > right) {
        std::swap(left, right);
    }
    return {left, right};
}

bool final_call_event_envelopes_overlap(const FinalCall& a, const FinalCall& b) {
    if (a.tid != b.tid || a.chrom != b.chrom) {
        return false;
    }
    const auto [a_left, a_right] = final_call_event_envelope(a);
    const auto [b_left, b_right] = final_call_event_envelope(b);
    if (a_left < 0 || a_right < 0 || b_left < 0 || b_right < 0) {
        return false;
    }
    return a_right >= b_left && b_right >= a_left;
}

bool final_call_has_explanation_posteriors(const FinalCall& call) {
    const double total = call.te_posterior + call.non_te_posterior + call.artifact_posterior;
    return std::isfinite(total) && total > 0.50;
}

double final_call_non_artifact_posterior(const FinalCall& call) {
    return std::clamp(call.te_posterior, 0.0, 1.0) +
           std::clamp(call.non_te_posterior, 0.0, 1.0);
}

bool final_call_is_te_specific_event(const FinalCall& call) {
    const bool legacy_concrete_label =
        call.family != "UNKNOWN" && call.family != "NA" && !call.family.empty();
    return starts_with_text(call.final_qc, "PASS_TE") &&
           (call.family_committed || legacy_concrete_label);
}

double combine_independent_local_fdr(double lhs, double rhs);

double final_call_statistical_te_local_fdr(const FinalCall& call) {
    double local_fdr = 1.0;
    if (call.lfdr_qc == "PASS_TE_LFDR" || call.worst_case_lfdr <= kDefaultFinalConformalFdr) {
        local_fdr = std::min(local_fdr, std::clamp(call.worst_case_lfdr, 0.0, 1.0));
    }
    if (call.posterior_qc == "PASS_TE_POSTERIOR" || call.te_posterior >= 1.0 - kDefaultFinalConformalFdr) {
        local_fdr = std::min(
            local_fdr,
            1.0 - std::clamp(call.te_posterior, 0.0, 1.0));
    }
    return local_fdr;
}

double final_call_structure_te_local_fdr(const FinalCall& call) {
    if (call.te_structure_path.empty() || call.te_structure_path == "NA") {
        return 1.0;
    }
    return 1.0 - std::clamp(call.te_structure_path_confidence, 0.0, 1.0);
}

bool final_call_has_te_origin_certificate(const FinalCall& call) {
    if (!final_call_is_te_specific_event(call)) {
        return false;
    }
    const double statistical_local_fdr = final_call_statistical_te_local_fdr(call);
    const double structure_local_fdr = final_call_structure_te_local_fdr(call);
    return combine_independent_local_fdr(statistical_local_fdr, structure_local_fdr) <=
           kDefaultFinalConformalFdr;
}

bool promoted_event_envelopes_should_compete(const FinalCall& a, const FinalCall& b) {
    if (!is_promoted_event_cluster_call(a) || !is_promoted_event_cluster_call(b)) {
        return false;
    }
    const bool a_te_specific = final_call_is_te_specific_event(a);
    const bool b_te_specific = final_call_is_te_specific_event(b);
    if (a_te_specific == b_te_specific) {
        return false;
    }
    return final_call_event_envelopes_overlap(a, b);
}

bool prefer_promoted_event_envelope_representative(
    const FinalCall& candidate,
    const FinalCall& incumbent) {
    const bool candidate_te_origin = final_call_has_te_origin_certificate(candidate);
    const bool incumbent_te_origin = final_call_has_te_origin_certificate(incumbent);
    if (candidate_te_origin != incumbent_te_origin) {
        return candidate_te_origin;
    }
    const bool candidate_has_posteriors = final_call_has_explanation_posteriors(candidate);
    const bool incumbent_has_posteriors = final_call_has_explanation_posteriors(incumbent);
    if (candidate_has_posteriors && incumbent_has_posteriors) {
        const double candidate_non_artifact = final_call_non_artifact_posterior(candidate);
        const double incumbent_non_artifact = final_call_non_artifact_posterior(incumbent);
        if (candidate_non_artifact != incumbent_non_artifact) {
            return candidate_non_artifact > incumbent_non_artifact;
        }
        if (candidate.artifact_posterior != incumbent.artifact_posterior) {
            return candidate.artifact_posterior < incumbent.artifact_posterior;
        }
    } else if (candidate_has_posteriors != incumbent_has_posteriors) {
        return candidate_has_posteriors;
    }

    const bool candidate_te_specific = final_call_is_te_specific_event(candidate);
    const bool incumbent_te_specific = final_call_is_te_specific_event(incumbent);
    if (candidate_te_specific != incumbent_te_specific) {
        return candidate_te_specific;
    }
    return prefer_new_call(candidate, incumbent);
}

double log_sum_exp_pair(double lhs, double rhs) {
    if (!std::isfinite(lhs)) {
        return rhs;
    }
    if (!std::isfinite(rhs)) {
        return lhs;
    }
    const double max_value = std::max(lhs, rhs);
    return max_value + std::log(std::exp(lhs - max_value) + std::exp(rhs - max_value));
}

double odds_to_local_fdr(double log_odds) {
    if (!std::isfinite(log_odds)) {
        return 1.0;
    }
    if (log_odds >= 60.0) {
        return 0.0;
    }
    if (log_odds <= -60.0) {
        return 1.0;
    }
    return 1.0 / (1.0 + std::exp(log_odds));
}

double local_fdr_to_log_odds(double local_fdr) {
    const double fdr = std::clamp(local_fdr, 0.0, 1.0);
    if (fdr <= 0.0) {
        return 60.0;
    }
    if (fdr >= 1.0) {
        return -std::numeric_limits<double>::infinity();
    }
    return std::log1p(-fdr) - std::log(fdr);
}

double combine_independent_local_fdr(double lhs, double rhs) {
    const double lhs_log_odds = local_fdr_to_log_odds(lhs);
    const double rhs_log_odds = local_fdr_to_log_odds(rhs);
    if (!std::isfinite(lhs_log_odds)) {
        return rhs;
    }
    if (!std::isfinite(rhs_log_odds)) {
        return lhs;
    }
    return odds_to_local_fdr(lhs_log_odds + rhs_log_odds);
}

double ledger_row_statistical_te_local_fdr(const EvidenceLedgerRow& row) {
    double local_fdr = 1.0;
    if (row.lfdr_qc == "PASS_TE_LFDR" || row.worst_case_lfdr <= 0.10) {
        local_fdr = std::min(local_fdr, std::clamp(row.worst_case_lfdr, 0.0, 1.0));
    }
    if (row.posterior_qc == "PASS_TE_POSTERIOR" || row.te_posterior >= 0.95) {
        local_fdr = std::min(
            local_fdr,
            1.0 - std::clamp(row.te_posterior, 0.0, 1.0));
    }
    return local_fdr;
}

double ledger_row_structure_te_local_fdr(const EvidenceLedgerRow& row) {
    if (row.te_structure_path.empty() || row.te_structure_path == "NA") {
        return 1.0;
    }
    const double competing_structure = log_sum_exp_pair(
        row.nonte_structure_log_evidence,
        row.artifact_structure_log_evidence);
    return odds_to_local_fdr(row.te_structure_log_evidence - competing_structure);
}

bool has_te_structure_origin_certificate(const EvidenceLedgerRow& row) {
    if (!has_conformal_feature_certificate(row)) {
        return false;
    }
    const bool direct_insert_observation =
        row.raw_cigar_insert_reads > 0 ||
        row.max_raw_cigar_insert_len > 0 ||
        row.event_consensus_len > 0;
    const double statistical_local_fdr = ledger_row_statistical_te_local_fdr(row);
    if (row.mechanistic_ref_conflict_signal > 0.50 &&
        !direct_insert_observation &&
        statistical_local_fdr > kDefaultFinalConformalFdr) {
        return false;
    }
    const double structure_local_fdr = ledger_row_structure_te_local_fdr(row);
    return combine_independent_local_fdr(statistical_local_fdr, structure_local_fdr) <=
           kDefaultFinalConformalFdr;
}

bool has_pass_te_gate(const EvidenceLedgerRow& row) {
    const bool robust_pass =
        row.robust_mechanistic_qc == "PASS_TE_LFDR" &&
        row.robust_mechanistic_worst_case_lfdr <= 0.10;
    const bool structure_origin_pass = has_te_structure_origin_certificate(row);
    if (!robust_pass && !structure_origin_pass) {
        return false;
    }
    return row.posterior_qc == "PASS_TE_POSTERIOR" ||
           row.lfdr_qc == "PASS_TE_LFDR" ||
           row.robust_mechanistic_qc == "PASS_TE_LFDR" ||
           structure_origin_pass;
}

bool has_calibrated_te_statistical_support(const EvidenceLedgerRow& row) {
    const bool robust_pass =
        row.robust_mechanistic_qc == "PASS_TE_LFDR" &&
        row.robust_mechanistic_worst_case_lfdr <= 0.10;
    const bool structure_origin_pass = has_te_structure_origin_certificate(row);
    if (!robust_pass && !structure_origin_pass) {
        return false;
    }
    return row.posterior_qc == "PASS_TE_POSTERIOR" ||
           row.lfdr_qc == "PASS_TE_LFDR" ||
           robust_pass ||
           structure_origin_pass ||
           row.te_posterior >= 0.95 ||
           row.worst_case_lfdr <= 0.10;
}

double ledger_row_mechanistic_log_e_value(const EvidenceLedgerRow& row) {
    if (!has_conformal_feature_certificate(row)) {
        return -std::numeric_limits<double>::infinity();
    }
    if (row.mechanistic_ref_conflict_signal > 0.50) {
        return -std::numeric_limits<double>::infinity();
    }
    if (row.ref_span_reads > std::max(2, row.alt_struct_reads / 2)) {
        return -std::numeric_limits<double>::infinity();
    }
    const double lower_log_bf = std::min(
        row.mechanistic_lower_log_bf_te_vs_artifact,
        row.mechanistic_lower_log_bf_te_vs_non_te);
    return lower_log_bf - std::max(0.0, row.mechanistic_ambiguity_width);
}

bool has_event_evalue_support(const EvidenceLedgerRow& row) {
    const double log_e = ledger_row_mechanistic_log_e_value(row);
    return std::isfinite(log_e) && log_e > 0.0;
}

double ledger_row_event_existence_log_evidence(const EvidenceLedgerRow& row) {
    if (!has_conformal_feature_certificate(row)) {
        return -std::numeric_limits<double>::infinity();
    }
    const int32_t alt = std::max(0, row.alt_struct_reads);
    const int32_t ref = std::max(0, row.ref_span_reads);
    if (alt <= 0 || alt <= ref) {
        return -std::numeric_limits<double>::infinity();
    }
    if (row.mechanistic_ref_conflict_signal > 0.65) {
        return -std::numeric_limits<double>::infinity();
    }
    const double support_log_odds =
        std::log(static_cast<double>(alt + 1)) -
        std::log(static_cast<double>(ref + 1));
    const double structural_log_bf =
        row.mechanistic_lower_log_bf_te_vs_artifact -
        std::max(0.0, row.mechanistic_ambiguity_width);
    return support_log_odds +
           (0.50 * structural_log_bf) -
           (1.20 * std::clamp(row.mechanistic_ref_conflict_signal, 0.0, 1.0));
}

double low_allele_fraction_event_log_bf(int32_t alt, int32_t ref);

double competing_nonte_long_insertion_log_evidence(
    int32_t alt_struct_reads,
    int32_t ref_span_reads,
    int32_t raw_cigar_insert_reads,
    int32_t max_raw_cigar_insert_len,
    int32_t event_consensus_len,
    double non_te_posterior,
    double artifact_posterior,
    double ref_conflict_signal,
    double ambiguity_width) {
    const int32_t alt = std::max(0, alt_struct_reads);
    const int32_t ref = std::max(0, ref_span_reads);
    if (alt <= 0) {
        return -std::numeric_limits<double>::infinity();
    }
    const bool has_long_insert_observation =
        raw_cigar_insert_reads > 0 ||
        max_raw_cigar_insert_len > 0 ||
        event_consensus_len > 0;
    if (!has_long_insert_observation) {
        return -std::numeric_limits<double>::infinity();
    }
    const double posterior_log_odds =
        std::log(std::max(std::clamp(non_te_posterior, 0.0, 1.0), 1e-6)) -
        std::log(std::max(std::clamp(artifact_posterior, 0.0, 1.0), 1e-6));
    const double support_log_odds =
        std::log(static_cast<double>(alt + 1)) -
        std::log(static_cast<double>(ref + 1));
    const double raw_insert_log_bf = low_allele_fraction_event_log_bf(
        std::max(0, raw_cigar_insert_reads),
        ref);
    const double direct_insert_evidence = std::isfinite(raw_insert_log_bf)
        ? std::max(0.0, raw_insert_log_bf)
        : 0.0;
    const double assembled_insert_evidence =
        0.05 * std::log1p(static_cast<double>(std::max(0, max_raw_cigar_insert_len))) +
        0.03 * std::log1p(static_cast<double>(std::max(0, event_consensus_len)));
    return posterior_log_odds + support_log_odds + direct_insert_evidence +
           assembled_insert_evidence -
           std::max(0.0, ambiguity_width) -
           (0.35 * std::clamp(ref_conflict_signal, 0.0, 1.0));
}

bool has_competing_nonte_long_insertion_support(const EvidenceLedgerRow& row) {
    if (!has_conformal_feature_certificate(row)) {
        return false;
    }
    const double te_structure_local_fdr = ledger_row_structure_te_local_fdr(row);
    const bool te_origin_dominates =
        row.te_posterior > row.non_te_posterior &&
        row.te_posterior > row.artifact_posterior &&
        te_structure_local_fdr <= kDefaultFinalConformalFdr;
    if (te_origin_dominates) {
        return false;
    }
    const double log_evidence = competing_nonte_long_insertion_log_evidence(
        row.alt_struct_reads,
        row.ref_span_reads,
        row.raw_cigar_insert_reads,
        row.max_raw_cigar_insert_len,
        row.event_consensus_len,
        row.non_te_posterior,
        row.artifact_posterior,
        row.mechanistic_ref_conflict_signal,
        row.mechanistic_ambiguity_width);
    return std::isfinite(log_evidence) && log_evidence > 2.0;
}

double heterozygous_balance_log_evidence(int32_t alt, int32_t ref) {
    if (alt <= 0 || ref <= 0) {
        return 0.0;
    }
    const double depth = static_cast<double>(alt + ref);
    const double af = static_cast<double>(alt) / depth;
    const double balance = 1.0 - std::min(1.0, std::abs(af - 0.5) / 0.5);
    const double count_signal = std::log1p(static_cast<double>(std::min(alt, ref)));
    return balance * count_signal;
}

double log_choose_count(int32_t n, int32_t k) {
    if (k < 0 || n < 0 || k > n) {
        return -std::numeric_limits<double>::infinity();
    }
    return std::lgamma(static_cast<double>(n) + 1.0) -
           std::lgamma(static_cast<double>(k) + 1.0) -
           std::lgamma(static_cast<double>(n - k) + 1.0);
}

double beta_binomial_log_pmf(int32_t alt, int32_t total, double alpha, double beta) {
    if (alt < 0 || total < 0 || alt > total || alpha <= 0.0 || beta <= 0.0) {
        return -std::numeric_limits<double>::infinity();
    }
    return log_choose_count(total, alt) +
           std::lgamma(static_cast<double>(alt) + alpha) +
           std::lgamma(static_cast<double>(total - alt) + beta) -
           std::lgamma(static_cast<double>(total) + alpha + beta) +
           std::lgamma(alpha + beta) -
           std::lgamma(alpha) -
           std::lgamma(beta);
}

double binomial_log_pmf(int32_t alt, int32_t total, double p) {
    if (alt < 0 || total < 0 || alt > total) {
        return -std::numeric_limits<double>::infinity();
    }
    const double clamped_p = std::clamp(p, 1e-6, 1.0 - 1e-6);
    return log_choose_count(total, alt) +
           (static_cast<double>(alt) * std::log(clamped_p)) +
           (static_cast<double>(total - alt) * std::log1p(-clamped_p));
}

double low_allele_fraction_event_log_bf(int32_t alt, int32_t ref) {
    const int32_t total = std::max(0, alt) + std::max(0, ref);
    if (alt <= 0 || total <= 0) {
        return -std::numeric_limits<double>::infinity();
    }
    constexpr double kAltFractionPriorAlpha = 1.0;
    constexpr double kAltFractionPriorBeta = 9.0;
    constexpr double kArtifactErrorRate = 0.02;
    return beta_binomial_log_pmf(
               alt,
               total,
               kAltFractionPriorAlpha,
               kAltFractionPriorBeta) -
           binomial_log_pmf(alt, total, kArtifactErrorRate);
}

bool has_low_allele_fraction_event_certificate(const EvidenceLedgerRow& row) {
    if (!has_conformal_feature_certificate(row)) {
        return false;
    }
    const int32_t alt = std::max(0, row.alt_struct_reads);
    const int32_t ref = std::max(0, row.ref_span_reads);
    if (alt < 2 || ref <= 0) {
        return false;
    }
    const double af = static_cast<double>(alt) /
                      static_cast<double>(std::max(1, alt + ref));
    if (af >= 0.35) {
        return false;
    }
    const bool sequence_or_insert_signal =
        row.best_te_query_coverage >= 0.50 ||
        row.best_te_identity >= 0.60 ||
        row.cross_family_margin >= 0.015 ||
        row.final_qc == "TE_AMBIGUOUS" ||
        is_final_event_qc_for_conformal(row.final_qc);
    if (!sequence_or_insert_signal) {
        return false;
    }
    const double log_bf = low_allele_fraction_event_log_bf(alt, ref);
    return std::isfinite(log_bf) && log_bf >= 6.0;
}

bool has_self_stable_low_af_interval_certificate(const EvidenceLedgerRow& row) {
    if (!has_low_allele_fraction_event_certificate(row)) {
        return false;
    }
    if (row.bp_left < 0 || row.bp_right < 0) {
        return false;
    }
    const int32_t left = std::min(row.bp_left, row.bp_right);
    const int32_t right = std::max(row.bp_left, row.bp_right);
    if (left < 0 || right < 0 || right <= left) {
        return false;
    }
    if (row.mechanistic_blocks.find("boundary") == std::string::npos) {
        return false;
    }
    return right - left >= kFinalCallDedupDistanceBp;
}

bool has_ref_unopposed_bilateral_partial_anchor_certificate(
    const EvidenceLedgerRow& row) {
    if (!has_conformal_feature_certificate(row)) {
        return false;
    }
    if (row.final_qc != "NO_CALL_INCOMPLETE") {
        return false;
    }
    if (row.ref_span_reads != 0 || row.mechanistic_ref_conflict_signal > 1e-9) {
        return false;
    }
    if (row.alt_split_reads > 0 || row.alt_indel_reads > 0) {
        return false;
    }
    if (row.alt_struct_reads < 2 || row.support_qnames.size() < 2) {
        return false;
    }
    const bool bilateral_clip_evidence =
        row.alt_left_clip_reads > 0 && row.alt_right_clip_reads > 0;
    const bool bilateral_consensus_anchor =
        row.left_anchor_input_reads > 0 && row.right_anchor_input_reads > 0;
    if (!bilateral_clip_evidence || !bilateral_consensus_anchor) {
        return false;
    }
    if (row.partial_context_input_reads < 2 || row.input_event_reads < 2) {
        return false;
    }
    if (row.best_te_identity > 0.0 || row.best_te_query_coverage > 0.0 ||
        row.cross_family_margin > 0.0) {
        return false;
    }
    const int32_t left = std::min(row.bp_left, row.bp_right);
    const int32_t right = std::max(row.bp_left, row.bp_right);
    return left >= 0 && right >= left;
}

bool has_low_allele_fraction_event_certificate(const FinalCall& call) {
    if (!starts_with_text(call.final_qc, "PASS_STRUCTURAL_INSERTION")) {
        return false;
    }
    if (call.mechanistic_blocks.empty() || call.mechanistic_blocks == "NA") {
        return false;
    }
    const int32_t alt = std::max(0, call.alt_struct_reads);
    const int32_t ref = std::max(0, call.ref_span_reads);
    if (alt < 2 || ref <= 0) {
        return false;
    }
    const double af = static_cast<double>(alt) /
                      static_cast<double>(std::max(1, alt + ref));
    if (af >= 0.35) {
        return false;
    }
    const bool sequence_or_insert_signal =
        call.best_te_query_coverage >= 0.50 ||
        call.best_te_identity >= 0.60 ||
        call.cross_family_margin >= 0.015;
    if (!sequence_or_insert_signal) {
        return false;
    }
    const double log_bf = low_allele_fraction_event_log_bf(alt, ref);
    return std::isfinite(log_bf) && log_bf >= 6.0;
}

bool has_aggregated_low_allele_fraction_event_certificate(const FinalCall& call) {
    if (!starts_with_text(call.final_qc, "PASS_STRUCTURAL_INSERTION") ||
        call.final_qc.find("EVENT_COMMUNITY_STABLE") == std::string::npos ||
        call.final_qc.find("EVENT_COMMUNITY_AGGREGATED") == std::string::npos) {
        return false;
    }
    if (call.mechanistic_blocks.empty() || call.mechanistic_blocks == "NA" ||
        call.mechanistic_blocks.find("community") == std::string::npos ||
        call.mechanistic_blocks.find("read_overlap") == std::string::npos) {
        return false;
    }
    const int32_t alt = std::max(0, call.alt_struct_reads);
    const int32_t ref = std::max(0, call.ref_span_reads);
    if (alt < 2 || ref <= 0) {
        return false;
    }
    const double af = static_cast<double>(alt) /
                      static_cast<double>(std::max(1, alt + ref));
    if (af >= 0.35) {
        return false;
    }
    const double log_bf = low_allele_fraction_event_log_bf(alt, ref);
    return std::isfinite(log_bf) && log_bf >= 6.0;
}

bool has_self_stable_low_af_interval_certificate(const FinalCall& call) {
    if (!starts_with_text(call.final_qc, "PASS_STRUCTURAL_INSERTION") ||
        call.final_qc.find("EVENT_INTERVAL_STABLE") == std::string::npos) {
        return false;
    }
    if (call.mechanistic_blocks.empty() || call.mechanistic_blocks == "NA" ||
        call.mechanistic_blocks.find("boundary") == std::string::npos) {
        return false;
    }
    const int32_t left = std::min(call.bp_left, call.bp_right);
    const int32_t right = std::max(call.bp_left, call.bp_right);
    if (left < 0 || right <= left ||
        right - left < kFinalCallDedupDistanceBp) {
        return false;
    }
    const int32_t alt = std::max(0, call.alt_struct_reads);
    const int32_t ref = std::max(0, call.ref_span_reads);
    if (alt < 2 || ref <= 0) {
        return false;
    }
    const double af = static_cast<double>(alt) /
                      static_cast<double>(std::max(1, alt + ref));
    if (af >= 0.35) {
        return false;
    }
    const double log_bf = low_allele_fraction_event_log_bf(alt, ref);
    return std::isfinite(log_bf) && log_bf >= 6.0;
}

size_t sorted_support_intersection_size(
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
    return intersect;
}

std::vector<std::string> normalized_support_qnames(
    const std::vector<std::string>& values) {
    std::vector<std::string> normalized = values;
    std::sort(normalized.begin(), normalized.end());
    normalized.erase(
        std::unique(normalized.begin(), normalized.end()),
        normalized.end());
    return normalized;
}

double sorted_support_jaccard(
    const std::vector<std::string>& lhs,
    const std::vector<std::string>& rhs) {
    const size_t intersect = sorted_support_intersection_size(lhs, rhs);
    const size_t union_size = lhs.size() + rhs.size() - intersect;
    return union_size == 0
        ? 0.0
        : static_cast<double>(intersect) / static_cast<double>(union_size);
}

bool ledger_rows_share_support_edge(
    const EvidenceLedgerRow& a,
    const EvidenceLedgerRow& b) {
    if (a.support_qnames.empty() || b.support_qnames.empty()) {
        return false;
    }
    const size_t overlap = sorted_support_intersection_size(
        a.support_qnames,
        b.support_qnames);
    if (overlap == 0) {
        return false;
    }
    if (overlap >= 2) {
        return true;
    }
    return sorted_support_jaccard(a.support_qnames, b.support_qnames) >= 0.50;
}

bool final_calls_share_support_edge(const FinalCall& a, const FinalCall& b) {
    if (a.support_qnames.empty() || b.support_qnames.empty()) {
        return false;
    }
    const auto a_support = normalized_support_qnames(a.support_qnames);
    const auto b_support = normalized_support_qnames(b.support_qnames);
    const size_t overlap = sorted_support_intersection_size(
        a_support,
        b_support);
    if (overlap == 0) {
        return false;
    }
    if (overlap >= 2) {
        return true;
    }
    return sorted_support_jaccard(a_support, b_support) >= 0.50;
}

bool final_calls_share_identical_support(
    const FinalCall& a,
    const FinalCall& b) {
    if (a.support_qnames.empty() || b.support_qnames.empty()) {
        return false;
    }
    const auto a_support = normalized_support_qnames(a.support_qnames);
    const auto b_support = normalized_support_qnames(b.support_qnames);
    const size_t overlap = sorted_support_intersection_size(
        a_support,
        b_support);
    return overlap == a_support.size() && overlap == b_support.size();
}

bool final_calls_share_support_containment_edge(
    const FinalCall& a,
    const FinalCall& b) {
    if (a.support_qnames.empty() || b.support_qnames.empty()) {
        return false;
    }
    const auto a_support = normalized_support_qnames(a.support_qnames);
    const auto b_support = normalized_support_qnames(b.support_qnames);
    const size_t overlap = sorted_support_intersection_size(
        a_support,
        b_support);
    const size_t smaller_support = std::min(a_support.size(), b_support.size());
    return smaller_support > 0 && overlap * 2 >= smaller_support;
}

bool final_calls_have_overlapping_event_context(
    const FinalCall& a,
    const FinalCall& b) {
    if (a.tid != b.tid || a.chrom != b.chrom) {
        return false;
    }
    if (a.window_start < 0 || a.window_end < 0 ||
        b.window_start < 0 || b.window_end < 0) {
        return false;
    }
    const int32_t a_left = std::min(a.window_start, a.window_end);
    const int32_t a_right = std::max(a.window_start, a.window_end);
    const int32_t b_left = std::min(b.window_start, b.window_end);
    const int32_t b_right = std::max(b.window_start, b.window_end);
    return a_right >= b_left && b_right >= a_left;
}

bool final_calls_have_overlapping_breakpoint_intervals(
    const FinalCall& a,
    const FinalCall& b) {
    if (a.tid != b.tid || a.chrom != b.chrom) {
        return false;
    }
    const int32_t a_left = std::min(a.bp_left, a.bp_right);
    const int32_t a_right = std::max(a.bp_left, a.bp_right);
    const int32_t b_left = std::min(b.bp_left, b.bp_right);
    const int32_t b_right = std::max(b.bp_left, b.bp_right);
    if (a_left < 0 || a_right < 0 || b_left < 0 || b_right < 0) {
        return false;
    }
    return a_right >= b_left && b_right >= a_left;
}

bool ledger_row_is_structural_neighbor_evidence(const EvidenceLedgerRow& row) {
    if (!has_conformal_feature_certificate(row)) {
        return false;
    }
    if (row.alt_struct_reads <= 0) {
        return false;
    }
    const bool sequence_or_insert_signal =
        row.best_te_query_coverage >= 0.50 ||
        row.best_te_identity >= 0.60 ||
        row.cross_family_margin >= 0.015 ||
        row.final_qc == "TE_AMBIGUOUS" ||
        is_final_event_qc_for_conformal(row.final_qc);
    if (!sequence_or_insert_signal) {
        return false;
    }
    return row.mechanistic_ref_conflict_signal <= 0.95;
}

bool mechanistic_blocks_include(
    const EvidenceLedgerRow& row,
    const char* token) {
    return row.mechanistic_blocks.find(token) != std::string::npos;
}

bool ledger_row_is_aggregatable_event_community_evidence(
    const EvidenceLedgerRow& row) {
    if (!has_conformal_feature_certificate(row)) {
        return false;
    }
    if (row.pos < 0 && row.bp_left < 0 && row.bp_right < 0) {
        return false;
    }
    if (row.alt_struct_reads <= 0 || row.support_qnames.empty()) {
        return false;
    }
    if (is_final_te_qc_for_conformal(row.final_qc)) {
        return false;
    }
    const bool interpretable_event =
        mechanistic_blocks_include(row, "event") ||
        row.final_qc == "TE_AMBIGUOUS" ||
        row.final_qc == "REFERENCE_OR_ARTIFACT";
    if (!interpretable_event) {
        return false;
    }
    return row.final_qc == "REFERENCE_OR_ARTIFACT" ||
           row.final_qc == "TE_AMBIGUOUS" ||
           row.final_qc == "PASS_NONTE_INSERTION";
}

int32_t ledger_row_breakpoint_event_start(const EvidenceLedgerRow& row);
int32_t ledger_row_breakpoint_event_end(const EvidenceLedgerRow& row);

bool ledger_rows_are_spatial_event_neighbors(
    const EvidenceLedgerRow& a,
    const EvidenceLedgerRow& b) {
    if (a.tid != b.tid || a.chrom != b.chrom) {
        return false;
    }
    const int32_t a_start = ledger_row_breakpoint_event_start(a);
    const int32_t a_end = ledger_row_breakpoint_event_end(a);
    const int32_t b_start = ledger_row_breakpoint_event_start(b);
    const int32_t b_end = ledger_row_breakpoint_event_end(b);
    if (a_start < 0 || a_end < 0 || b_start < 0 || b_end < 0) {
        return false;
    }
    if (a_end < b_start) {
        return b_start - a_end <= kFinalCallDedupDistanceBp;
    }
    if (b_end < a_start) {
        return a_start - b_end <= kFinalCallDedupDistanceBp;
    }
    return true;
}

bool event_cluster_has_stable_low_af_support(
    size_t focal_index,
    const std::vector<EvidenceLedgerRow>& rows) {
    const auto& focal = rows[focal_index];
    if (!has_low_allele_fraction_event_certificate(focal)) {
        return false;
    }
    for (size_t idx = 0; idx < rows.size(); ++idx) {
        if (idx == focal_index) {
            continue;
        }
        const auto& neighbor = rows[idx];
        if (!ledger_rows_are_spatial_event_neighbors(focal, neighbor)) {
            continue;
        }
        if (!ledger_row_is_structural_neighbor_evidence(neighbor)) {
            continue;
        }
        if (ledger_rows_share_support_edge(focal, neighbor)) {
            return true;
        }
    }
    return false;
}

bool has_event_existence_support(const EvidenceLedgerRow& row) {
    const double log_evidence = ledger_row_event_existence_log_evidence(row);
    return std::isfinite(log_evidence) && log_evidence > 2.0;
}

bool has_balanced_heterozygous_event_support(const EvidenceLedgerRow& row) {
    if (!has_conformal_feature_certificate(row)) {
        return false;
    }
    const int32_t alt = std::max(0, row.alt_struct_reads);
    const int32_t ref = std::max(0, row.ref_span_reads);
    if (alt < 16 || alt <= ref) {
        return false;
    }
    const double balance = heterozygous_balance_log_evidence(alt, ref);
    if (balance < 1.75 || row.mechanistic_ref_conflict_signal > 0.65) {
        return false;
    }
    const bool event_interpretable =
        row.final_qc == "TE_AMBIGUOUS" ||
        row.final_qc == "PASS_NONTE_INSERTION" ||
        is_final_event_qc_for_conformal(row.final_qc);
    if (!event_interpretable) {
        return false;
    }
    const bool sequence_or_event_signal =
        row.best_te_query_coverage >= 0.80 ||
        row.best_te_identity >= 0.54 ||
        row.cross_family_margin >= 0.015 ||
        row.mechanistic_lower_log_bf_te_vs_artifact > -1.0;
    return sequence_or_event_signal;
}

bool ledger_row_is_final_te(const EvidenceLedgerRow& row) {
    return is_final_te_qc_for_conformal(row.final_qc);
}

bool ledger_row_is_promotable_event_evidence(const EvidenceLedgerRow& row) {
    if (!has_conformal_feature_certificate(row)) {
        return false;
    }
    if (row.pos < 0 && row.bp_left < 0 && row.bp_right < 0) {
        return false;
    }
    if (row.alt_struct_reads <= 0) {
        return false;
    }
    const bool low_reference_conflict =
        row.mechanistic_ref_conflict_signal <= 0.50 &&
        row.ref_span_reads <= std::max(2, row.alt_struct_reads / 2);
    const bool balanced_heterozygous_event =
        has_balanced_heterozygous_event_support(row);
    const bool low_allele_fraction_event =
        has_low_allele_fraction_event_certificate(row);
    const bool bilateral_partial_anchor_event =
        has_ref_unopposed_bilateral_partial_anchor_certificate(row);
    const bool competing_nonte_event =
        has_competing_nonte_long_insertion_support(row);
    if (!low_reference_conflict &&
        !balanced_heterozygous_event &&
        !low_allele_fraction_event &&
        !bilateral_partial_anchor_event &&
        !competing_nonte_event &&
        !has_calibrated_te_statistical_support(row)) {
        return false;
    }
    const bool statistical_te_signal =
        is_final_event_qc_for_conformal(row.final_qc) ||
        has_pass_te_gate(row) ||
        row.te_posterior >= 0.50 ||
        row.lfdr <= 0.50 ||
        row.worst_case_lfdr <= 0.50 ||
        row.robust_mechanistic_worst_case_lfdr <= 0.50 ||
        has_event_evalue_support(row);
    const bool structural_event_signal =
        has_event_existence_support(row) ||
        balanced_heterozygous_event ||
        low_allele_fraction_event ||
        bilateral_partial_anchor_event ||
        competing_nonte_event;
    const bool sequence_signal =
        row.best_te_identity > 0.0 ||
        row.best_te_query_coverage > 0.0 ||
        row.cross_family_margin > 0.0 ||
        row.family == "UNKNOWN" ||
        row.subfamily == "UNKNOWN";
    return (statistical_te_signal && sequence_signal) ||
           structural_event_signal;
}

bool ledger_row_requires_stable_event_community(const EvidenceLedgerRow& row) {
    return has_low_allele_fraction_event_certificate(row) &&
           !has_self_stable_low_af_interval_certificate(row) &&
           !has_balanced_heterozygous_event_support(row) &&
           !has_event_existence_support(row) &&
           !has_ref_unopposed_bilateral_partial_anchor_certificate(row) &&
           !has_calibrated_te_statistical_support(row);
}

// L4: discrete local-context bucket for context-conditioned null calibration.
// One-sided events (exactly one flank anchored) are structurally artifact-prone
// and are calibrated against one-sided null controls; two-sided (closed) or
// flank-less events use the other bucket. Flank lengths are the context signal
// available on both the evidence-ledger null controls and the final calls.
int32_t conformal_context_bucket(int32_t left_flank_align_len, int32_t right_flank_align_len) {
    const bool has_left = left_flank_align_len > 0;
    const bool has_right = right_flank_align_len > 0;
    return (has_left != has_right) ? 1 : 0;
}

ConformalFeatureVector conformal_feature_from_call(const FinalCall& call, size_t index) {
    ConformalFeatureVector feature;
    feature.id = call.chrom + ":" + std::to_string(call.pos) + ":" +
                 std::to_string(call.bp_left) + ":" +
                 std::to_string(call.bp_right) + ":" + std::to_string(index);
    feature.pro_te = {
        static_cast<double>(std::max(0, call.alt_struct_reads)),
        std::max(0.0, call.best_te_identity),
        std::max(0.0, call.best_te_query_coverage),
        std::max(0.0, call.cross_family_margin),
    };
    feature.ref_span_reads = static_cast<double>(std::max(0, call.ref_span_reads));
    feature.context = conformal_context_bucket(
        call.left_flank_align_len, call.right_flank_align_len);
    return feature;
}

ConformalFeatureVector conformal_feature_from_row(const EvidenceLedgerRow& row, size_t index) {
    ConformalFeatureVector feature;
    feature.id = row.chrom + ":" + std::to_string(row.pos) + ":" +
                 std::to_string(row.bp_left) + ":" +
                 std::to_string(row.bp_right) + ":" + std::to_string(index);
    feature.pro_te = {
        static_cast<double>(std::max(0, row.alt_struct_reads)),
        std::max(0.0, row.best_te_identity),
        std::max(0.0, row.best_te_query_coverage),
        std::max(0.0, row.cross_family_margin),
    };
    feature.ref_span_reads = static_cast<double>(std::max(0, row.ref_span_reads));
    feature.context = conformal_context_bucket(
        row.left_flank_align_len, row.right_flank_align_len);
    return feature;
}

int32_t ledger_row_breakpoint_event_start(const EvidenceLedgerRow& row) {
    if (row.bp_left >= 0 && row.bp_right >= 0) {
        return std::min(row.bp_left, row.bp_right);
    }
    if (row.pos >= 0) {
        return row.pos;
    }
    return -1;
}

int32_t ledger_row_breakpoint_event_end(const EvidenceLedgerRow& row) {
    if (row.bp_left >= 0 && row.bp_right >= 0) {
        return std::max(row.bp_left, row.bp_right);
    }
    if (row.pos >= 0) {
        return row.pos;
    }
    return -1;
}

int32_t ledger_row_breakpoint_event_center(const EvidenceLedgerRow& row) {
    const int32_t start = ledger_row_breakpoint_event_start(row);
    const int32_t end = ledger_row_breakpoint_event_end(row);
    if (start >= 0 && end >= 0) {
        return start + ((end - start) / 2);
    }
    return row.pos;
}

int32_t ledger_row_context_start(const EvidenceLedgerRow& row) {
    if (row.coverage_left >= 0 && row.coverage_right >= 0) {
        return std::min(row.coverage_left, row.coverage_right);
    }
    return ledger_row_breakpoint_event_start(row);
}

int32_t ledger_row_context_end(const EvidenceLedgerRow& row) {
    if (row.coverage_left >= 0 && row.coverage_right >= 0) {
        return std::max(row.coverage_left, row.coverage_right);
    }
    return ledger_row_breakpoint_event_end(row);
}

bool is_sample_local_null_control(const EvidenceLedgerRow& row) {
    if (!has_conformal_feature_certificate(row)) {
        return false;
    }
    if (ledger_row_is_promotable_event_evidence(row)) {
        return false;
    }
    return !is_final_te_qc_for_conformal(row.final_qc);
}

bool ledger_rows_share_event_cluster(
    const EvidenceLedgerRow& a,
    const EvidenceLedgerRow& b) {
    if (a.tid != b.tid || a.chrom != b.chrom) {
        return false;
    }
    const int32_t a_start = ledger_row_breakpoint_event_start(a);
    const int32_t a_end = ledger_row_breakpoint_event_end(a);
    const int32_t b_start = ledger_row_breakpoint_event_start(b);
    const int32_t b_end = ledger_row_breakpoint_event_end(b);
    if (a_start < 0 || a_end < 0 || b_start < 0 || b_end < 0) {
        return false;
    }
    if (a_end < b_start) {
        return b_start - a_end <= kFinalCallDedupDistanceBp;
    }
    if (b_end < a_start) {
        return a_start - b_end <= kFinalCallDedupDistanceBp;
    }
    return true;
}

bool ledger_rows_have_overlapping_breakpoint_intervals(
    const EvidenceLedgerRow& a,
    const EvidenceLedgerRow& b) {
    if (a.tid != b.tid || a.chrom != b.chrom) {
        return false;
    }
    const int32_t a_start = ledger_row_breakpoint_event_start(a);
    const int32_t a_end = ledger_row_breakpoint_event_end(a);
    const int32_t b_start = ledger_row_breakpoint_event_start(b);
    const int32_t b_end = ledger_row_breakpoint_event_end(b);
    if (a_start < 0 || a_end < 0 || b_start < 0 || b_end < 0) {
        return false;
    }
    return a_end >= b_start && b_end >= a_start;
}

bool ledger_rows_share_promotable_event_cluster(
    const EvidenceLedgerRow& a,
    const EvidenceLedgerRow& b) {
    if (!ledger_rows_share_event_cluster(a, b)) {
        return false;
    }
    const bool a_stable_interval = has_self_stable_low_af_interval_certificate(a);
    const bool b_stable_interval = has_self_stable_low_af_interval_certificate(b);
    if (!a_stable_interval && !b_stable_interval) {
        return true;
    }
    return ledger_rows_have_overlapping_breakpoint_intervals(a, b) ||
           ledger_rows_share_support_edge(a, b);
}

double ledger_row_promotion_score(const EvidenceLedgerRow& row) {
    double score = 0.0;
    score += static_cast<double>(std::max(0, row.alt_struct_reads));
    score += std::clamp(row.te_posterior, 0.0, 1.0) * 8.0;
    score += (1.0 - std::clamp(row.worst_case_lfdr, 0.0, 1.0)) * 4.0;
    score += (1.0 - std::clamp(row.robust_mechanistic_worst_case_lfdr, 0.0, 1.0)) * 6.0;
    score += std::max(0.0, row.best_te_identity) * 3.0;
    score += std::max(0.0, row.best_te_query_coverage) * 3.0;
    score += std::max(0.0, row.cross_family_margin) * 2.0;
    score -= static_cast<double>(std::max(0, row.ref_span_reads)) * 1.5;
    score -= std::clamp(row.mechanistic_ref_conflict_signal, 0.0, 1.0) * 10.0;
    if (row.robust_mechanistic_qc == "PASS_TE_LFDR") {
        score += 6.0;
    }
    if (row.lfdr_qc == "PASS_TE_LFDR") {
        score += 4.0;
    }
    if (row.posterior_qc == "PASS_TE_POSTERIOR") {
        score += 2.0;
    }
    return score;
}

bool prefer_ledger_row_for_promotion(
    const EvidenceLedgerRow& candidate,
    const EvidenceLedgerRow& incumbent) {
    const double candidate_score = ledger_row_promotion_score(candidate);
    const double incumbent_score = ledger_row_promotion_score(incumbent);
    if (candidate_score != incumbent_score) {
        return candidate_score > incumbent_score;
    }
    if (candidate.alt_struct_reads != incumbent.alt_struct_reads) {
        return candidate.alt_struct_reads > incumbent.alt_struct_reads;
    }
    if (candidate.robust_mechanistic_worst_case_lfdr !=
        incumbent.robust_mechanistic_worst_case_lfdr) {
        return candidate.robust_mechanistic_worst_case_lfdr <
               incumbent.robust_mechanistic_worst_case_lfdr;
    }
    if (candidate.ref_span_reads != incumbent.ref_span_reads) {
        return candidate.ref_span_reads < incumbent.ref_span_reads;
    }
    return ledger_row_breakpoint_event_center(candidate) <
           ledger_row_breakpoint_event_center(incumbent);
}

bool final_call_overlaps_ledger_event(
    const FinalCall& call,
    const EvidenceLedgerRow& row) {
    if (call.tid != row.tid || call.chrom != row.chrom) {
        return false;
    }
    const int32_t call_left = std::min(call.bp_left, call.bp_right);
    const int32_t call_right = std::max(call.bp_left, call.bp_right);
    const int32_t row_left = ledger_row_breakpoint_event_start(row);
    const int32_t row_right = ledger_row_breakpoint_event_end(row);
    if (call_left < 0 || call_right < 0 || row_left < 0 || row_right < 0) {
        return call.pos >= 0 && ledger_row_breakpoint_event_center(row) >= 0 &&
               std::abs(call.pos - ledger_row_breakpoint_event_center(row)) <=
                   kFinalCallDedupDistanceBp;
    }
    if (call_right < row_left) {
        return row_left - call_right <= kFinalCallDedupDistanceBp;
    }
    if (row_right < call_left) {
        return call_left - row_right <= kFinalCallDedupDistanceBp;
    }
    return true;
}

FinalCall promoted_call_from_ledger_row(const EvidenceLedgerRow& row) {
    FinalCall call;
    const bool bilateral_partial_anchor_event =
        has_ref_unopposed_bilateral_partial_anchor_certificate(row);
    call.chrom = row.chrom;
    call.tid = row.tid;
    call.bp_left = row.bp_left >= 0 ? row.bp_left : ledger_row_breakpoint_event_center(row);
    call.bp_right = row.bp_right >= 0 ? row.bp_right : call.bp_left;
    if (call.bp_left > call.bp_right) {
        std::swap(call.bp_left, call.bp_right);
    }
    call.pos = ledger_row_breakpoint_event_center(row);
    if (call.pos < 0) {
        call.pos = call.bp_left;
    }
    call.window_start = ledger_row_context_start(row);
    call.window_end = ledger_row_context_end(row);
    if (!ledger_row_is_final_te(row) &&
        row.owner_context_left >= 0 && row.owner_context_right >= 0) {
        call.window_start = std::min(row.owner_context_left, row.owner_context_right);
        call.window_end = std::max(row.owner_context_left, row.owner_context_right);
    }
    if (!ledger_row_is_final_te(row) && row.event_consensus_len > 0 && call.pos >= 0) {
        const int32_t uncertainty = row.event_consensus_len / 2;
        call.window_start = std::min(
            call.window_start >= 0 ? call.window_start : call.pos,
            std::max(0, call.pos - uncertainty));
        call.window_end = std::max(
            call.window_end >= 0 ? call.window_end : call.pos,
            call.pos + uncertainty);
    }
    call.family = row.family.empty() ? "UNKNOWN" : row.family;
    if (call.family == "NA") {
        call.family = "UNKNOWN";
    }
    call.subfamily = row.subfamily.empty() ? "UNKNOWN" : row.subfamily;
    if (call.subfamily == "NA") {
        call.subfamily = call.family;
    }
    call.te_name = call.subfamily != "NA" ? call.subfamily : call.family;
    call.family_committed =
        row.family_alignment_resolved &&
        !(call.family == "UNKNOWN" && call.subfamily == "UNKNOWN");
    call.alt_struct_reads = row.alt_struct_reads;
    call.support_reads = row.alt_struct_reads;
    call.raw_cigar_insert_reads = row.raw_cigar_insert_reads;
    call.max_raw_cigar_insert_len = row.max_raw_cigar_insert_len;
    call.ref_span_reads = row.ref_span_reads;
    call.support_qnames = normalized_support_qnames(row.support_qnames);
    call.best_te_identity = row.best_te_identity;
    call.best_te_query_coverage = row.best_te_query_coverage;
    call.cross_family_margin = row.cross_family_margin;
    call.te_structure_path = row.te_structure_path;
    call.te_structure_log_evidence = row.te_structure_log_evidence;
    call.nonte_structure_log_evidence = row.nonte_structure_log_evidence;
    call.artifact_structure_log_evidence = row.artifact_structure_log_evidence;
    call.te_structure_path_confidence = row.te_structure_path_confidence;
    call.polyA_posterior = row.polyA_posterior;
    call.transduction_posterior = row.transduction_posterior;
    call.te_core_coverage = row.te_core_coverage;
    call.unexplained_high_complexity_bp = row.unexplained_high_complexity_bp;
    call.event_consensus_len = row.event_consensus_len;
    call.left_flank_align_len = row.left_flank_align_len;
    call.right_flank_align_len = row.right_flank_align_len;
    call.insert_seq = row.insert_seq;
    call.final_qc = "PASS_TE_IMPRECISE|EVENT_CLUSTER_PROMOTED";
    if (!ledger_row_is_final_te(row) &&
        (has_event_existence_support(row) ||
         has_balanced_heterozygous_event_support(row) ||
         has_low_allele_fraction_event_certificate(row) ||
         has_competing_nonte_long_insertion_support(row) ||
         bilateral_partial_anchor_event)) {
        call.family = "UNKNOWN";
        call.subfamily = "UNKNOWN";
        call.te_name = "UNKNOWN";
        call.family_committed = false;
        call.final_qc = "PASS_STRUCTURAL_INSERTION|EVENT_CLUSTER_PROMOTED";
    }
    if (bilateral_partial_anchor_event) {
        const int32_t left = std::min(row.bp_left, row.bp_right);
        const int32_t right = std::max(row.bp_left, row.bp_right);
        const int32_t uncertainty = std::max(
            0,
            std::max(row.event_consensus_len, right - left));
        call.bp_left = std::max(0, left - uncertainty);
        call.bp_right = right + uncertainty;
        call.final_qc += "|EVENT_BILATERAL_PARTIAL_ANCHOR";
    }
    call.te_qc = "EVENT_CLUSTER_LEDGER_EVIDENCE";
    call.best_explanation = "EVENT_CLUSTER";
    call.explanation_residual = "PROMOTED_FROM_LEDGER";
    call.explanation_path = "NA";
    call.te_posterior = row.te_posterior;
    call.non_te_posterior = row.non_te_posterior;
    call.artifact_posterior = row.artifact_posterior;
    call.lfdr = row.lfdr;
    call.worst_case_lfdr = row.worst_case_lfdr;
    call.posterior_qc = row.posterior_qc;
    call.lfdr_qc = row.lfdr_qc;
    call.mechanistic_lower_log_bf_te_vs_artifact =
        row.mechanistic_lower_log_bf_te_vs_artifact;
    call.mechanistic_lower_log_bf_te_vs_non_te =
        row.mechanistic_lower_log_bf_te_vs_non_te;
    call.mechanistic_ref_conflict_signal = row.mechanistic_ref_conflict_signal;
    call.mechanistic_ambiguity_width = row.mechanistic_ambiguity_width;
    call.mechanistic_blocks = row.mechanistic_blocks;
    call.robust_mechanistic_lfdr = row.robust_mechanistic_lfdr;
    call.robust_mechanistic_worst_case_lfdr =
        row.robust_mechanistic_worst_case_lfdr;
    call.robust_mechanistic_qc = row.robust_mechanistic_qc;
    return call;
}

void annotate_stable_event_community_certificate(FinalCall& call) {
    if (call.final_qc.find("EVENT_COMMUNITY_STABLE") == std::string::npos) {
        call.final_qc += "|EVENT_COMMUNITY_STABLE";
    }
}

void annotate_stable_event_interval_certificate(FinalCall& call) {
    if (call.final_qc.find("EVENT_INTERVAL_STABLE") == std::string::npos) {
        call.final_qc += "|EVENT_INTERVAL_STABLE";
    }
}

void annotate_aggregated_event_community_certificate(FinalCall& call) {
    if (call.final_qc.find("EVENT_COMMUNITY_AGGREGATED") == std::string::npos) {
        call.final_qc += "|EVENT_COMMUNITY_AGGREGATED";
    }
}

std::vector<std::string> merged_support_qnames_for_component(
    const std::vector<size_t>& component,
    const std::vector<EvidenceLedgerRow>& rows) {
    std::vector<std::string> support;
    for (const size_t idx : component) {
        support.insert(
            support.end(),
            rows[idx].support_qnames.begin(),
            rows[idx].support_qnames.end());
    }
    std::sort(support.begin(), support.end());
    support.erase(std::unique(support.begin(), support.end()), support.end());
    return support;
}

int32_t median_ref_span_for_component(
    const std::vector<size_t>& component,
    const std::vector<EvidenceLedgerRow>& rows) {
    std::vector<int32_t> ref_spans;
    ref_spans.reserve(component.size());
    for (const size_t idx : component) {
        ref_spans.push_back(std::max(0, rows[idx].ref_span_reads));
    }
    if (ref_spans.empty()) {
        return 0;
    }
    std::sort(ref_spans.begin(), ref_spans.end());
    return ref_spans[ref_spans.size() / 2];
}

EvidenceLedgerRow aggregate_event_community_row(
    const std::vector<size_t>& component,
    const std::vector<EvidenceLedgerRow>& rows) {
    EvidenceLedgerRow aggregate = rows[component.front()];
    const auto support = merged_support_qnames_for_component(component, rows);
    int32_t left = std::numeric_limits<int32_t>::max();
    int32_t right = -1;
    int32_t context_left = std::numeric_limits<int32_t>::max();
    int32_t context_right = -1;
    int32_t best_center = ledger_row_breakpoint_event_center(rows[component.front()]);
    double best_identity = 0.0;
    double best_coverage = 0.0;
    double best_margin = 0.0;
    double best_te_posterior = 0.0;
    double best_lfdr = 1.0;
    double best_worst_case_lfdr = 1.0;
    double min_ref_conflict = 1.0;
    double max_ambiguity = 0.0;
    int32_t raw_cigar_insert_reads = 0;
    int32_t max_raw_cigar_insert_len = 0;
    int32_t max_left_flank_align_len = 0;
    int32_t max_right_flank_align_len = 0;
    std::string representative_insert_seq;

    for (const size_t idx : component) {
        const auto& row = rows[idx];
        const int32_t row_left = ledger_row_breakpoint_event_start(row);
        const int32_t row_right = ledger_row_breakpoint_event_end(row);
        if (row_left >= 0) {
            left = std::min(left, row_left);
        }
        if (row_right >= 0) {
            right = std::max(right, row_right);
        }
        const int32_t row_context_left = ledger_row_context_start(row);
        const int32_t row_context_right = ledger_row_context_end(row);
        if (row_context_left >= 0) {
            context_left = std::min(context_left, row_context_left);
        }
        if (row_context_right >= 0) {
            context_right = std::max(context_right, row_context_right);
        }
        best_identity = std::max(best_identity, row.best_te_identity);
        best_coverage = std::max(best_coverage, row.best_te_query_coverage);
        best_margin = std::max(best_margin, row.cross_family_margin);
        best_te_posterior = std::max(best_te_posterior, row.te_posterior);
        best_lfdr = std::min(best_lfdr, row.lfdr);
        best_worst_case_lfdr = std::min(best_worst_case_lfdr, row.worst_case_lfdr);
        min_ref_conflict =
            std::min(min_ref_conflict, row.mechanistic_ref_conflict_signal);
        raw_cigar_insert_reads += std::max(0, row.raw_cigar_insert_reads);
        max_raw_cigar_insert_len = std::max(
            max_raw_cigar_insert_len,
            row.max_raw_cigar_insert_len);
        max_left_flank_align_len = std::max(
            max_left_flank_align_len,
            row.left_flank_align_len);
        max_right_flank_align_len = std::max(
            max_right_flank_align_len,
            row.right_flank_align_len);
        if (row.insert_seq.size() > representative_insert_seq.size()) {
            representative_insert_seq = row.insert_seq;
        }
        max_ambiguity =
            std::max(max_ambiguity, row.mechanistic_ambiguity_width);
    }

    if (left == std::numeric_limits<int32_t>::max()) {
        left = best_center;
    }
    if (right < 0) {
        right = left;
    }
    if (context_left == std::numeric_limits<int32_t>::max()) {
        context_left = left;
    }
    if (context_right < 0) {
        context_right = right;
    }
    best_center = left + ((right - left) / 2);

    aggregate.pos = best_center;
    aggregate.bp_left = std::max(0, left - kFinalCallDedupDistanceBp);
    aggregate.bp_right = right + kFinalCallDedupDistanceBp;
    aggregate.coverage_left = std::min(context_left, aggregate.bp_left);
    aggregate.coverage_right = std::max(context_right, aggregate.bp_right);
    aggregate.family = "UNKNOWN";
    aggregate.subfamily = "UNKNOWN";
    aggregate.final_qc =
        "PASS_STRUCTURAL_INSERTION|EVENT_CLUSTER_PROMOTED|"
        "EVENT_COMMUNITY_STABLE|EVENT_COMMUNITY_AGGREGATED";
    aggregate.posterior_qc = "TE_POSTERIOR_LOW";
    aggregate.lfdr_qc = "TE_LFDR_HIGH";
    aggregate.candidate_retention_reason = "EVENT_COMMUNITY_AGGREGATE";
    aggregate.alt_struct_reads =
        static_cast<int32_t>(std::min<size_t>(
            support.size(),
            static_cast<size_t>(std::numeric_limits<int32_t>::max())));
    aggregate.raw_cigar_insert_reads = raw_cigar_insert_reads;
    aggregate.max_raw_cigar_insert_len = max_raw_cigar_insert_len;
    aggregate.left_flank_align_len = max_left_flank_align_len;
    aggregate.right_flank_align_len = max_right_flank_align_len;
    aggregate.insert_seq = representative_insert_seq;
    aggregate.ref_span_reads = median_ref_span_for_component(component, rows);
    aggregate.support_qnames = support;
    aggregate.best_te_identity = best_identity;
    aggregate.best_te_query_coverage = best_coverage;
    aggregate.cross_family_margin = best_margin;
    aggregate.te_posterior = best_te_posterior;
    aggregate.lfdr = best_lfdr;
    aggregate.worst_case_lfdr = best_worst_case_lfdr;
    aggregate.mechanistic_lower_log_bf_te_vs_artifact = 0.0;
    aggregate.mechanistic_lower_log_bf_te_vs_non_te = 0.0;
    aggregate.mechanistic_ref_conflict_signal = min_ref_conflict;
    aggregate.mechanistic_ambiguity_width = max_ambiguity;
    aggregate.mechanistic_blocks =
        "event;community;read_overlap;boundary;ref_conflict";
    aggregate.robust_mechanistic_lfdr = 1.0;
    aggregate.robust_mechanistic_worst_case_lfdr = 1.0;
    aggregate.robust_mechanistic_qc = "TE_LFDR_HIGH";
    return aggregate;
}

bool aggregate_event_community_has_low_af_certificate(
    const EvidenceLedgerRow& row) {
    if (row.final_qc.find("EVENT_COMMUNITY_AGGREGATED") == std::string::npos ||
        row.final_qc.find("EVENT_COMMUNITY_STABLE") == std::string::npos) {
        return false;
    }
    return has_low_allele_fraction_event_certificate(row);
}

FinalCall aggregated_event_community_call_from_row(
    const EvidenceLedgerRow& row) {
    FinalCall call = promoted_call_from_ledger_row(row);
    call.family = "UNKNOWN";
    call.subfamily = "UNKNOWN";
    call.te_name = "UNKNOWN";
    call.final_qc =
        "PASS_STRUCTURAL_INSERTION|EVENT_CLUSTER_PROMOTED|"
        "EVENT_COMMUNITY_STABLE|EVENT_COMMUNITY_AGGREGATED";
    call.te_qc = "EVENT_COMMUNITY_AGGREGATE";
    call.best_explanation = "EVENT_COMMUNITY";
    call.explanation_residual = "AGGREGATED_WEAK_LEDGER_EVIDENCE";
    call.explanation_path = "READ_OVERLAP_COMMUNITY";
    return call;
}

void append_aggregated_event_community_calls(PipelineResult& result) {
    std::vector<size_t> candidate_indices;
    candidate_indices.reserve(result.evidence_ledger.size());
    for (size_t i = 0; i < result.evidence_ledger.size(); ++i) {
        if (ledger_row_is_aggregatable_event_community_evidence(
                result.evidence_ledger[i])) {
            candidate_indices.push_back(i);
        }
    }
    if (candidate_indices.size() < 2) {
        return;
    }

    std::sort(candidate_indices.begin(), candidate_indices.end(), [&](size_t lhs, size_t rhs) {
        const auto& a = result.evidence_ledger[lhs];
        const auto& b = result.evidence_ledger[rhs];
        if (a.tid != b.tid) {
            return a.tid < b.tid;
        }
        if (a.chrom != b.chrom) {
            return a.chrom < b.chrom;
        }
        const int32_t a_start = ledger_row_breakpoint_event_start(a);
        const int32_t b_start = ledger_row_breakpoint_event_start(b);
        if (a_start != b_start) {
            return a_start < b_start;
        }
        return lhs < rhs;
    });

    std::vector<size_t> parent(candidate_indices.size());
    for (size_t i = 0; i < parent.size(); ++i) {
        parent[i] = i;
    }
    const auto find_root = [&parent](size_t value) {
        size_t root = value;
        while (parent[root] != root) {
            root = parent[root];
        }
        while (parent[value] != value) {
            const size_t next = parent[value];
            parent[value] = root;
            value = next;
        }
        return root;
    };
    auto unite = [&parent, &find_root](size_t lhs, size_t rhs) {
        const size_t lhs_root = find_root(lhs);
        const size_t rhs_root = find_root(rhs);
        if (lhs_root != rhs_root) {
            parent[rhs_root] = lhs_root;
        }
    };

    for (size_t i = 0; i < candidate_indices.size(); ++i) {
        const auto& a = result.evidence_ledger[candidate_indices[i]];
        const int32_t a_end = ledger_row_breakpoint_event_end(a);
        for (size_t j = i + 1; j < candidate_indices.size(); ++j) {
            const auto& b = result.evidence_ledger[candidate_indices[j]];
            if (a.tid != b.tid || a.chrom != b.chrom) {
                break;
            }
            const int32_t b_start = ledger_row_breakpoint_event_start(b);
            if (a_end >= 0 && b_start >= 0 &&
                b_start - a_end > kFinalCallDedupDistanceBp) {
                break;
            }
            if (ledger_rows_are_spatial_event_neighbors(a, b) &&
                ledger_rows_share_support_edge(a, b)) {
                unite(i, j);
            }
        }
    }

    std::unordered_map<size_t, std::vector<size_t>> components;
    components.reserve(candidate_indices.size());
    for (size_t i = 0; i < candidate_indices.size(); ++i) {
        components[find_root(i)].push_back(candidate_indices[i]);
    }

    for (const auto& [unused_root, component] : components) {
        (void)unused_root;
        if (component.size() < 2) {
            continue;
        }
        bool contains_promotable_row = false;
        for (const size_t idx : component) {
            if (ledger_row_is_promotable_event_evidence(
                    result.evidence_ledger[idx])) {
                contains_promotable_row = true;
                break;
            }
        }
        if (contains_promotable_row) {
            continue;
        }

        const EvidenceLedgerRow aggregate =
            aggregate_event_community_row(component, result.evidence_ledger);
        if (!aggregate_event_community_has_low_af_certificate(aggregate)) {
            continue;
        }
        result.final_calls.push_back(
            aggregated_event_community_call_from_row(aggregate));
    }
}

void append_promoted_event_cluster_calls(PipelineResult& result) {
    std::vector<size_t> promotable_indices;
    promotable_indices.reserve(result.evidence_ledger.size());
    for (size_t i = 0; i < result.evidence_ledger.size(); ++i) {
        const auto& row = result.evidence_ledger[i];
        if (!ledger_row_is_promotable_event_evidence(row)) {
            continue;
        }
        promotable_indices.push_back(i);
    }
    if (promotable_indices.empty()) {
        return;
    }

    std::sort(promotable_indices.begin(), promotable_indices.end(), [&](size_t lhs, size_t rhs) {
        const auto& a = result.evidence_ledger[lhs];
        const auto& b = result.evidence_ledger[rhs];
        if (a.tid != b.tid) {
            return a.tid < b.tid;
        }
        const int32_t a_start = ledger_row_breakpoint_event_start(a);
        const int32_t b_start = ledger_row_breakpoint_event_start(b);
        if (a_start != b_start) {
            return a_start < b_start;
        }
        return lhs < rhs;
    });

    size_t cluster_start = 0;
    while (cluster_start < promotable_indices.size()) {
        size_t best = promotable_indices[cluster_start];
        std::vector<size_t> cluster_indices;
        cluster_indices.push_back(promotable_indices[cluster_start]);
        size_t cluster_end = cluster_start + 1;
        while (cluster_end < promotable_indices.size() &&
               ledger_rows_share_promotable_event_cluster(
                   result.evidence_ledger[promotable_indices[cluster_start]],
                   result.evidence_ledger[promotable_indices[cluster_end]])) {
            const size_t idx = promotable_indices[cluster_end];
            cluster_indices.push_back(idx);
            if (prefer_ledger_row_for_promotion(
                    result.evidence_ledger[idx],
                    result.evidence_ledger[best])) {
                best = idx;
            }
            ++cluster_end;
        }
        const bool requires_stable_community =
            ledger_row_requires_stable_event_community(result.evidence_ledger[best]);
        const bool has_stable_community =
            event_cluster_has_stable_low_af_support(
                best,
                result.evidence_ledger);
        if (!requires_stable_community || has_stable_community) {
            FinalCall call = promoted_call_from_ledger_row(result.evidence_ledger[best]);
            if (has_stable_community) {
                annotate_stable_event_community_certificate(call);
            }
            if (has_self_stable_low_af_interval_certificate(
                    result.evidence_ledger[best])) {
                annotate_stable_event_interval_certificate(call);
            }
            result.final_calls.push_back(std::move(call));
        }
        cluster_start = cluster_end;
    }
}

bool conformal_feature_dominates(
    const ConformalFeatureVector& lhs,
    const ConformalFeatureVector& rhs) {
    if (lhs.pro_te.empty() || lhs.pro_te.size() != rhs.pro_te.size()) {
        return false;
    }
    bool strictly_better = lhs.ref_span_reads < rhs.ref_span_reads;
    for (size_t i = 0; i < lhs.pro_te.size(); ++i) {
        if (lhs.pro_te[i] < rhs.pro_te[i]) {
            return false;
        }
        if (lhs.pro_te[i] > rhs.pro_te[i]) {
            strictly_better = true;
        }
    }
    if (lhs.ref_span_reads > rhs.ref_span_reads) {
        return false;
    }
    return strictly_better;
}

void remove_pareto_dominated_final_calls(PipelineResult& result) {
    if (result.final_calls.size() < 2) {
        return;
    }
    std::vector<ConformalFeatureVector> features;
    features.reserve(result.final_calls.size());
    for (size_t i = 0; i < result.final_calls.size(); ++i) {
        features.push_back(conformal_feature_from_call(result.final_calls[i], i));
    }
    std::vector<bool> dominated(result.final_calls.size(), false);
    for (size_t i = 0; i < features.size(); ++i) {
        for (size_t j = 0; j < features.size(); ++j) {
            if (i == j) {
                continue;
            }
            if (!same_call_locus(result.final_calls[i], result.final_calls[j])) {
                continue;
            }
            if (is_promoted_event_cluster_call(result.final_calls[j]) &&
                !is_promoted_event_cluster_call(result.final_calls[i])) {
                continue;
            }
            if (conformal_feature_dominates(features[j], features[i])) {
                dominated[i] = true;
                break;
            }
        }
    }
    size_t write = 0;
    for (size_t read = 0; read < result.final_calls.size(); ++read) {
        if (!dominated[read]) {
            if (write != read) {
                result.final_calls[write] = std::move(result.final_calls[read]);
            }
            ++write;
        }
    }
    result.final_calls.resize(write);
}

void resolve_overlapping_promoted_event_envelopes(PipelineResult& result) {
    if (result.final_calls.size() < 2) {
        return;
    }

    std::vector<bool> removed(result.final_calls.size(), false);
    for (size_t i = 0; i < result.final_calls.size(); ++i) {
        if (removed[i] || !is_promoted_event_cluster_call(result.final_calls[i])) {
            continue;
        }
        size_t best = i;
        for (size_t j = i + 1; j < result.final_calls.size(); ++j) {
            if (removed[j]) {
                continue;
            }
            if (result.final_calls[best].tid != result.final_calls[j].tid) {
                break;
            }
            if (!promoted_event_envelopes_should_compete(result.final_calls[best],
                                                         result.final_calls[j])) {
                continue;
            }
            if (prefer_promoted_event_envelope_representative(
                    result.final_calls[j],
                    result.final_calls[best])) {
                removed[best] = true;
                best = j;
            } else {
                removed[j] = true;
            }
        }
    }

    size_t write = 0;
    for (size_t read = 0; read < result.final_calls.size(); ++read) {
        if (!removed[read]) {
            if (write != read) {
                result.final_calls[write] = std::move(result.final_calls[read]);
            }
            ++write;
        }
    }
    result.final_calls.resize(write);
}

bool final_call_is_collapsible_event_fragment(const FinalCall& call) {
    return is_final_event_qc_for_conformal(call.final_qc);
}

bool final_call_is_structural_event_context_representative(
    const FinalCall& call) {
    return starts_with_text(call.final_qc, "PASS_STRUCTURAL_INSERTION") &&
           (is_promoted_event_cluster_call(call) ||
            call.final_qc.find("EVENT_FRAGMENT_COLLAPSED") != std::string::npos ||
            call.final_qc.find("EVENT_COMMUNITY") != std::string::npos);
}

bool final_call_is_te_event_context_representative(const FinalCall& call) {
    return starts_with_text(call.final_qc, "PASS_TE") &&
           (final_call_has_te_origin_certificate(call) ||
            is_promoted_event_cluster_call(call));
}

bool final_calls_are_mixed_te_structural_event_context(
    const FinalCall& a,
    const FinalCall& b) {
    if (a.event_consensus_len <= 0 || b.event_consensus_len <= 0) {
        return false;
    }
    const bool a_structural = final_call_is_structural_event_context_representative(a);
    const bool b_structural = final_call_is_structural_event_context_representative(b);
    const bool a_te = final_call_is_te_event_context_representative(a);
    const bool b_te = final_call_is_te_event_context_representative(b);
    if (!((a_structural && b_te) || (b_structural && a_te))) {
        return false;
    }
    return final_calls_have_overlapping_event_context(a, b);
}

bool final_calls_are_same_supported_event_fragment(
    const FinalCall& a,
    const FinalCall& b) {
    if (!final_call_is_collapsible_event_fragment(a) ||
        !final_call_is_collapsible_event_fragment(b)) {
        return false;
    }
    const bool structural_owner_context_pair =
        starts_with_text(a.final_qc, "PASS_STRUCTURAL_INSERTION") &&
        starts_with_text(b.final_qc, "PASS_STRUCTURAL_INSERTION") &&
        a.event_consensus_len > 0 &&
        b.event_consensus_len > 0 &&
        final_calls_have_overlapping_event_context(a, b);
    if (structural_owner_context_pair) {
        return true;
    }
    if (final_calls_are_mixed_te_structural_event_context(a, b)) {
        return true;
    }
    if (!final_calls_have_overlapping_breakpoint_intervals(a, b)) {
        if (final_call_event_envelopes_overlap(a, b) &&
            final_calls_have_overlapping_event_context(a, b) &&
            ((starts_with_text(a.final_qc, "PASS_STRUCTURAL_INSERTION") &&
              starts_with_text(b.final_qc, "PASS_STRUCTURAL_INSERTION")) ||
             final_calls_share_support_containment_edge(a, b))) {
            return true;
        }
        return final_calls_share_identical_support(a, b);
    }
    return final_calls_share_support_edge(a, b);
}

bool prefer_fragment_component_representative(
    const FinalCall& candidate,
    const FinalCall& incumbent) {
    if (final_calls_are_mixed_te_structural_event_context(candidate, incumbent)) {
        return prefer_promoted_event_envelope_representative(candidate, incumbent);
    }
    if (promoted_event_envelopes_should_compete(candidate, incumbent)) {
        return prefer_promoted_event_envelope_representative(candidate, incumbent);
    }
    if (candidate.support_qnames.size() != incumbent.support_qnames.size()) {
        return candidate.support_qnames.size() > incumbent.support_qnames.size();
    }
    return prefer_new_call(candidate, incumbent);
}

FinalCall collapsed_fragment_component_call(
    const std::vector<size_t>& component,
    const std::vector<FinalCall>& calls) {
    size_t best = component.front();
    int32_t merged_left = std::numeric_limits<int32_t>::max();
    int32_t merged_right = -1;
    for (const size_t idx : component) {
        const auto [call_left, call_right] = final_call_event_envelope(calls[idx]);
        if (call_left >= 0) {
            merged_left = std::min(merged_left, call_left);
        }
        if (call_right >= 0) {
            merged_right = std::max(merged_right, call_right);
        }
        if (prefer_fragment_component_representative(calls[idx], calls[best])) {
            best = idx;
        }
    }

    FinalCall merged = calls[best];
    if (merged_left != std::numeric_limits<int32_t>::max()) {
        merged.bp_left = merged_left;
    }
    if (merged_right >= 0) {
        merged.bp_right = merged_right;
    }
    if (merged.bp_left > merged.bp_right) {
        std::swap(merged.bp_left, merged.bp_right);
    }
    std::vector<std::string> support;
    for (const size_t idx : component) {
        support.insert(
            support.end(),
            calls[idx].support_qnames.begin(),
            calls[idx].support_qnames.end());
        merged.raw_cigar_insert_reads = std::max(
            merged.raw_cigar_insert_reads,
            calls[idx].raw_cigar_insert_reads);
        merged.max_raw_cigar_insert_len = std::max(
            merged.max_raw_cigar_insert_len,
            calls[idx].max_raw_cigar_insert_len);
        merged.ref_span_reads = std::max(
            merged.ref_span_reads,
            calls[idx].ref_span_reads);
        merged.low_mapq_ref_span_reads = std::max(
            merged.low_mapq_ref_span_reads,
            calls[idx].low_mapq_ref_span_reads);
    }
    merged.support_qnames = normalized_support_qnames(support);
    const int32_t unique_support =
        static_cast<int32_t>(std::min<size_t>(
            merged.support_qnames.size(),
            static_cast<size_t>(std::numeric_limits<int32_t>::max())));
    merged.support_reads = std::max(merged.support_reads, unique_support);
    merged.alt_struct_reads = std::max(merged.alt_struct_reads, unique_support);
    if (merged.final_qc.find("EVENT_FRAGMENT_COLLAPSED") == std::string::npos) {
        merged.final_qc += "|EVENT_FRAGMENT_COLLAPSED";
    }
    if (merged.mechanistic_blocks.empty() || merged.mechanistic_blocks == "NA") {
        merged.mechanistic_blocks = "event;read_overlap";
    } else if (merged.mechanistic_blocks.find("read_overlap") == std::string::npos) {
        merged.mechanistic_blocks += ";read_overlap";
    }
    return merged;
}

void collapse_shared_support_final_call_fragments(PipelineResult& result) {
    if (result.final_calls.size() < 2) {
        return;
    }

    std::vector<size_t> parent(result.final_calls.size());
    for (size_t i = 0; i < parent.size(); ++i) {
        parent[i] = i;
    }
    const auto find_root = [&parent](size_t value) {
        size_t root = value;
        while (parent[root] != root) {
            root = parent[root];
        }
        while (parent[value] != value) {
            const size_t next = parent[value];
            parent[value] = root;
            value = next;
        }
        return root;
    };
    auto unite = [&parent, &find_root](size_t lhs, size_t rhs) {
        const size_t lhs_root = find_root(lhs);
        const size_t rhs_root = find_root(rhs);
        if (lhs_root != rhs_root) {
            parent[rhs_root] = lhs_root;
        }
    };

    for (size_t i = 0; i < result.final_calls.size(); ++i) {
        for (size_t j = i + 1; j < result.final_calls.size(); ++j) {
            if (result.final_calls[i].tid != result.final_calls[j].tid) {
                break;
            }
            if (final_calls_are_same_supported_event_fragment(
                    result.final_calls[i],
                    result.final_calls[j])) {
                unite(i, j);
            }
        }
    }

    std::unordered_map<size_t, std::vector<size_t>> components;
    components.reserve(result.final_calls.size());
    for (size_t i = 0; i < result.final_calls.size(); ++i) {
        components[find_root(i)].push_back(i);
    }

    std::vector<FinalCall> collapsed;
    collapsed.reserve(result.final_calls.size());
    for (size_t i = 0; i < result.final_calls.size(); ++i) {
        const size_t root = find_root(i);
        auto found = components.find(root);
        if (found == components.end()) {
            continue;
        }
        if (found->second.size() == 1) {
            collapsed.push_back(result.final_calls[i]);
        } else {
            collapsed.push_back(
                collapsed_fragment_component_call(found->second, result.final_calls));
        }
        components.erase(found);
    }

    result.final_calls = std::move(collapsed);
}

double event_mechanistic_e_value(const FinalCall& call) {
    if (call.mechanistic_blocks.empty() || call.mechanistic_blocks == "NA") {
        return 0.0;
    }
    if (call.mechanistic_ref_conflict_signal > 0.50) {
        return 0.0;
    }
    if (call.ref_span_reads > std::max(2, call.alt_struct_reads / 2)) {
        return 0.0;
    }
    const double lower_log_bf = std::min(
        call.mechanistic_lower_log_bf_te_vs_artifact,
        call.mechanistic_lower_log_bf_te_vs_non_te);
    const double log_e = lower_log_bf - std::max(0.0, call.mechanistic_ambiguity_width);
    if (!std::isfinite(log_e) || log_e <= 0.0) {
        return 0.0;
    }
    return std::exp(std::min(log_e, 40.0));
}

double low_allele_fraction_event_e_value(const FinalCall& call) {
    if (!(has_low_allele_fraction_event_certificate(call) ||
          has_aggregated_low_allele_fraction_event_certificate(call) ||
          has_self_stable_low_af_interval_certificate(call)) ||
        (call.final_qc.find("EVENT_COMMUNITY_STABLE") == std::string::npos &&
         call.final_qc.find("EVENT_INTERVAL_STABLE") == std::string::npos)) {
        return 0.0;
    }
    const int32_t alt = std::max(0, call.alt_struct_reads);
    const int32_t ref = std::max(0, call.ref_span_reads);
    const double log_bf = low_allele_fraction_event_log_bf(alt, ref);
    if (!std::isfinite(log_bf) || log_bf < 6.0) {
        return 0.0;
    }
    return std::exp(std::min(log_bf, 40.0));
}

double ref_unopposed_bilateral_partial_anchor_e_value(const FinalCall& call) {
    if (call.final_qc.find("EVENT_BILATERAL_PARTIAL_ANCHOR") == std::string::npos) {
        return 0.0;
    }
    if (!starts_with_text(call.final_qc, "PASS_STRUCTURAL_INSERTION") ||
        call.mechanistic_blocks.empty() || call.mechanistic_blocks == "NA") {
        return 0.0;
    }
    const int32_t alt = std::max(0, call.alt_struct_reads);
    const int32_t ref = std::max(0, call.ref_span_reads);
    if (alt < 2 || ref != 0 || call.mechanistic_ref_conflict_signal > 1e-9) {
        return 0.0;
    }
    const double log_bf = low_allele_fraction_event_log_bf(alt, ref);
    if (!std::isfinite(log_bf) || log_bf <= 0.0) {
        return 0.0;
    }
    return std::exp(std::min(log_bf, 40.0));
}

void apply_event_ebh_fallback(PipelineResult& result, double target_fdr) {
    if (result.final_calls.empty()) {
        return;
    }

    std::vector<std::pair<double, size_t>> evalues;
    evalues.reserve(result.final_calls.size());
    for (size_t i = 0; i < result.final_calls.size(); ++i) {
        const auto& call = result.final_calls[i];
        if (call.conformal_qc == "PASS_CONFORMAL_FDR" ||
            call.conformal_qc == "PASS_EVENT_LFDR_FALLBACK") {
            continue;
        }
        const double evalue = std::max({
            event_mechanistic_e_value(call),
            low_allele_fraction_event_e_value(call),
            ref_unopposed_bilateral_partial_anchor_e_value(call),
        });
        if (evalue > 0.0 && std::isfinite(evalue)) {
            evalues.emplace_back(evalue, i);
        }
    }
    if (evalues.empty()) {
        return;
    }

    std::sort(evalues.begin(), evalues.end(), [](const auto& lhs, const auto& rhs) {
        if (lhs.first != rhs.first) {
            return lhs.first > rhs.first;
        }
        return lhs.second < rhs.second;
    });

    const double q = std::clamp(target_fdr, 0.0, 1.0);
    size_t selected_prefix = 0;
    const double m = static_cast<double>(evalues.size());
    for (size_t rank = 1; rank <= evalues.size(); ++rank) {
        const double threshold = m / (q * static_cast<double>(rank));
        if (q > 0.0 && evalues[rank - 1].first >= threshold) {
            selected_prefix = rank;
        }
    }
    for (size_t rank = 0; rank < selected_prefix; ++rank) {
        auto& call = result.final_calls[evalues[rank].second];
        call.conformal_qc = "PASS_EVENT_EBH";
        call.conformal_by_threshold = evalues[rank].first;
    }
}

void apply_event_bayesian_fdr_fallback(PipelineResult& result, double target_fdr) {
    if (result.final_calls.empty()) {
        return;
    }

    std::vector<std::pair<double, size_t>> local_fdrs;
    local_fdrs.reserve(result.final_calls.size());
    for (size_t i = 0; i < result.final_calls.size(); ++i) {
        const auto& call = result.final_calls[i];
        if (call.conformal_qc == "PASS_CONFORMAL_FDR" ||
            call.conformal_qc == "PASS_EVENT_LFDR_FALLBACK" ||
            call.conformal_qc == "PASS_EVENT_EBH") {
            continue;
        }
        const bool calibrated =
            call.lfdr_qc == "PASS_TE_LFDR" ||
            call.posterior_qc == "PASS_TE_POSTERIOR" ||
            call.robust_mechanistic_qc == "PASS_TE_LFDR" ||
            call.worst_case_lfdr <= 0.10 ||
            call.robust_mechanistic_worst_case_lfdr <= 0.10 ||
            call.te_posterior >= 0.95;
        if (!calibrated) {
            continue;
        }
        double local_fdr = 1.0;
        if (call.lfdr_qc == "PASS_TE_LFDR" || call.worst_case_lfdr <= 0.10) {
            local_fdr = std::min(local_fdr, std::clamp(call.worst_case_lfdr, 0.0, 1.0));
        }
        if (call.robust_mechanistic_qc == "PASS_TE_LFDR" ||
            call.robust_mechanistic_worst_case_lfdr <= 0.10) {
            local_fdr = std::min(
                local_fdr,
                std::clamp(call.robust_mechanistic_worst_case_lfdr, 0.0, 1.0));
        }
        if (call.posterior_qc == "PASS_TE_POSTERIOR") {
            local_fdr = std::min(local_fdr, 0.10);
        } else if (call.te_posterior >= 0.95) {
            local_fdr = std::min(local_fdr, 1.0 - std::clamp(call.te_posterior, 0.0, 1.0));
        }
        local_fdrs.emplace_back(local_fdr, i);
    }
    if (local_fdrs.empty()) {
        return;
    }

    std::sort(local_fdrs.begin(), local_fdrs.end(), [](const auto& lhs, const auto& rhs) {
        if (lhs.first != rhs.first) {
            return lhs.first < rhs.first;
        }
        return lhs.second < rhs.second;
    });

    const double q = std::clamp(target_fdr, 0.0, 1.0);
    double cumulative = 0.0;
    size_t selected_prefix = 0;
    for (size_t rank = 1; rank <= local_fdrs.size(); ++rank) {
        cumulative += local_fdrs[rank - 1].first;
        const double mean_lfdr = cumulative / static_cast<double>(rank);
        if (mean_lfdr <= q) {
            selected_prefix = rank;
        }
    }
    for (size_t rank = 0; rank < selected_prefix; ++rank) {
        auto& call = result.final_calls[local_fdrs[rank].second];
        call.conformal_qc = "PASS_EVENT_BFDR";
        call.conformal_by_threshold = local_fdrs[rank].first;
    }
}

bool is_structural_insertion_call(const FinalCall& call) {
    return starts_with_text(call.final_qc, "PASS_STRUCTURAL_INSERTION");
}

double structural_event_local_fdr(const FinalCall& call) {
    if (!is_structural_insertion_call(call)) {
        return 1.0;
    }
    const int32_t alt = std::max(0, call.alt_struct_reads);
    const int32_t ref = std::max(0, call.ref_span_reads);
    if (alt <= 0) {
        return 1.0;
    }
    if (call.mechanistic_blocks.empty() || call.mechanistic_blocks == "NA") {
        return 1.0;
    }
    const double heterozygous_balance =
        heterozygous_balance_log_evidence(alt, ref);
    const double raw_insert_log_bf = low_allele_fraction_event_log_bf(
        std::max(0, call.raw_cigar_insert_reads),
        ref);
    const double direct_insert_evidence = std::isfinite(raw_insert_log_bf)
        ? std::max(0.0, raw_insert_log_bf)
        : 0.0;
    const double assembled_insert_evidence =
        0.05 * std::log1p(static_cast<double>(std::max(0, call.max_raw_cigar_insert_len))) +
        0.03 * std::log1p(static_cast<double>(std::max(0, call.event_consensus_len)));
    const bool has_structure_model =
        call.mechanistic_blocks.find("structure") != std::string::npos;
    const double event_structure_log_bf = has_structure_model
        ? std::max(call.te_structure_log_evidence, call.nonte_structure_log_evidence) -
              call.artifact_structure_log_evidence
        : -std::numeric_limits<double>::infinity();
    const double event_structure_evidence = std::isfinite(event_structure_log_bf)
        ? event_structure_log_bf
        : 0.0;
    const bool structure_supports_event = event_structure_evidence > 0.0;
    const bool balanced_heterozygous =
        heterozygous_balance >= 1.75 &&
        alt >= 8 &&
        call.mechanistic_ref_conflict_signal <= 0.65;
    const bool direct_raw_insertion = direct_insert_evidence > 0.0;
    if (call.mechanistic_ref_conflict_signal > 0.50 &&
        !balanced_heterozygous &&
        !direct_raw_insertion &&
        !structure_supports_event) {
        return 1.0;
    }

    const double artifact_posterior =
        std::clamp(call.artifact_posterior, 0.0, 1.0);
    const double non_te_posterior =
        std::clamp(call.non_te_posterior, 0.0, 1.0);
    const double event_posterior =
        std::clamp(non_te_posterior + std::clamp(call.te_posterior, 0.0, 1.0), 0.0, 1.0);
    const double support_log_evidence =
        std::log(static_cast<double>(alt + 1)) -
        std::log(static_cast<double>(ref + 1));
    const double ref_conflict_penalty = balanced_heterozygous
        ? 0.25 * std::clamp(call.mechanistic_ref_conflict_signal, 0.0, 1.0)
        : 1.20 * std::clamp(call.mechanistic_ref_conflict_signal, 0.0, 1.0);
    const double log_event_odds =
        std::log(std::max(event_posterior, 1e-6)) -
        std::log(std::max(artifact_posterior, 1e-6)) +
        (0.35 * support_log_evidence) +
        event_structure_evidence +
        direct_insert_evidence +
        assembled_insert_evidence +
        (0.80 * heterozygous_balance) -
        ref_conflict_penalty;
    const double event_odds = std::exp(std::clamp(log_event_odds, -60.0, 60.0));
    return 1.0 / (1.0 + event_odds);
}

void apply_structural_event_bayesian_fdr_fallback(PipelineResult& result, double target_fdr) {
    if (result.final_calls.empty()) {
        return;
    }

    std::vector<std::pair<double, size_t>> local_fdrs;
    local_fdrs.reserve(result.final_calls.size());
    for (size_t i = 0; i < result.final_calls.size(); ++i) {
        const auto& call = result.final_calls[i];
        if (call.conformal_qc == "PASS_CONFORMAL_FDR" ||
            call.conformal_qc == "PASS_EVENT_LFDR_FALLBACK" ||
            call.conformal_qc == "PASS_EVENT_EBH" ||
            call.conformal_qc == "PASS_EVENT_BFDR") {
            continue;
        }
        const double local_fdr = structural_event_local_fdr(call);
        if (local_fdr <= 0.25) {
            local_fdrs.emplace_back(local_fdr, i);
        }
    }
    if (local_fdrs.empty()) {
        return;
    }

    std::sort(local_fdrs.begin(), local_fdrs.end(), [](const auto& lhs, const auto& rhs) {
        if (lhs.first != rhs.first) {
            return lhs.first < rhs.first;
        }
        return lhs.second < rhs.second;
    });

    const double q = std::clamp(target_fdr, 0.0, 1.0);
    double cumulative = 0.0;
    size_t selected_prefix = 0;
    for (size_t rank = 1; rank <= local_fdrs.size(); ++rank) {
        cumulative += local_fdrs[rank - 1].first;
        const double mean_lfdr = cumulative / static_cast<double>(rank);
        if (mean_lfdr <= q) {
            selected_prefix = rank;
        }
    }
    for (size_t rank = 0; rank < selected_prefix; ++rank) {
        auto& call = result.final_calls[local_fdrs[rank].second];
        call.conformal_qc = "PASS_EVENT_EXISTENCE_BFDR";
        call.conformal_by_threshold = local_fdrs[rank].first;
    }
}

void apply_sample_local_conformal_selector(PipelineResult& result, double target_fdr) {
    if (result.final_calls.empty() || result.evidence_ledger.empty()) {
        return;
    }

    ConformalNullSelector selector;
    for (size_t i = 0; i < result.evidence_ledger.size(); ++i) {
        const auto& row = result.evidence_ledger[i];
        if (is_sample_local_null_control(row)) {
            selector.add_null_control(conformal_feature_from_row(row, i));
        }
    }

    std::vector<ConformalFeatureVector> candidates;
    candidates.reserve(result.final_calls.size());
    for (size_t i = 0; i < result.final_calls.size(); ++i) {
        candidates.push_back(conformal_feature_from_call(result.final_calls[i], i));
    }

    const auto decisions = selector.select(candidates, target_fdr);
    std::unordered_map<std::string, ConformalSelectionResult> decision_by_id;
    decision_by_id.reserve(decisions.size());
    for (const auto& decision : decisions) {
        decision_by_id.emplace(decision.id, decision);
    }

    for (size_t i = 0; i < result.final_calls.size(); ++i) {
        auto& call = result.final_calls[i];
        const auto feature = conformal_feature_from_call(call, i);
        const auto found = decision_by_id.find(feature.id);
        if (found == decision_by_id.end()) {
            continue;
        }
        const auto& decision = found->second;
        call.conformal_null_p = decision.conformal_p;
        call.conformal_by_threshold = decision.by_threshold;
        call.conformal_dominated_nulls =
            static_cast<int32_t>(std::min<size_t>(
                decision.dominated_null_count,
                static_cast<size_t>(std::numeric_limits<int32_t>::max())));
        call.conformal_null_count =
            static_cast<int32_t>(std::min<size_t>(
                decision.null_count,
                static_cast<size_t>(std::numeric_limits<int32_t>::max())));
        call.conformal_qc = decision.qc;
    }

    const auto passes_event_lfdr_fallback = [target_fdr](const FinalCall& call) {
        const double q = std::clamp(target_fdr, 0.0, 1.0);
        if (call.mechanistic_blocks.empty() || call.mechanistic_blocks == "NA") {
            return false;
        }
        if (call.mechanistic_ref_conflict_signal > 0.50) {
            return false;
        }
        if (call.ref_span_reads > std::max(2, call.alt_struct_reads / 2)) {
            return false;
        }
        const bool robust_pass =
            call.robust_mechanistic_qc == "PASS_TE_LFDR" &&
            call.robust_mechanistic_worst_case_lfdr <= std::max(q, 0.10);
        return robust_pass;
    };
    const bool conformal_null_too_small = selector.null_count() < 20;
    if (conformal_null_too_small) {
        for (auto& call : result.final_calls) {
            if (call.conformal_qc != "PASS_CONFORMAL_FDR" &&
                passes_event_lfdr_fallback(call)) {
                call.conformal_qc = "PASS_EVENT_LFDR_FALLBACK";
                call.conformal_null_count =
                    static_cast<int32_t>(std::min<size_t>(
                        selector.null_count(),
                        static_cast<size_t>(std::numeric_limits<int32_t>::max())));
            }
        }
    }

    apply_event_ebh_fallback(result, target_fdr);
    apply_event_bayesian_fdr_fallback(result, target_fdr);
    apply_structural_event_bayesian_fdr_fallback(result, target_fdr);

    result.final_calls.erase(
        std::remove_if(
            result.final_calls.begin(),
            result.final_calls.end(),
            [](const FinalCall& call) {
                return !final_call_has_reportable_conformal_certificate(call);
            }),
        result.final_calls.end());
    remove_pareto_dominated_final_calls(result);
}

void apply_final_report_mode(PipelineResult& result, FinalReportMode report_mode) {
    if (report_mode == FinalReportMode::Legacy) {
        return;
    }
    result.final_calls.erase(
        std::remove_if(
            result.final_calls.begin(),
            result.final_calls.end(),
            [](const FinalCall& call) {
                return !final_call_passes_te_calibrated_report_gate(call);
            }),
        result.final_calls.end());
}

}  // namespace

void finalize_final_calls(PipelineResult& result) {
    finalize_final_calls(result, kDefaultFinalConformalFdr);
}

void finalize_final_calls(PipelineResult& result, double target_fdr) {
    FinalCallFilterConfig filter_config;
    finalize_final_calls(result, target_fdr, filter_config);
}

// L3 (per-sample step): estimate the beta-binomial overdispersion from the run's
// own alt/ref count distribution (all candidate events) and recompute the final
// calls' reported genotype quality with it, so GQ/AF reflect the sample's actual
// count overdispersion rather than a fixed default. Deterministic: runs once over
// the fully-populated evidence ledger in the single-threaded finalization pass.
void apply_sample_overdispersion_calibration(PipelineResult& result) {
    std::vector<AltDepthObservation> observations;
    observations.reserve(result.evidence_ledger.size());
    for (const EvidenceLedgerRow& row : result.evidence_ledger) {
        const int32_t alt = std::max(0, row.alt_struct_reads);
        const int32_t depth = alt + std::max(0, row.ref_span_reads);
        observations.push_back({alt, depth});
    }
    const double rho = estimate_alt_depth_overdispersion(
        observations, result.estimated_overdispersion);
    result.estimated_overdispersion = rho;

    for (FinalCall& call : result.final_calls) {
        const int32_t depth =
            std::max(0, call.alt_struct_reads) + std::max(0, call.ref_span_reads);
        if (depth <= 0) {
            continue;
        }
        EventGenotypeInput input;
        input.alt_struct_reads = call.alt_struct_reads;
        input.ref_span_reads = call.ref_span_reads;
        input.overdispersion = rho;
        const EventGenotypeDecision decision = genotype_event_from_alt_vs_ref(input);
        call.gq = decision.gq;
        call.af = decision.allele_fraction;
    }
}

void apply_sequence_family_commitments(PipelineResult& result) {
    for (FinalCall& call : result.final_calls) {
        if (!call.sequence_family_commit_eligible ||
            call.family_committed ||
            !starts_with_text(call.final_qc, "PASS_TE")) {
            continue;
        }
        call.family = call.sequence_family_candidate;
        call.subfamily = call.sequence_subfamily_candidate;
        call.te_name = call.subfamily;
        call.family_committed = true;
        append_qc_token(call.final_qc, "SEQUENCE_FAMILY_COMMITTED");
    }
}

void finalize_final_calls(
    PipelineResult& result,
    double target_fdr,
    const FinalCallFilterConfig& filter_config) {
    append_aggregated_event_community_calls(result);
    append_promoted_event_cluster_calls(result);
    apply_sample_overdispersion_calibration(result);
    std::sort(result.final_calls.begin(), result.final_calls.end(), final_call_sort_less);

    std::vector<FinalCall> deduped;
    deduped.reserve(result.final_calls.size());

    size_t cluster_start = 0;
    while (cluster_start < result.final_calls.size()) {
        const FinalCall* best = &result.final_calls[cluster_start];
        size_t cluster_end = cluster_start + 1;
        while (cluster_end < result.final_calls.size() &&
               same_call_locus(result.final_calls[cluster_start], result.final_calls[cluster_end])) {
            if (prefer_new_call(result.final_calls[cluster_end], *best)) {
                best = &result.final_calls[cluster_end];
            }
            ++cluster_end;
        }
        deduped.push_back(*best);
        cluster_start = cluster_end;
    }

    result.final_calls = std::move(deduped);
    std::sort(result.final_calls.begin(), result.final_calls.end(), final_call_sort_less);
    resolve_overlapping_promoted_event_envelopes(result);
    collapse_shared_support_final_call_fragments(result);
    std::sort(result.final_calls.begin(), result.final_calls.end(), final_call_sort_less);
    filter_short_final_insertions(result, filter_config);
    apply_sample_local_conformal_selector(result, target_fdr);
    apply_final_report_mode(result, filter_config.report_mode);
    apply_sequence_family_commitments(result);

    result.final_pass_calls = static_cast<int64_t>(result.final_calls.size());
}

void finalize_final_calls(PipelineResult& result, const TSDDetector& reference) {
    finalize_final_calls(result, reference, kDefaultFinalConformalFdr);
}

void finalize_final_calls(PipelineResult& result, const TSDDetector& reference, double target_fdr) {
    FinalCallFilterConfig filter_config;
    finalize_final_calls(result, reference, target_fdr, filter_config);
}

void finalize_final_calls(
    PipelineResult& result,
    const TSDDetector& reference,
    double target_fdr,
    const FinalCallFilterConfig& filter_config) {
    if (!reference.can_fetch_reference() || result.final_calls.empty()) {
        finalize_final_calls(result, target_fdr, filter_config);
        return;
    }

    constexpr int32_t kFinalCallReferenceNContextBp = 250;
    const auto call_overlaps_n_rich_reference_context =
        [&reference](const FinalCall& call) {
            if (call.chrom.empty()) {
                return false;
            }
            const int32_t left = std::min(call.bp_left, call.bp_right);
            const int32_t right = std::max(call.bp_left, call.bp_right);
            const int32_t anchor = call.pos >= 0
                ? call.pos
                : ((left >= 0 && right >= 0) ? left + ((right - left) / 2) : -1);
            const int32_t interval_start = (left >= 0)
                ? std::max(0, left - kFinalCallReferenceNContextBp)
                : ((anchor >= 0) ? std::max(0, anchor - kFinalCallReferenceNContextBp) : -1);
            const int32_t interval_end = (right >= 0)
                ? right + kFinalCallReferenceNContextBp + 1
                : ((anchor >= 0) ? anchor + kFinalCallReferenceNContextBp + 1 : -1);
            if (interval_start >= 0 && interval_end > interval_start &&
                reference.reference_interval_is_n_rich(
                    call.chrom,
                    interval_start,
                    interval_end)) {
                return true;
            }
            return (anchor >= 0 && reference.reference_position_is_poly_n(call.chrom, anchor)) ||
                (call.bp_left >= 0 && reference.reference_position_is_poly_n(call.chrom, call.bp_left)) ||
                (call.bp_right >= 0 && reference.reference_position_is_poly_n(call.chrom, call.bp_right));
        };

    result.final_calls.erase(
        std::remove_if(
            result.final_calls.begin(),
            result.final_calls.end(),
            call_overlaps_n_rich_reference_context),
        result.final_calls.end());
    finalize_final_calls(result, target_fdr, filter_config);
}
