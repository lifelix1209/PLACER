#include "te_sequence_explainer.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <sstream>

namespace placer {
namespace {

std::string normalized_label(const std::string& value, const std::string& fallback) {
    if (value.empty() || value == "NA") {
        return fallback;
    }
    return value;
}

int32_t rounded_bases(int32_t insert_len, double fraction) {
    if (insert_len <= 0 || !std::isfinite(fraction)) {
        return 0;
    }
    const double clamped = std::clamp(fraction, 0.0, 1.0);
    return static_cast<int32_t>(std::lround(static_cast<double>(insert_len) * clamped));
}

double clamp01(double value) {
    return std::clamp(value, 0.0, 1.0);
}

double logistic(double value) {
    if (value >= 40.0) {
        return 1.0;
    }
    if (value <= -40.0) {
        return 0.0;
    }
    return 1.0 / (1.0 + std::exp(-value));
}

int32_t terminal_poly_at_run(const std::string& seq) {
    if (seq.empty()) {
        return 0;
    }
    const char base = seq.back();
    if (base != 'A' && base != 'T') {
        return 0;
    }
    int32_t run = 0;
    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        if (*it != base) {
            break;
        }
        ++run;
    }
    return run;
}

double interval_entropy_norm(
    const std::string& seq,
    int32_t start,
    int32_t end) {
    start = std::max(0, std::min(start, static_cast<int32_t>(seq.size())));
    end = std::max(start, std::min(end, static_cast<int32_t>(seq.size())));
    const int32_t len = end - start;
    if (len <= 0) {
        return 0.0;
    }
    std::array<int32_t, 4> counts{0, 0, 0, 0};
    for (int32_t i = start; i < end; ++i) {
        switch (seq[static_cast<size_t>(i)]) {
            case 'A': counts[0] += 1; break;
            case 'C': counts[1] += 1; break;
            case 'G': counts[2] += 1; break;
            case 'T': counts[3] += 1; break;
            default: break;
        }
    }
    double entropy = 0.0;
    for (const int32_t count : counts) {
        if (count <= 0) {
            continue;
        }
        const double p = static_cast<double>(count) / static_cast<double>(len);
        entropy -= p * std::log(p);
    }
    return clamp01(entropy / std::log(4.0));
}

TeAnnotationStatus status_from_qc(
    const std::string& qc_reason,
    const std::string& family,
    const std::string& subfamily) {
    if (qc_reason == "TE_LIBRARY_UNAVAILABLE") {
        return TeAnnotationStatus::kUnavailable;
    }
    if (qc_reason == "EMPTY_INSERT_SEQUENCE") {
        return TeAnnotationStatus::kEmpty;
    }
    if (qc_reason == "INSERT_SEQ_TOO_SHORT") {
        return TeAnnotationStatus::kTooShort;
    }
    if (qc_reason == "NO_TE_ALIGNMENT_SHORTLIST" ||
        qc_reason == "NO_TE_ALIGNMENT_MATCH" ||
        qc_reason == "NO_TE_ALIGNMENT") {
        return TeAnnotationStatus::kNoCandidate;
    }
    if (qc_reason == "PASS_INSERT_TE_ALIGNMENT") {
        return TeAnnotationStatus::kResolved;
    }
    if (qc_reason == "PASS_INSERT_TE_ALIGNMENT_FAMILY_ONLY") {
        return TeAnnotationStatus::kFamilyOnly;
    }
    if (qc_reason == "PASS_INSERT_TE_ALIGNMENT_UNKNOWN") {
        return TeAnnotationStatus::kUnknownTe;
    }
    if (qc_reason == "TE_ALIGNMENT_LOW_IDENTITY" ||
        qc_reason == "TE_ALIGNMENT_LOW_QUERY_COVERAGE") {
        return TeAnnotationStatus::kNonTeLike;
    }
    if (family.empty() || family == "UNKNOWN" ||
        subfamily.empty() || subfamily == "UNKNOWN") {
        return TeAnnotationStatus::kAmbiguous;
    }
    return TeAnnotationStatus::kNoCandidate;
}

bool has_te_path(TeAnnotationStatus status) {
    return status == TeAnnotationStatus::kResolved ||
           status == TeAnnotationStatus::kFamilyOnly ||
           status == TeAnnotationStatus::kUnknownTe ||
           status == TeAnnotationStatus::kAmbiguous;
}

std::string te_path_label(
    TeAnnotationStatus status,
    const std::string& family,
    const std::string& subfamily) {
    if (status == TeAnnotationStatus::kUnknownTe) {
        return "UNKNOWN/UNKNOWN";
    }
    const std::string normalized_family = normalized_label(family, "UNKNOWN");
    const std::string normalized_subfamily =
        status == TeAnnotationStatus::kFamilyOnly
            ? "UNKNOWN"
            : normalized_label(subfamily, "UNKNOWN");
    return normalized_family + "/" + normalized_subfamily;
}

}  // namespace

TeSequenceExplanation explain_te_sequence_structure(
    const std::string& insert_seq,
    const std::string& qc_reason,
    const std::string& family,
    const std::string& subfamily,
    double best_identity,
    double effective_query_coverage,
    double annotation_residual_fraction,
    double annotation_masked_fraction,
    double cross_family_margin,
    double second_score,
    const std::string& sequence_model_label,
    double sequence_model_score) {
    const int32_t insert_len = static_cast<int32_t>(insert_seq.size());
    TeSequenceExplanation explanation;
    explanation.status = status_from_qc(qc_reason, family, subfamily);
    explanation.family = explanation.status == TeAnnotationStatus::kUnknownTe
        ? "UNKNOWN"
        : normalized_label(family, "UNKNOWN");
    explanation.subfamily =
        explanation.status == TeAnnotationStatus::kFamilyOnly ||
            explanation.status == TeAnnotationStatus::kUnknownTe
            ? "UNKNOWN"
            : normalized_label(subfamily, "UNKNOWN");

    const double high_complexity_residual =
        std::max(0.0, annotation_residual_fraction - annotation_masked_fraction);
    explanation.residual.unexplained_high_complexity_bases =
        rounded_bases(insert_len, high_complexity_residual);
    explanation.residual.low_complexity_only_bases =
        rounded_bases(insert_len, annotation_masked_fraction);
    explanation.residual.edit_distance =
        rounded_bases(insert_len, std::max(0.0, 1.0 - best_identity) *
                                   std::max(0.0, effective_query_coverage));
    explanation.residual.family_conflicts =
        (second_score > 0.0 && explanation.family == "UNKNOWN") ? 1 : 0;
    explanation.residual.subfamily_conflicts =
        (explanation.status == TeAnnotationStatus::kFamilyOnly ||
         explanation.status == TeAnnotationStatus::kUnknownTe ||
         explanation.subfamily == "UNKNOWN") ? 1 : 0;
    explanation.residual.segment_breaks =
        explanation.residual.unexplained_high_complexity_bases > 0 ? 1 : 0;
    explanation.unexplained_high_complexity_bp =
        explanation.residual.unexplained_high_complexity_bases;
    explanation.te_core_coverage = clamp01(effective_query_coverage);

    const int32_t te_core_bases = insert_len > 0 && has_te_path(explanation.status)
        ? std::max(1, std::min(insert_len, rounded_bases(insert_len, effective_query_coverage)))
        : 0;
    const int32_t residual_start = te_core_bases;
    const int32_t residual_len = std::max(0, insert_len - residual_start);
    const int32_t poly_bases = std::min(terminal_poly_at_run(insert_seq), residual_len);
    const int32_t residual_before_poly = std::max(0, residual_len - poly_bases);
    const double residual_entropy = interval_entropy_norm(
        insert_seq,
        residual_start,
        residual_start + residual_before_poly);

    // L2: likelihood-ratio structural decode of the insert tail -- a minimal,
    // duration-aware HSMM over the left-to-right state chain
    // TE_CORE -> [TRANSDUCTION] -> [POLYA] -> residual. Each optional state is
    // "on" when its summed per-base emission log-LR versus a generic-residual
    // null, plus a log prior-odds to open the state (its duration prior), is
    // positive; the state posterior is the logistic of that margin. This replaces
    // the heuristic feature-logistic posteriors and the hard 0.50 inclusion
    // cutoffs with an explicit emission-vs-null decode.
    //
    // POLYA emission: per-base homopolymer-purity LR ~ log(0.90 / 0.25) for an
    // A/T-dominant terminal run vs random sequence; the open-odds enforces a
    // minimum plausible run length before the state is opened.
    constexpr double kPolyPerBaseLogLr = 1.28;
    constexpr double kPolyOpenLogOdds = -3.0;
    const double poly_state_score =
        (static_cast<double>(poly_bases) * kPolyPerBaseLogLr) + kPolyOpenLogOdds;
    explanation.polyA_posterior = logistic(poly_state_score);

    // TRANSDUCTION emission: a high-complexity genomic 3' transduction between the
    // TE core and the polyA. Per-base emission rewards residual entropy above a
    // random-tail baseline; the open-odds is the transduction prior. A resolved
    // TE core (high identity) makes an adjacent transduction more credible.
    constexpr double kTransBaselineEntropy = 0.55;
    constexpr double kTransOpenLogOdds = -1.6;
    const double transduction_state_score =
        (static_cast<double>(residual_before_poly) * 0.06 *
         (residual_entropy - kTransBaselineEntropy)) +
        kTransOpenLogOdds +
        (0.6 * clamp01(best_identity));
    explanation.transduction_posterior = logistic(transduction_state_score);

    const double identity = clamp01(best_identity);
    const double coverage = explanation.te_core_coverage;
    const double margin_signal = clamp01(cross_family_margin * 4.0);
    const double residual_fraction = clamp01(annotation_residual_fraction);
    const double high_complexity_residual_fraction = clamp01(high_complexity_residual);
    const double masked_fraction = clamp01(annotation_masked_fraction);
    const double ambiguity_penalty =
        explanation.status == TeAnnotationStatus::kUnknownTe ? 0.45 :
        explanation.status == TeAnnotationStatus::kFamilyOnly ? 0.25 : 0.0;
    double model_support = 0.0;
    if (sequence_model_label == "TE_MODEL_IN_DISTRIBUTION") {
        model_support = 0.35 + (0.25 * clamp01(sequence_model_score));
    } else if (sequence_model_label == "TE_MODEL_EDGE") {
        model_support = -0.15;
    } else if (sequence_model_label == "TE_MODEL_OUTLIER") {
        model_support = -0.85;
    }

    explanation.te_structure_log_evidence =
        -1.10 +
        (2.00 * identity) +
        (2.45 * coverage) +
        (0.55 * margin_signal) +
        (0.55 * explanation.polyA_posterior) +
        (0.35 * explanation.transduction_posterior) +
        model_support -
        (1.25 * high_complexity_residual_fraction) -
        (0.55 * masked_fraction) -
        ambiguity_penalty;
    if (!has_te_path(explanation.status)) {
        explanation.te_structure_log_evidence -= 2.0;
    }
    explanation.nonte_structure_log_evidence =
        -0.25 +
        (1.65 * high_complexity_residual_fraction) +
        (0.75 * residual_fraction) -
        (1.25 * coverage) -
        (0.85 * identity) -
        (0.25 * explanation.polyA_posterior);
    explanation.artifact_structure_log_evidence =
        -0.10 +
        (1.10 * masked_fraction) +
        (sequence_model_label == "TE_MODEL_OUTLIER" ? 1.25 : 0.0) +
        (qc_reason == "TE_ALIGNMENT_LOW_IDENTITY" ? 0.65 : 0.0) -
        (0.95 * coverage) -
        (0.75 * identity);
    const double best_null = std::max(
        explanation.nonte_structure_log_evidence,
        explanation.artifact_structure_log_evidence);
    explanation.structure_path_confidence = logistic(
        explanation.te_structure_log_evidence - best_null);

    if (has_te_path(explanation.status) && insert_len > 0) {
        explanation.path.push_back({
            "TE_CORE",
            0,
            te_core_bases,
            te_path_label(explanation.status, explanation.family, explanation.subfamily),
        });
        int32_t cursor = te_core_bases;
        const int32_t poly_start = insert_len - poly_bases;
        // Viterbi include/exclude: open the optional state iff its decode margin
        // (emission log-LR + open log-odds) is positive.
        const bool add_transduction =
            residual_before_poly > 0 && transduction_state_score > 0.0;
        const bool add_poly = poly_bases > 0 && poly_state_score > 0.0;
        if (add_transduction && cursor < poly_start) {
            explanation.path.push_back({
                "TRANSDUCTION",
                cursor,
                poly_start,
                "HIGH_COMPLEXITY_TAIL",
            });
            cursor = poly_start;
        }
        if (add_poly && cursor < insert_len) {
            explanation.path.push_back({
                "POLYA",
                std::max(cursor, poly_start),
                insert_len,
                insert_seq.empty() ? "POLY_N" : std::string("POLY_") + insert_seq.back(),
            });
            cursor = insert_len;
        }
        if (masked_fraction > 0.20 && cursor < insert_len) {
            explanation.path.push_back({
                "LOW_COMPLEXITY",
                cursor,
                insert_len,
                "MASKED_RESIDUAL",
            });
            cursor = insert_len;
        }
        if (cursor < insert_len) {
            explanation.path.push_back({
                "UNEXPLAINED",
                cursor,
                insert_len,
                "HIGH_COMPLEXITY_RESIDUAL",
            });
        }
    }

    return explanation;
}

TeSequenceExplanation explain_te_alignment_shadow(
    int32_t insert_len,
    const std::string& qc_reason,
    const std::string& family,
    const std::string& subfamily,
    double best_identity,
    double effective_query_coverage,
    double annotation_residual_fraction,
    double annotation_masked_fraction,
    double cross_family_margin,
    double second_score) {
    return explain_te_sequence_structure(
        std::string(static_cast<size_t>(std::max(0, insert_len)), 'N'),
        qc_reason,
        family,
        subfamily,
        best_identity,
        effective_query_coverage,
        annotation_residual_fraction,
        annotation_masked_fraction,
        cross_family_margin,
        second_score,
        "TE_MODEL_UNAVAILABLE",
        0.0);
}

const char* te_annotation_status_name(TeAnnotationStatus status) {
    switch (status) {
        case TeAnnotationStatus::kUnavailable: return "UNAVAILABLE";
        case TeAnnotationStatus::kEmpty: return "EMPTY";
        case TeAnnotationStatus::kTooShort: return "TOO_SHORT";
        case TeAnnotationStatus::kNoCandidate: return "NO_CANDIDATE";
        case TeAnnotationStatus::kResolved: return "RESOLVED";
        case TeAnnotationStatus::kFamilyOnly: return "FAMILY_ONLY";
        case TeAnnotationStatus::kUnknownTe: return "UNKNOWN_TE";
        case TeAnnotationStatus::kAmbiguous: return "AMBIGUOUS";
        case TeAnnotationStatus::kNonTeLike: return "NON_TE_LIKE";
    }
    return "UNAVAILABLE";
}

std::string serialize_te_path_residual(const TePathResidual& residual) {
    std::ostringstream out;
    out << "unexplained=" << residual.unexplained_high_complexity_bases
        << ";edit_distance=" << residual.edit_distance
        << ";family_conflicts=" << residual.family_conflicts
        << ";subfamily_conflicts=" << residual.subfamily_conflicts
        << ";orientation_conflicts=" << residual.orientation_conflicts
        << ";segment_breaks=" << residual.segment_breaks
        << ";low_complexity_only=" << residual.low_complexity_only_bases;
    return out.str();
}

std::string serialize_te_sequence_path(const TeSequenceExplanation& explanation) {
    if (explanation.path.empty()) {
        return "NA";
    }
    std::ostringstream out;
    for (size_t i = 0; i < explanation.path.size(); ++i) {
        if (i > 0) {
            out << ",";
        }
        const TeSequenceSegment& segment = explanation.path[i];
        out << segment.query_start << "-" << segment.query_end
            << ":" << segment.kind
            << ":" << segment.label;
    }
    return out.str();
}

}  // namespace placer
