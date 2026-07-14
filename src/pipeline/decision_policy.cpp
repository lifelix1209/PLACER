#include "decision_policy.h"
#include "decision_thresholds.h"
#include "mechanistic_evidence.h"
#include "pipeline.h"
#include "te_alignment_qc.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <string>

namespace placer {

namespace {

using DT = DecisionThresholds;

double safe_log(double p) {
    return std::log(std::max(1e-12, p));
}

double phred_from_error_probability(double p_error) {
    const double bounded = std::clamp(p_error, 1e-12, 1.0);
    return -10.0 * std::log10(bounded);
}

double logsumexp3(double a, double b, double c) {
    const double m = std::max(a, std::max(b, c));
    return m + std::log(std::exp(a - m) + std::exp(b - m) + std::exp(c - m));
}

int32_t median_i32(std::vector<int32_t> values) {
    if (values.empty()) {
        return 0;
    }
    std::sort(values.begin(), values.end());
    return values[values.size() / 2];
}

int32_t infer_event_length(const EventGenotypeInput& input) {
    if (input.event_length > 0) {
        return input.event_length;
    }
    return median_i32(input.alt_observed_lengths);
}

double beta_binomial_log_pmf(int32_t alt, int32_t total, double alpha, double beta);

double alt_length_support_probability(
    int32_t observed_length,
    int32_t event_length,
    double error_rate) {
    if (event_length <= 0 || observed_length <= 0) {
        return std::clamp(1.0 - (error_rate * 2.0), 0.55, 0.95);
    }
    const double norm_factor = (error_rate <= 0.03) ? 20.0 : 30.0;
    const double sigma = std::max(1.0, static_cast<double>(event_length) / norm_factor);
    const double z =
        (static_cast<double>(observed_length) - static_cast<double>(event_length)) / sigma;
    return std::max(1e-6, std::exp(-0.5 * z * z));
}

// Mean length-support probability of the precise alt reads (1.0 when no lengths
// are observed). Low values mean the alt reads do not match the event length --
// i.e. they look like a different event or artifact rather than genuine support.
double length_concordance_factor(const EventGenotypeInput& input) {
    if (input.alt_observed_lengths.empty()) {
        return 1.0;
    }
    const int32_t event_length = infer_event_length(input);
    double acc = 0.0;
    for (int32_t observed_length : input.alt_observed_lengths) {
        acc += alt_length_support_probability(observed_length, event_length, input.error_rate);
    }
    return std::clamp(
        acc / static_cast<double>(input.alt_observed_lengths.size()), 1e-3, 1.0);
}

// Overdispersed genotype likelihood: alt | depth, GT ~ BetaBinomial(depth,
// alpha_GT, beta_GT), where the per-genotype expected alt-read fraction is
// mu_00 = error rate (alt only from error/artifact), mu_01 = alt_copy_fraction
// (~0.5), mu_11 = 1 - error rate. The overdispersion rho (intra-class
// correlation) sets the concentration kappa = (1-rho)/rho; rho -> 0 recovers the
// independent-read binomial limit. A length-concordance factor down-weights the
// non-reference genotypes when the precise alt reads do not match the event
// length. This replaces the previous independent-per-read Bernoulli product,
// which understated count variance in repeats.
double genotype_log_likelihood(
    const EventGenotypeInput& input,
    double alt_copy_fraction) {
    const int32_t alt = std::max(0, input.alt_struct_reads);
    const int32_t ref = std::max(0, input.ref_span_reads);
    const int32_t depth = alt + ref;
    if (depth <= 0) {
        return 0.0;
    }

    const double err = std::clamp(input.error_rate, 1e-4, 0.2);
    double mu;
    if (alt_copy_fraction <= 0.0) {
        mu = err;
    } else if (alt_copy_fraction >= 1.0) {
        mu = 1.0 - err;
    } else {
        mu = std::clamp(alt_copy_fraction, err, 1.0 - err);
    }

    const double rho = std::clamp(input.overdispersion, 0.0, 0.99);
    const double kappa = rho <= 1e-9 ? 1e9 : (1.0 - rho) / rho;
    const double alpha = std::max(1e-6, mu * kappa);
    const double beta = std::max(1e-6, (1.0 - mu) * kappa);

    double ll = beta_binomial_log_pmf(alt, depth, alpha, beta);

    // Non-reference genotypes claim the alt reads are genuine; discount that claim
    // by the per-read log length-concordance so length-discordant "alt" reads do
    // not inflate the het/hom likelihoods.
    if (alt_copy_fraction > 0.0 && alt > 0) {
        ll += static_cast<double>(alt) * std::log(length_concordance_factor(input));
    }
    return ll;
}

double clamp_score(double value, double lo = -3.0, double hi = 3.0) {
    return std::clamp(value, lo, hi);
}

double joint_event_existence_score(const EventExistenceEvidence& existence) {
    constexpr int32_t kUnopposedAltMinReads = 3;
    if (existence.ref_span_reads == 0 &&
        existence.alt_struct_reads >= kUnopposedAltMinReads) {
        return std::max(0.0, existence.score);
    }
    return existence.score;
}

bool is_one_sided_segmentation_pass(const EventSegmentationEvidence& segmentation) {
    return segmentation.pair_valid &&
           segmentation.has_insert_seq &&
           (segmentation.has_left_flank != segmentation.has_right_flank);
}

bool has_structural_breakdown(const EventExistenceEvidence& existence) {
    return existence.alt_split_reads >= 0 ||
           existence.alt_indel_reads >= 0 ||
           existence.alt_left_clip_reads >= 0 ||
           existence.alt_right_clip_reads >= 0;
}

int32_t precise_structural_reads(const EventExistenceEvidence& existence) {
    if (!has_structural_breakdown(existence)) {
        return existence.alt_struct_reads;
    }
    return std::max(0, existence.alt_split_reads) +
           std::max(0, existence.alt_indel_reads);
}

int32_t bilateral_clip_support_reads(const EventExistenceEvidence& existence) {
    if (!has_structural_breakdown(existence)) {
        return existence.alt_struct_reads;
    }
    return std::min(
        std::max(0, existence.alt_left_clip_reads),
        std::max(0, existence.alt_right_clip_reads));
}

double adjusted_segmentation_score(const EventSegmentationEvidence& segmentation);

bool annotation_quality_allows_te_decision(
    const TEAlignmentEvidence& te_alignment) {
    constexpr double kLowAnnotationMaxResidual = 0.50;
    constexpr double kLowAnnotationMaxMasked = 0.65;
    constexpr double kMediumAnnotationMaxResidual = 0.65;
    constexpr double kMediumAnnotationMaxMasked = 0.80;

    if (te_alignment.annotation_confidence == "LOW") {
        return te_alignment.annotation_residual_fraction <= kLowAnnotationMaxResidual &&
               te_alignment.annotation_masked_fraction <= kLowAnnotationMaxMasked;
    }
    if (te_alignment.annotation_confidence == "MEDIUM") {
        return te_alignment.annotation_residual_fraction <= kMediumAnnotationMaxResidual &&
               te_alignment.annotation_masked_fraction <= kMediumAnnotationMaxMasked;
    }
    return true;
}

double adjusted_nonref_existence_score(
    const EventExistenceEvidence& existence,
    const EventSegmentationEvidence& segmentation,
    double te_model_support) {
    double adjusted = joint_event_existence_score(existence);
    const bool has_strong_structural_te_support =
        segmentation.pair_valid &&
        segmentation.has_insert_seq &&
        te_model_support >= 0.75;
    const bool alt_not_weaker_than_ref =
        existence.alt_struct_reads > 0 &&
        existence.alt_struct_reads >= existence.ref_span_reads;
    if (has_strong_structural_te_support && alt_not_weaker_than_ref) {
        adjusted = std::max(0.0, adjusted);
    }
    return adjusted;
}

double adjusted_segmentation_score(const EventSegmentationEvidence& segmentation) {
    if (is_one_sided_segmentation_pass(segmentation)) {
        return std::max(0.25, segmentation.score);
    }
    return segmentation.score;
}

double adjusted_te_boundary_score(
    const EventSegmentationEvidence& segmentation,
    const BoundaryEvidence& boundary) {
    if (!boundary.geometry_defined && is_one_sided_segmentation_pass(segmentation)) {
        return 0.0;
    }
    return boundary.score;
}

double te_sequence_model_score(const TEAlignmentEvidence& te_alignment) {
    if (te_alignment.sequence_model_label == "TE_MODEL_UNAVAILABLE") {
        return 0.0;
    }
    if (te_alignment.sequence_model_label == "TE_MODEL_IN_DISTRIBUTION") {
        return 0.75 + clamp_score(te_alignment.sequence_model_score, 0.0, 0.50);
    }
    if (te_alignment.sequence_model_label == "TE_MODEL_EDGE") {
        return 0.0;
    }
    if (te_alignment.sequence_model_label == "TE_MODEL_OUTLIER") {
        return clamp_score(te_alignment.sequence_model_score, -1.0, -0.50);
    }
    return 0.0;
}

JointHypothesisScore make_hypothesis(FinalHypothesisKind kind) {
    JointHypothesisScore out;
    out.kind = kind;
    return out;
}

struct LocalHypothesisPosterior {
    double log_te = -1e9;
    double log_non_te = -1e9;
    double log_artifact = -1e9;
    double te = 0.0;
    double non_te = 0.0;
    double artifact = 0.0;
    double te_vs_artifact_log_odds = 0.0;
    double te_vs_non_te_log_odds = 0.0;
    std::string qc = "TE_POSTERIOR_LOW";
};

constexpr std::array<const char*, 6> kLatentMechanismStates = {
    "active_tprt_te",
    "cut_paste_te",
    "ltr_complex_te",
    "degraded_unknown_te",
    "non_te_insert",
    "artifact_reference"};

bool is_te_latent_state(const std::string& state) {
    return state == "active_tprt_te" ||
           state == "cut_paste_te" ||
           state == "ltr_complex_te" ||
           state == "degraded_unknown_te";
}

double count_signal(int32_t count, double scale) {
    if (count <= 0) {
        return 0.0;
    }
    return std::clamp(
        1.0 - std::exp(-static_cast<double>(count) / std::max(1e-6, scale)),
        0.0,
        1.0);
}

double beta_feature_loglike(double x, double mean, double strength) {
    constexpr double kEps = 1e-6;
    const double bounded_x = std::clamp(x, kEps, 1.0 - kEps);
    const double bounded_mean = std::clamp(mean, kEps, 1.0 - kEps);
    return strength *
           ((bounded_mean * std::log(bounded_x)) +
            ((1.0 - bounded_mean) * std::log(1.0 - bounded_x)));
}

double logsumexp_values(const std::array<double, 6>& values) {
    const double top = *std::max_element(values.begin(), values.end());
    double sum = 0.0;
    for (double value : values) {
        sum += std::exp(value - top);
    }
    return top + std::log(sum);
}

double logsumexp_pair(double lhs, double rhs) {
    const double top = std::max(lhs, rhs);
    return top + std::log(std::exp(lhs - top) + std::exp(rhs - top));
}

double log_choose_count(int32_t n, int32_t k) {
    if (k < 0 || n < 0 || k > n) {
        return -1e300;
    }
    return std::lgamma(static_cast<double>(n) + 1.0) -
           std::lgamma(static_cast<double>(k) + 1.0) -
           std::lgamma(static_cast<double>(n - k) + 1.0);
}

double beta_binomial_log_pmf(int32_t alt, int32_t total, double alpha, double beta) {
    if (alt < 0 || total < 0 || alt > total || alpha <= 0.0 || beta <= 0.0) {
        return -1e300;
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
        return -1e300;
    }
    const double clamped_p = std::clamp(p, 1e-6, 1.0 - 1e-6);
    return log_choose_count(total, alt) +
           (static_cast<double>(alt) * std::log(clamped_p)) +
           (static_cast<double>(total - alt) * std::log1p(-clamped_p));
}

double low_allele_fraction_insertion_log_bf(int32_t alt, int32_t ref) {
    const int32_t total = std::max(0, alt) + std::max(0, ref);
    if (alt <= 0 || total <= 0) {
        return -1e300;
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

double logistic_logit(double x) {
    constexpr double kEps = 1e-6;
    const double bounded = std::clamp(x, kEps, 1.0 - kEps);
    return std::log(bounded / (1.0 - bounded));
}

std::string lower_ascii(std::string value) {
    for (char& ch : value) {
        ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    }
    return value;
}

bool contains_token(const std::string& text, const char* token) {
    return text.find(token) != std::string::npos;
}

std::string family_kind(const TEAlignmentEvidence& te_alignment) {
    const std::string family = lower_ascii(te_alignment.best_family);
    const std::string subfamily = lower_ascii(te_alignment.best_subfamily);
    const std::string name = family + ":" + subfamily;
    const bool family_unknown =
        family.empty() ||
        family == "unknown" ||
        family == "na" ||
        family == "none";
    const bool subfamily_unknown =
        subfamily.empty() ||
        subfamily == "unknown" ||
        subfamily == "na" ||
        subfamily == "none";
    if (family_unknown && subfamily_unknown) {
        return "unknown";
    }
    if (contains_token(name, "ltr") ||
        contains_token(name, "erv") ||
        contains_token(name, "gypsy") ||
        contains_token(name, "copia") ||
        contains_token(name, "bel-pao") ||
        contains_token(name, "dirs")) {
        return "ltr";
    }
    if (contains_token(name, "dna") ||
        contains_token(name, "hat") ||
        contains_token(name, "pif") ||
        contains_token(name, "harbinger") ||
        contains_token(name, "tc1") ||
        contains_token(name, "mariner") ||
        contains_token(name, "piggybac") ||
        contains_token(name, "mutator") ||
        contains_token(name, "merlin") ||
        contains_token(name, "cmc") ||
        contains_token(name, "dong")) {
        return "dna";
    }
    if (contains_token(name, "l1") ||
        contains_token(name, "l2") ||
        contains_token(name, "line") ||
        contains_token(name, "sine") ||
        contains_token(name, "alu") ||
        contains_token(name, "sva") ||
        contains_token(name, "rte") ||
        contains_token(name, "r2") ||
        contains_token(name, "cr1") ||
        contains_token(name, "penelope") ||
        contains_token(name, "5s-deu-l2")) {
        return "retro";
    }
    return "other";
}

struct LatentFeatureVector {
    double te_identity = 0.0;
    double event_signal = 0.0;
    double independent_signal = 0.0;
    double quality_signal = 0.0;
    double ref_conflict_signal = 0.0;
    double artifact_context_signal = 0.0;
    double family_activity_prior = 0.0;
    std::string family_kind = "unknown";
};

LatentFeatureVector build_latent_feature_vector(
    const EventExistenceEvidence& existence,
    const EventSegmentationEvidence& segmentation,
    const TEAlignmentEvidence& te_alignment,
    const BoundaryEvidence& boundary,
    const ClipInsertConcordanceEvidence* clip_insert_concordance) {
    LatentFeatureVector features;
    const int32_t alt = std::max(0, existence.alt_struct_reads);
    const int32_t ref = std::max(0, existence.ref_span_reads);
    const int32_t precise = precise_structural_reads(existence);
    const int32_t bilateral_clip = bilateral_clip_support_reads(existence);
    int32_t mechanistic_reads = precise + bilateral_clip;
    if (clip_insert_concordance != nullptr && clip_insert_concordance->pass) {
        mechanistic_reads += std::max(0, clip_insert_concordance->full_insert_reads);
        mechanistic_reads += std::min(
            std::max(0, clip_insert_concordance->left_clip_reads),
            std::max(0, clip_insert_concordance->right_clip_reads));
    }

    features.te_identity = std::clamp(te_alignment.best_identity, 0.0, 1.0);
    const double independent_count_signal = count_signal(mechanistic_reads, 4.0);
    const double independent_fraction_signal = alt > 0
        ? std::clamp(static_cast<double>(mechanistic_reads) / static_cast<double>(alt), 0.0, 1.0)
        : 0.0;
    features.independent_signal = std::sqrt(
        independent_count_signal * independent_fraction_signal);
    features.quality_signal = std::clamp(static_cast<double>(existence.gq) / 60.0, 0.0, 1.0);
    const double support_signal = count_signal(alt, 8.0);
    const double segmentation_signal = std::clamp(
        (adjusted_segmentation_score(segmentation) + 1.0) / 3.0,
        0.0,
        1.0);
    features.event_signal = std::clamp(
        (0.35 * support_signal) +
        (0.30 * features.independent_signal) +
        (0.20 * segmentation_signal) +
        (0.15 * std::clamp(existence.af, 0.0, 1.0)),
        0.0,
        1.0);

    const double ref_fraction = (alt + ref) > 0
        ? static_cast<double>(ref) / static_cast<double>(alt + ref)
        : 0.0;
    features.ref_conflict_signal = std::max(
        ref_fraction,
        0.60 * count_signal(ref, 8.0));
    if (is_one_sided_segmentation_pass(segmentation) && ref > 0) {
        features.ref_conflict_signal = std::max(features.ref_conflict_signal, 0.35);
    }
    features.ref_conflict_signal = std::clamp(features.ref_conflict_signal, 0.0, 1.0);

    double artifact_context = 0.0;
    if (!segmentation.has_insert_seq) {
        artifact_context = std::max(artifact_context, 1.0);
    }
    if (te_alignment.sequence_model_label == "TE_MODEL_OUTLIER") {
        artifact_context = std::max(artifact_context, 1.0);
    }
    if (te_alignment.sequence_model_label == "TE_MODEL_EDGE") {
        artifact_context = std::max(artifact_context, 0.45);
    }
    if (!boundary.geometry_defined && !is_one_sided_segmentation_pass(segmentation)) {
        artifact_context = std::max(artifact_context, 0.35);
    }
    if (boundary.geometry_defined && !boundary.canonical_pass && !boundary.evidence_consistent) {
        artifact_context = std::max(artifact_context, 0.50);
    }
    if (te_alignment.annotation_confidence == "LOW") {
        artifact_context = std::max(artifact_context, 0.35);
    }
    artifact_context = std::max(artifact_context, 0.80 * features.ref_conflict_signal);
    features.artifact_context_signal = std::clamp(artifact_context, 0.0, 1.0);

    features.family_kind = family_kind(te_alignment);
    const double family_known = features.family_kind == "unknown" ? 0.0 : 1.0;
    features.family_activity_prior = std::clamp(
        (0.30 * features.te_identity) +
        (0.25 * features.event_signal) +
        (0.20 * features.independent_signal) +
        (0.15 * features.quality_signal) +
        (0.10 * family_known) -
        (0.30 * features.artifact_context_signal) -
        (0.10 * features.ref_conflict_signal),
        0.02,
        0.98);
    return features;
}

double family_state_compatibility(
    const LatentFeatureVector& features,
    const std::string& state) {
    const std::string& kind = features.family_kind;
    if (state == "active_tprt_te") {
        if (kind == "retro") {
            return 0.75;
        }
        if (kind == "ltr") {
            return 0.20;
        }
        if (kind == "dna") {
            return -0.55;
        }
        if (kind == "unknown") {
            return -0.35;
        }
    }
    if (state == "cut_paste_te") {
        if (kind == "dna") {
            return 0.75;
        }
        if (kind == "retro" || kind == "ltr") {
            return -0.40;
        }
        if (kind == "unknown") {
            return -0.20;
        }
    }
    if (state == "ltr_complex_te") {
        if (kind == "ltr") {
            return 0.75;
        }
        if (kind == "retro") {
            return -0.10;
        }
        if (kind == "dna") {
            return -0.50;
        }
        if (kind == "unknown") {
            return -0.30;
        }
    }
    if (state == "degraded_unknown_te") {
        if (kind == "unknown") {
            return 0.40;
        }
        if (features.te_identity < 0.75) {
            return 0.20;
        }
    }
    return 0.0;
}

double latent_state_log_prior(
    const LatentFeatureVector& features,
    const std::string& state) {
    double base = -0.75;
    if (state == "active_tprt_te" || state == "cut_paste_te") {
        base = -1.35;
    } else if (state == "ltr_complex_te") {
        base = -1.50;
    } else if (state == "degraded_unknown_te") {
        base = -1.80;
    } else if (state == "non_te_insert") {
        base = -0.85;
    }
    const double activity_logit = logistic_logit(features.family_activity_prior);
    if (is_te_latent_state(state)) {
        base += 0.55 * activity_logit;
    } else {
        base -= 0.20 * activity_logit;
    }
    return base + family_state_compatibility(features, state);
}

double latent_state_feature_log_likelihood(
    const LatentFeatureVector& features,
    const std::string& state) {
    struct Params {
        double te_mean;
        double te_strength;
        double event_mean;
        double event_strength;
        double independent_mean;
        double independent_strength;
        double quality_mean;
        double quality_strength;
        double artifact_mean;
        double artifact_strength;
        double ref_mean;
        double ref_strength;
    };
    Params params{0.38, 3.5, 0.22, 4.0, 0.12, 4.0, 0.35, 3.0, 0.72, 5.0, 0.55, 3.0};
    if (state == "active_tprt_te") {
        params = {0.94, 12.0, 0.86, 5.5, 0.80, 6.0, 0.88, 3.0, 0.06, 4.5, 0.08, 2.5};
    } else if (state == "cut_paste_te") {
        params = {0.92, 11.0, 0.84, 5.5, 0.76, 6.0, 0.86, 3.0, 0.07, 4.5, 0.08, 2.5};
    } else if (state == "ltr_complex_te") {
        params = {0.90, 10.0, 0.80, 5.0, 0.70, 5.5, 0.84, 3.0, 0.10, 4.0, 0.12, 2.5};
    } else if (state == "degraded_unknown_te") {
        params = {0.72, 3.5, 0.62, 4.0, 0.50, 5.0, 0.72, 2.5, 0.22, 3.0, 0.20, 2.0};
    } else if (state == "non_te_insert") {
        params = {0.30, 5.5, 0.82, 5.0, 0.68, 4.0, 0.82, 3.0, 0.12, 3.5, 0.12, 2.0};
    }
    double log_like = 0.0;
    log_like += beta_feature_loglike(features.te_identity, params.te_mean, params.te_strength);
    log_like += beta_feature_loglike(features.event_signal, params.event_mean, params.event_strength);
    log_like += beta_feature_loglike(
        features.independent_signal,
        params.independent_mean,
        params.independent_strength);
    log_like += beta_feature_loglike(features.quality_signal, params.quality_mean, params.quality_strength);
    log_like += beta_feature_loglike(
        features.artifact_context_signal,
        params.artifact_mean,
        params.artifact_strength);
    log_like += beta_feature_loglike(
        features.ref_conflict_signal,
        params.ref_mean,
        params.ref_strength);
    return log_like;
}

bool is_te_interpretable_without_threshold_gate(
    const TEAlignmentEvidence& te_alignment) {
    if (!te_alignment.pass) {
        return false;
    }
    if (te_alignment.sequence_model_label == "TE_MODEL_OUTLIER") {
        return false;
    }
    return is_pass_insert_te_alignment(
        classify_te_alignment_qc(te_alignment.qc_reason));
}

bool has_closed_te_breakpoints(
    const EventSegmentationEvidence& segmentation,
    const BoundaryEvidence& boundary) {
    return segmentation.pair_valid &&
           segmentation.has_left_flank &&
           segmentation.has_right_flank &&
           boundary.geometry_defined &&
           (boundary.canonical_pass || boundary.evidence_consistent);
}

double te_alignment_log_support(const TEAlignmentEvidence& te_alignment) {
    const double identity = std::clamp(te_alignment.best_identity, 0.0, 1.0);
    const double coverage = std::clamp(te_alignment.best_query_coverage, 0.0, 1.0);
    const double margin = std::clamp(te_alignment.cross_family_margin, 0.0, 1.0);

    double score = -1.5;
    switch (classify_te_alignment_qc(te_alignment.qc_reason)) {
        case TeAlignmentQc::kPassInsert:
            score = 1.2 +
                    (3.0 * (identity - 0.78)) +
                    (1.4 * (coverage - 0.72)) +
                    std::min(1.0, 2.0 * margin);
            break;
        case TeAlignmentQc::kPassInsertFamilyOnly:
            score = 0.65 +
                    (2.2 * (identity - 0.68)) +
                    (1.0 * (coverage - 0.55)) +
                    std::min(0.7, 1.5 * margin);
            break;
        case TeAlignmentQc::kPassInsertUnknown:
            score = -0.45 +
                    (1.9 * (identity - 0.55)) +
                    (0.9 * (coverage - 0.70)) +
                    std::min(0.45, 1.0 * margin);
            break;
        case TeAlignmentQc::kLowIdentity:
            score = -0.85 +
                    (1.4 * (identity - 0.50)) +
                    (0.6 * (coverage - 0.40));
            break;
        case TeAlignmentQc::kNoTeAlignment:
        case TeAlignmentQc::kOther:
            if (!te_alignment.pass) {
                score = -2.0;
            }
            break;
    }

    if (te_alignment.sequence_model_label == "TE_MODEL_IN_DISTRIBUTION") {
        score += 0.55 + clamp_score(te_alignment.sequence_model_score, 0.0, 0.50);
    } else if (te_alignment.sequence_model_label == "TE_MODEL_EDGE") {
        score -= 0.25;
    } else if (te_alignment.sequence_model_label == "TE_MODEL_OUTLIER") {
        score -= 1.0;
    }

    if (te_alignment.annotation_confidence == "LOW") {
        score -= 0.65;
    } else if (te_alignment.annotation_confidence == "MEDIUM") {
        score -= 0.20;
    } else if (te_alignment.annotation_confidence == "HIGH") {
        score += 0.25;
    }
    score -= std::clamp(te_alignment.annotation_residual_fraction, 0.0, 1.0) * 0.8;
    score -= std::clamp(te_alignment.annotation_masked_fraction, 0.0, 1.0) * 0.3;
    return clamp_score(score, -3.0, 3.0);
}

LocalHypothesisPosterior evaluate_local_hypothesis_posterior(
    const EventExistenceEvidence& existence,
    const EventSegmentationEvidence& segmentation,
    const TEAlignmentEvidence& te_alignment,
    const BoundaryEvidence& boundary) {
    LocalHypothesisPosterior posterior;
    const int32_t alt = std::max(0, existence.alt_struct_reads);
    const int32_t ref = std::max(0, existence.ref_span_reads);
    const int32_t precise = precise_structural_reads(existence);
    const int32_t bilateral_clip = bilateral_clip_support_reads(existence);
    const bool one_sided = is_one_sided_segmentation_pass(segmentation);
    const bool closed = has_closed_te_breakpoints(segmentation, boundary);
    const bool unknown_te =
        classify_te_alignment_qc(te_alignment.qc_reason) ==
            TeAlignmentQc::kPassInsertUnknown ||
        te_alignment.best_family == "UNKNOWN" ||
        te_alignment.best_subfamily == "UNKNOWN";

    const double count_balance =
        clamp_score(std::log(static_cast<double>(alt + 1)) -
                    std::log(static_cast<double>(ref + 1)), -2.0, 2.0);
    const double precise_balance =
        clamp_score(std::log(static_cast<double>(precise + bilateral_clip + 1)), 0.0, 2.0);
    const double genotype_support =
        clamp_score((static_cast<double>(existence.gq) - 20.0) / 40.0, 0.0, 1.5);
    const double structural_support =
        clamp_score(0.65 * count_balance + 0.45 * precise_balance + genotype_support, -2.0, 3.0);
    const double segmentation_support = clamp_score(adjusted_segmentation_score(segmentation), -2.0, 2.0);
    const double sequence_support = te_alignment_log_support(te_alignment);
    const double te_identity = std::clamp(te_alignment.best_identity, 0.0, 1.0);
    const double te_coverage = std::clamp(te_alignment.best_query_coverage, 0.0, 1.0);
    const double te_margin = std::clamp(te_alignment.cross_family_margin, 0.0, 1.0);

    double boundary_support = -1.0;
    if (closed) {
        boundary_support = 1.0;
    } else if (one_sided && segmentation.has_insert_seq) {
        boundary_support = -0.35;
    } else if (boundary.geometry_defined && boundary.evidence_consistent) {
        boundary_support = 0.25;
    }

    double artifact_context = 0.0;
    if (!segmentation.has_insert_seq) {
        artifact_context += 1.6;
    }
    if (one_sided) {
        artifact_context += 1.1;
    }
    if (!boundary.geometry_defined || (!boundary.canonical_pass && !boundary.evidence_consistent)) {
        artifact_context += 0.9;
    }
    if (unknown_te) {
        const double identity_deficit =
            clamp_score((0.62 - te_identity) / 0.16, 0.0, 1.0);
        const double coverage_deficit =
            clamp_score((0.90 - te_coverage) / 0.40, 0.0, 1.0);
        const double margin_deficit =
            clamp_score((0.05 - te_margin) / 0.05, 0.0, 1.0);
        if (te_alignment.annotation_confidence == "LOW") {
            artifact_context += 0.50;
        } else if (te_alignment.annotation_confidence == "MEDIUM") {
            artifact_context += 0.20;
        }
        artifact_context += (0.45 * identity_deficit) +
                            (0.20 * coverage_deficit) +
                            (0.15 * margin_deficit);
    }
    if (precise == 0 && bilateral_clip == 0) {
        artifact_context += 0.5;
    }
    artifact_context += std::min(1.0, static_cast<double>(ref) * 0.08);

    posterior.log_te =
        -0.25 +
        structural_support +
        (0.75 * segmentation_support) +
        sequence_support +
        boundary_support -
        (0.35 * artifact_context);
    posterior.log_non_te =
        0.10 +
        (0.85 * structural_support) +
        (0.60 * segmentation_support) -
        (0.75 * sequence_support) +
        (segmentation.has_insert_seq ? 0.25 : -0.75);
    posterior.log_artifact =
        -0.10 +
        artifact_context -
        (0.25 * structural_support) -
        (0.15 * segmentation_support) -
        (0.35 * std::max(sequence_support, 0.0));

    const double denom = logsumexp3(
        posterior.log_te,
        posterior.log_non_te,
        posterior.log_artifact);
    posterior.te = std::exp(posterior.log_te - denom);
    posterior.non_te = std::exp(posterior.log_non_te - denom);
    posterior.artifact = std::exp(posterior.log_artifact - denom);
    posterior.te_vs_artifact_log_odds = posterior.log_te - posterior.log_artifact;
    posterior.te_vs_non_te_log_odds = posterior.log_te - posterior.log_non_te;

    const double min_te_posterior = closed ? 0.65 : 0.72;
    const double min_artifact_odds = closed ? 0.75 : 1.25;
    const double min_non_te_odds = unknown_te ? 0.55 : 0.25;
    if (posterior.te >= min_te_posterior &&
        posterior.te_vs_artifact_log_odds >= min_artifact_odds &&
        posterior.te_vs_non_te_log_odds >= min_non_te_odds) {
        posterior.qc = "PASS_TE_POSTERIOR";
    }
    return posterior;
}

// Emit a non-TE structural insertion (h1) call when the event is a closed
// structural event (both flanks + a usable boundary -- definitional, not a tuned
// ladder) and the assembled event-existence log-evidence clears a Bayes-factor
// cut-off while the artifact posterior stays low. This is a single log-evidence
// gate: no read-count / insert-length / segmentation-score ladders and no
// certificate fallback.
bool should_emit_structural_event_call(
    const EventExistenceEvidence& existence,
    const EventSegmentationEvidence& segmentation,
    const BoundaryEvidence& boundary,
    const LocalHypothesisPosterior& posterior,
    const MechanisticEvidenceCertificate& mechanistic_certificate,
    bool nonte_explanation) {
    const int32_t alt = std::max(0, existence.alt_struct_reads);
    const int32_t ref = std::max(0, existence.ref_span_reads);
    const int32_t precise = precise_structural_reads(existence);
    const int32_t bilateral_clip = bilateral_clip_support_reads(existence);
    const bool closed_structural_event =
        segmentation.pair_valid &&
        segmentation.has_insert_seq &&
        segmentation.has_left_flank &&
        segmentation.has_right_flank &&
        boundary.geometry_defined &&
        (boundary.canonical_pass || boundary.evidence_consistent);
    if (!closed_structural_event) {
        return false;
    }

    // Real-insertion-vs-artifact log odds: a structural insertion exists when the
    // event is a real insertion (TE *or* non-TE) rather than a reference/mapping
    // artifact. Using te+non_te (not non_te alone) is what lets a TE-like event
    // that abstained on the TE lFDR gate still be reported as a real insertion.
    const double real_insertion_posterior =
        std::clamp(posterior.te + posterior.non_te, 1e-6, 1.0);
    const double structural_log_odds =
        std::log(real_insertion_posterior) -
        std::log(std::max(std::clamp(posterior.artifact, 0.0, 1.0), 1e-6));
    // Count evidence for a real (non-artifact) insertion, valid across both
    // regimes: a dominant-allele log-odds for high-AF events, and a
    // beta-binomial insertion-vs-error Bayes factor for low-AF het/mosaic events
    // where alt < ref. Taking the max asks whether the counts are consistent
    // with a real allele under *either* model -- no allele-fraction ladder.
    const double dominant_allele_log_odds =
        std::log(static_cast<double>(alt + 1)) -
        std::log(static_cast<double>(ref + 1));
    const double low_af_log_bf = low_allele_fraction_insertion_log_bf(alt, ref);
    // Low-AF (het/mosaic) regime: the ALT allele is the minority (alt < ref) yet
    // the overdispersed insertion-vs-error Bayes factor still favours a real
    // allele. Otherwise the event is in the dominant-allele regime.
    const bool low_af_regime =
        alt < ref && std::isfinite(low_af_log_bf) &&
        low_af_log_bf > dominant_allele_log_odds;
    const double count_evidence =
        low_af_regime ? low_af_log_bf : dominant_allele_log_odds;
    // Low-AF het/mosaic insertions legitimately carry many reference-spanning
    // reads, so they are held to the looser artifact-posterior cap; dominant
    // -allele events use the strict cap.
    const double artifact_cap = low_af_regime
        ? DT::kArtifactPosteriorCap
        : DT::kArtifactPosteriorStrictCap;
    const double independent_support =
        std::log(static_cast<double>(precise + bilateral_clip + 1));
    const double boundary_support = boundary.canonical_pass ? 0.70 : 0.35;
    const double genotype_support =
        std::clamp((static_cast<double>(existence.gq) - 20.0) / 40.0, 0.0, 1.0);
    const double ref_conflict_penalty =
        1.20 * std::clamp(mechanistic_certificate.ref_conflict_signal, 0.0, 1.0);
    const double event_existence_log_evidence =
        structural_log_odds +
        (0.65 * count_evidence) +
        (0.35 * independent_support) +
        boundary_support +
        (0.40 * genotype_support) -
        ref_conflict_penalty;

    // The structural-insertion fallback catches real insertions the confident-TE
    // path misses: genuine non-TE insertions, and low-AF (het/mosaic) insertions
    // whose TE identity is uncertain only because ALT coverage is low. A high-AF
    // TE-like event that abstained on the TE lFDR gate is *not* downgraded to a
    // structural call -- it stays TE-evidence. Hence the regime guard.
    if (!nonte_explanation && !low_af_regime) {
        return false;
    }
    // Necessary conditions: decisive event-existence log-evidence, a low artifact
    // posterior, at least one precise/clip structural read, and -- crucially --
    // the event must genotype as a real variant (positive genotype support,
    // i.e. GQ beyond the Phred-20 confidence pivot). The genotype requirement
    // stops boundary/model priors alone from carrying a trivially-supported event
    // over the log-evidence bar.
    return event_existence_log_evidence > DT::kStructuralLogEvidence &&
           posterior.artifact <= artifact_cap &&
           (precise + bilateral_clip) > 0 &&
           genotype_support > 0.0;
}

TeSequenceExplanation te_structure_explanation_for_decision(
    const EventSegmentationEvidence& segmentation,
    const TEAlignmentEvidence& te_alignment) {
    if (te_alignment.te_sequence_explanation.status !=
        TeAnnotationStatus::kUnavailable) {
        return te_alignment.te_sequence_explanation;
    }
    return explain_te_alignment_shadow(
        segmentation.insert_len,
        te_alignment.qc_reason,
        te_alignment.best_family,
        te_alignment.best_subfamily,
        te_alignment.best_identity,
        te_alignment.best_query_coverage,
        te_alignment.annotation_residual_fraction,
        te_alignment.annotation_masked_fraction,
        te_alignment.cross_family_margin,
        te_alignment.second_score);
}

bool has_te_sequence_evidence(const TEAlignmentEvidence& te_alignment) {
    const TeAlignmentQc qc = classify_te_alignment_qc(te_alignment.qc_reason);
    return te_alignment.pass ||
           is_pass_insert_te_alignment(qc) ||
           qc == TeAlignmentQc::kLowIdentity ||
           te_alignment.best_identity > 0.0 ||
           te_alignment.best_query_coverage > 0.0;
}

int32_t positive_or_zero(int32_t value) {
    return std::max(0, value);
}

EventExplanation make_reference_explanation(
    const EventExistenceEvidence& existence,
    const EventSegmentationEvidence& segmentation) {
    EventExplanation out;
    out.kind = ExplanationKind::kReference;
    out.status = "REFERENCE";
    out.residual.structural_conflicts = segmentation.has_insert_seq ? 1 : 0;
    out.residual.unexplained_high_complexity_bases =
        segmentation.has_insert_seq ? std::max(1, segmentation.insert_len) : 0;
    out.residual.read_assignment_conflicts =
        positive_or_zero(existence.alt_struct_reads);
    out.residual.path_complexity = 0;
    return out;
}

EventExplanation make_non_te_explanation(
    const EventExistenceEvidence& existence,
    const EventSegmentationEvidence& segmentation,
    const TEAlignmentEvidence& te_alignment) {
    EventExplanation out;
    out.kind = ExplanationKind::kInsertionNonTe;
    out.status = "NON_TE_INSERTION";
    out.residual.missing_required_components = segmentation.has_insert_seq ? 0 : 1;
    out.residual.unexplained_high_complexity_bases =
        segmentation.has_insert_seq
            ? std::max(0, segmentation.insert_len / (te_alignment.pass ? 2 : 4))
            : 0;
    out.residual.read_assignment_conflicts =
        positive_or_zero(existence.ref_span_reads);
    out.residual.label_ambiguity = te_alignment.pass ? 1 : 0;
    out.residual.path_complexity = segmentation.has_insert_seq ? 1 : 0;
    return out;
}

EventExplanation make_te_explanation(
    const EventExistenceEvidence& existence,
    const EventSegmentationEvidence& segmentation,
    const TEAlignmentEvidence& te_alignment,
    const BoundaryEvidence& boundary) {
    EventExplanation out;
    out.kind = ExplanationKind::kTe;
    out.family = te_alignment.best_family.empty() ? "UNKNOWN" : te_alignment.best_family;
    out.subfamily = te_alignment.best_subfamily.empty() ? "NA" : te_alignment.best_subfamily;
    out.status = out.family == "UNKNOWN" ? "TE_UNKNOWN" : "TE_RESOLVED";
    out.residual.missing_required_components =
        (segmentation.has_insert_seq ? 0 : 1) +
        (is_te_interpretable_without_threshold_gate(te_alignment) ? 0 : 1);
    out.residual.structural_conflicts =
        segmentation.has_insert_seq ? 0 : 1;
    out.residual.unexplained_high_complexity_bases =
        is_te_interpretable_without_threshold_gate(te_alignment) ? 0 : std::max(1, segmentation.insert_len);
    out.residual.breakpoint_disagreement_bp =
        has_closed_te_breakpoints(segmentation, boundary) ? 0 : 25;
    out.residual.read_assignment_conflicts =
        positive_or_zero(existence.ref_span_reads);
    out.residual.artifact_evidence =
        te_alignment.sequence_model_label == "TE_MODEL_EDGE" ? 1 : 0;
    out.residual.label_ambiguity =
        out.family == "UNKNOWN" ? 1 : 0;
    out.residual.path_complexity = out.family == "UNKNOWN" ? 2 : 1;
    return out;
}

EventExplanation make_artifact_explanation(
    const EventExistenceEvidence& existence,
    const EventSegmentationEvidence& segmentation,
    const TEAlignmentEvidence& te_alignment,
    const BoundaryEvidence& boundary) {
    EventExplanation out;
    out.kind = ExplanationKind::kArtifact;
    out.status = "ARTIFACT";
    out.residual.unexplained_high_complexity_bases =
        segmentation.has_insert_seq ? std::max(0, segmentation.insert_len / 3) : 0;
    out.residual.read_assignment_conflicts =
        positive_or_zero(existence.alt_split_reads) +
        positive_or_zero(existence.alt_indel_reads);
    out.residual.reference_counterevidence =
        positive_or_zero(existence.ref_span_reads);
    out.residual.artifact_evidence =
        (te_alignment.sequence_model_label == "TE_MODEL_OUTLIER" ? 0 : 2) +
        (boundary.geometry_defined ? 1 : 0);
    out.residual.path_complexity = 1;
    return out;
}

}  // namespace

EventGenotypeDecision genotype_event_from_alt_vs_ref(
    const EventGenotypeInput& input) {
    EventGenotypeDecision decision;

    const int32_t alt = std::max(0, input.alt_struct_reads);
    const int32_t ref = std::max(0, input.ref_span_reads);
    const int32_t depth = alt + ref;
    decision.depth = depth;
    if (depth <= 0) {
        return decision;
    }

    decision.allele_fraction = std::clamp(
        static_cast<double>(alt) / static_cast<double>(depth),
        0.0,
        1.0);

    const double ll_00 = genotype_log_likelihood(input, 0.0) + safe_log(0.25);
    const double ll_01 = genotype_log_likelihood(input, 0.5) + safe_log(0.50);
    const double ll_11 = genotype_log_likelihood(input, 1.0) + safe_log(0.25);

    const bool het_is_best_nonref = ll_01 >= ll_11;
    const double best_nonref_ll = het_is_best_nonref ? ll_01 : ll_11;
    decision.best_gt = het_is_best_nonref ? "0/1" : "1/1";

    if (best_nonref_ll <= ll_00) {
        decision.best_gt = "0/0";
        decision.gq = 0;
        decision.pass = false;
        return decision;
    }

    const double delta_ll = best_nonref_ll - ll_00;
    decision.best_nonref_minus_ref_ll = delta_ll;
    const double best_total_ll = std::max(ll_00, std::max(ll_01, ll_11));
    const double posterior_best = std::exp(best_total_ll - logsumexp3(ll_00, ll_01, ll_11));
    const double gq = phred_from_error_probability(1.0 - posterior_best);
    decision.gq = std::max(0, std::min(99, static_cast<int32_t>(std::lround(gq))));
    decision.pass = decision.gq >= std::max(0, input.min_gq);
    return decision;
}

double estimate_alt_depth_overdispersion(
    const std::vector<AltDepthObservation>& observations,
    double fallback) {
    // Keep het-like informative sites (both alleles present at usable depth):
    // these carry the allele-balance overdispersion the genotyper cares about.
    std::vector<AltDepthObservation> obs;
    obs.reserve(observations.size());
    for (const AltDepthObservation& o : observations) {
        if (o.depth < 4 || o.alt < 0 || o.alt > o.depth) {
            continue;
        }
        const double p = static_cast<double>(o.alt) / static_cast<double>(o.depth);
        if (p > 0.15 && p < 0.85) {
            obs.push_back(o);
        }
    }
    const size_t m = obs.size();
    if (m < 20) {
        return fallback;  // too little information to estimate rho
    }

    double sum_x = 0.0;
    double sum_n = 0.0;
    double sum_n2 = 0.0;
    for (const AltDepthObservation& o : obs) {
        sum_x += o.alt;
        sum_n += o.depth;
        sum_n2 += static_cast<double>(o.depth) * static_cast<double>(o.depth);
    }
    const double p_hat = sum_x / sum_n;
    if (p_hat <= 1e-6 || p_hat >= 1.0 - 1e-6) {
        return fallback;
    }

    // One-way random-effects (Fleiss / ANOVA) estimator of the intra-class
    // correlation for proportions: rho = (MSB - MSW) / (MSB + (n0 - 1) MSW).
    double s_between = 0.0;
    double sum_within = 0.0;
    for (const AltDepthObservation& o : obs) {
        const double p = static_cast<double>(o.alt) / static_cast<double>(o.depth);
        s_between += static_cast<double>(o.depth) * (p - p_hat) * (p - p_hat);
        sum_within +=
            static_cast<double>(o.alt) * static_cast<double>(o.depth - o.alt) /
            static_cast<double>(o.depth);
    }
    const double N = sum_n;
    const double n0 = (N - (sum_n2 / N)) / (static_cast<double>(m) - 1.0);
    const double msb = s_between / (static_cast<double>(m) - 1.0);
    const double msw = sum_within / (N - static_cast<double>(m));
    const double denom = msb + ((n0 - 1.0) * msw);
    if (denom <= 1e-9) {
        return fallback;
    }
    return std::clamp((msb - msw) / denom, 0.0, 0.5);
}

FinalBoundaryDecision check_boundary_consistency(
    const FinalBoundaryInput& input) {
    FinalBoundaryDecision decision;

    if (input.left_ref_start < 0 ||
        input.left_ref_end < 0 ||
        input.right_ref_start < 0 ||
        input.right_ref_end < 0) {
        decision.qc = "REJECT_BOUNDARY_MISSING_REF_SEGMENTS";
        return decision;
    }
    if (input.left_ref_start >= input.left_ref_end ||
        input.right_ref_start >= input.right_ref_end) {
        decision.qc = "REJECT_BOUNDARY_INVALID_REF_SEGMENTS";
        return decision;
    }

    const int32_t tsd_min_len = std::max(1, input.tsd_min_len);
    const int32_t tsd_max_len = std::max(tsd_min_len, input.tsd_max_len);
    const int32_t link_delta = input.right_ref_start - input.left_ref_end;

    if (link_delta < 0) {
        const int32_t overlap = -link_delta;
        if (overlap < tsd_min_len || overlap > tsd_max_len) {
            decision.qc = "REJECT_BOUNDARY_TSD_RANGE";
            return decision;
        }
        decision.pass = true;
        decision.boundary_type = "TSD";
        decision.boundary_len = overlap;
        decision.qc = "PASS_BOUNDARY_TSD";
        return decision;
    }

    if (link_delta == 0) {
        decision.pass = true;
        decision.boundary_type = "BLUNT";
        decision.boundary_len = 0;
        decision.qc = "PASS_BOUNDARY_BLUNT";
        return decision;
    }

    if (link_delta <= tsd_max_len) {
        decision.pass = true;
        decision.boundary_type = "SMALL_DEL";
        decision.boundary_len = link_delta;
        decision.qc = "PASS_BOUNDARY_SMALL_DEL";
        return decision;
    }

    decision.qc = "REJECT_BOUNDARY_DEL_RANGE";
    return decision;
}

EventExistenceEvidence build_event_existence_evidence(
    const EventGenotypeInput& input) {
    EventExistenceEvidence evidence;
    const EventGenotypeDecision decision = genotype_event_from_alt_vs_ref(input);
    evidence.best_gt = decision.best_gt;
    evidence.af = decision.allele_fraction;
    evidence.gq = decision.gq;
    evidence.alt_struct_reads = std::max(0, input.alt_struct_reads);
    evidence.alt_split_reads = input.alt_split_reads;
    evidence.alt_indel_reads = input.alt_indel_reads;
    evidence.alt_left_clip_reads = input.alt_left_clip_reads;
    evidence.alt_right_clip_reads = input.alt_right_clip_reads;
    evidence.ref_span_reads = std::max(0, input.ref_span_reads);
    evidence.depth = decision.depth;
    evidence.best_nonref_minus_ref_ll = decision.best_nonref_minus_ref_ll;
    evidence.score = clamp_score((static_cast<double>(decision.gq) - 20.0) / 20.0);
    return evidence;
}

ExplanationDecision evaluate_event_explanations_for_test(
    const EventExistenceEvidence& existence,
    const EventSegmentationEvidence& segmentation,
    const TEAlignmentEvidence& te_alignment,
    const BoundaryEvidence& boundary,
    const ClipInsertConcordanceEvidence* /*clip_insert_concordance*/) {
    std::vector<EventExplanation> explanations;
    explanations.push_back(make_reference_explanation(existence, segmentation));
    explanations.push_back(make_non_te_explanation(existence, segmentation, te_alignment));
    explanations.push_back(make_te_explanation(existence, segmentation, te_alignment, boundary));
    explanations.push_back(make_artifact_explanation(existence, segmentation, te_alignment, boundary));
    return compare_event_explanations(
        explanations,
        has_closed_te_breakpoints(segmentation, boundary));
}

LatentMechanismEvidence evaluate_latent_mechanism_lfdr(
    const EventExistenceEvidence& existence,
    const EventSegmentationEvidence& segmentation,
    const TEAlignmentEvidence& te_alignment,
    const BoundaryEvidence& boundary,
    const ClipInsertConcordanceEvidence* clip_insert_concordance,
    double target_q) {
    LatentMechanismEvidence out;
    const LatentFeatureVector features = build_latent_feature_vector(
        existence,
        segmentation,
        te_alignment,
        boundary,
        clip_insert_concordance);
    out.family_activity_prior = features.family_activity_prior;

    std::array<double, 6> log_evidence{};
    for (size_t i = 0; i < kLatentMechanismStates.size(); ++i) {
        const std::string state = kLatentMechanismStates[i];
        log_evidence[i] =
            latent_state_log_prior(features, state) +
            latent_state_feature_log_likelihood(features, state);
    }

    const double denom = logsumexp_values(log_evidence);
    double best_log = -1e300;
    for (size_t i = 0; i < kLatentMechanismStates.size(); ++i) {
        const std::string state = kLatentMechanismStates[i];
        const double posterior = std::exp(log_evidence[i] - denom);
        if (is_te_latent_state(state)) {
            out.te_posterior += posterior;
        } else if (state == "non_te_insert") {
            out.non_te_posterior = posterior;
        } else if (state == "artifact_reference") {
            out.artifact_posterior = posterior;
        }
        if (log_evidence[i] > best_log) {
            best_log = log_evidence[i];
            out.latent_mechanism = state;
        }
    }
    out.lfdr = out.non_te_posterior + out.artifact_posterior;

    const double mechanistic_gap = 1.0 - features.independent_signal;
    const double te_ambiguity =
        0.95 + (2.60 * mechanistic_gap) + (0.75 * features.ref_conflict_signal);
    const double null_ambiguity =
        0.95 + (1.65 * mechanistic_gap) + (0.55 * features.ref_conflict_signal);
    double log_te_wc = -1e300;
    double log_null_wc = -1e300;
    for (size_t i = 0; i < kLatentMechanismStates.size(); ++i) {
        const std::string state = kLatentMechanismStates[i];
        if (is_te_latent_state(state)) {
            log_te_wc = logsumexp_pair(log_te_wc, log_evidence[i] - te_ambiguity);
        } else {
            log_null_wc = logsumexp_pair(log_null_wc, log_evidence[i] + null_ambiguity);
        }
    }
    out.worst_case_lfdr = std::exp(
        log_null_wc - logsumexp_pair(log_te_wc, log_null_wc));
    out.lfdr_qc =
        out.worst_case_lfdr <= std::clamp(target_q, 0.0, 1.0)
            ? "PASS_TE_LFDR"
            : "TE_LFDR_HIGH";
    return out;
}

// Structural-sanity veto for the TE-unknown hypothesis (h2): the hypothesis is
// ineligible for ranking only when the event definitionally cannot be a TE call
// -- no insert sequence, no TE alignment at all, or annotation quality too low to
// interpret. Whether an eligible TE-unknown hypothesis is actually emitted is
// decided downstream by the calibrated worst-case local FDR, not here.
// (`existence` / `boundary` are unused now that the read-count / one-sided
// ladders are gone, but are kept in the signature for symmetry with h3.)
bool compute_te_unknown_hard_veto(
    const EventExistenceEvidence& /*existence*/,
    const EventSegmentationEvidence& segmentation,
    const TEAlignmentEvidence& te_alignment,
    const BoundaryEvidence& /*boundary*/) {
    const TeAlignmentQc qc = classify_te_alignment_qc(te_alignment.qc_reason);
    return !segmentation.has_insert_seq ||
           qc == TeAlignmentQc::kNoTeAlignment ||
           !annotation_quality_allows_te_decision(te_alignment);
}

// Structural-sanity veto for the TE-resolved hypothesis (h3): as above, but a
// resolved call additionally requires a passing TE alignment with a resolved
// (non-"unknown") qc_reason. Emission remains gated by the worst-case local FDR.
bool compute_te_resolved_hard_veto(
    const EventExistenceEvidence& /*existence*/,
    const EventSegmentationEvidence& segmentation,
    const TEAlignmentEvidence& te_alignment,
    const BoundaryEvidence& /*boundary*/) {
    return !segmentation.has_insert_seq || !te_alignment.pass ||
           classify_te_alignment_qc(te_alignment.qc_reason) ==
               TeAlignmentQc::kPassInsertUnknown ||
           !annotation_quality_allows_te_decision(te_alignment);
}

// Copies the already-computed posterior / lFDR / mechanistic / structure
// signals into the diagnostic fields of the result. Pure plumbing, split out
// of evaluate_joint_hypotheses so the decision logic there stays legible.
void populate_joint_diagnostics(
    JointDecisionResult& result,
    const LocalHypothesisPosterior& posterior,
    const LatentMechanismEvidence& latent_lfdr,
    const MechanisticEvidenceCertificate& mechanistic_certificate,
    const RobustMechanisticLfdrResult& robust_mechanistic,
    const TeSequenceExplanation& te_structure) {
    result.te_posterior = posterior.te;
    result.non_te_posterior = posterior.non_te;
    result.artifact_posterior = posterior.artifact;
    result.te_vs_artifact_log_odds = posterior.te_vs_artifact_log_odds;
    result.te_vs_non_te_log_odds = posterior.te_vs_non_te_log_odds;
    result.posterior_qc = posterior.qc;
    result.latent_mechanism = latent_lfdr.latent_mechanism;
    result.family_activity_prior = latent_lfdr.family_activity_prior;
    result.lfdr = latent_lfdr.lfdr;
    result.worst_case_lfdr = latent_lfdr.worst_case_lfdr;
    result.lfdr_qc = latent_lfdr.lfdr_qc;
    result.mechanistic_lower_log_bf_te_vs_artifact =
        mechanistic_certificate.lower_log_bf_te_vs_artifact;
    result.mechanistic_lower_log_bf_te_vs_non_te =
        mechanistic_certificate.lower_log_bf_te_vs_non_te;
    result.mechanistic_ref_conflict_signal =
        mechanistic_certificate.ref_conflict_signal;
    result.mechanistic_ambiguity_width =
        mechanistic_certificate.ambiguity_width;
    result.mechanistic_blocks =
        serialize_mechanistic_blocks(mechanistic_certificate);
    result.robust_mechanistic_lfdr = robust_mechanistic.lfdr;
    result.robust_mechanistic_worst_case_lfdr =
        robust_mechanistic.worst_case_lfdr;
    result.robust_mechanistic_qc = robust_mechanistic.qc;
    result.te_structure_path = serialize_te_sequence_path(te_structure);
    result.te_structure_log_evidence = te_structure.te_structure_log_evidence;
    result.nonte_structure_log_evidence = te_structure.nonte_structure_log_evidence;
    result.artifact_structure_log_evidence =
        te_structure.artifact_structure_log_evidence;
    result.te_structure_path_confidence =
        te_structure.structure_path_confidence;
    result.polyA_posterior = te_structure.polyA_posterior;
    result.transduction_posterior = te_structure.transduction_posterior;
    result.te_core_coverage = te_structure.te_core_coverage;
    result.unexplained_high_complexity_bp =
        te_structure.unexplained_high_complexity_bp;
}

// L4: context-conditioned null. The prior over the TE / non-TE / artifact
// hypotheses is not universal -- it depends on the local context. In
// artifact-prone contexts (low-complexity / tandem inserts, one-sided or
// reference-conflicted junctions, edge/outlier sequence models -- summarised by
// the certificate's artifact_context_signal) the artifact/null prior is
// strengthened, so the worst-case lFDR gate demands more evidence there; in clean
// contexts it relaxes to the base prior. This adapts the null to context instead
// of applying one fixed prior everywhere, the decision-layer analogue of an
// empirical per-context background.
static PriorInterval context_conditioned_prior(double artifact_context_signal) {
    const double a = std::clamp(artifact_context_signal, 0.0, 1.0);
    PriorInterval prior;  // base (clean-context) prior
    // Lower the TE prior floor in artifact-prone contexts (raises worst-case
    // lFDR); leave it at the base 0.05 in clean contexts.
    prior.te_min = std::clamp(0.05 * (1.0 - (0.7 * a)), 5e-3, 0.05);
    prior.te_max = std::clamp(0.80 - (0.20 * a), 0.10, 0.80);
    prior.artifact_min = std::clamp(0.50 + (0.25 * a), 0.0, 0.95);
    prior.artifact_max = 0.95;
    prior.non_te_min = 1e-6;
    prior.non_te_max = 0.30;
    return prior;
}

JointDecisionResult evaluate_joint_hypotheses(
    const EventExistenceEvidence& existence,
    const EventSegmentationEvidence& segmentation,
    const TEAlignmentEvidence& te_alignment,
    const BoundaryEvidence& boundary,
    const ClipInsertConcordanceEvidence* clip_insert_concordance) {
    JointDecisionResult result;
    const double te_model_support = te_sequence_model_score(te_alignment);
    const double existence_support = adjusted_nonref_existence_score(
        existence,
        segmentation,
        te_model_support);
    const double segmentation_support = adjusted_segmentation_score(segmentation);
    const double te_boundary_support = adjusted_te_boundary_score(segmentation, boundary);
    const LocalHypothesisPosterior posterior =
        evaluate_local_hypothesis_posterior(
            existence,
            segmentation,
            te_alignment,
            boundary);
    const LatentMechanismEvidence latent_lfdr =
        evaluate_latent_mechanism_lfdr(
            existence,
            segmentation,
            te_alignment,
            boundary,
            clip_insert_concordance);
    const MechanisticEvidenceCertificate mechanistic_certificate =
        build_mechanistic_evidence_certificate(
            existence,
            segmentation,
            te_alignment,
            boundary,
            clip_insert_concordance);
    const TeSequenceExplanation te_structure =
        te_structure_explanation_for_decision(segmentation, te_alignment);
    const RobustMechanisticLfdrResult robust_mechanistic =
        evaluate_robust_mechanistic_lfdr(
            mechanistic_certificate,
            context_conditioned_prior(mechanistic_certificate.artifact_context_signal),
            DT::kTargetQ);

    JointHypothesisScore h0 = make_hypothesis(FinalHypothesisKind::kReference);
    h0.existence = DT::kScoreRefExistenceWeight * std::max(existence.score, 0.0);
    h0.segmentation = DT::kScoreRefSegmentationWeight * std::max(segmentation_support, 0.0);
    h0.te = 0.0;
    h0.total = h0.existence + h0.segmentation + h0.te;
    h0.reason = "REFERENCE";

    JointHypothesisScore h1 = make_hypothesis(FinalHypothesisKind::kInsertionNonTe);
    h1.existence = existence_support;
    h1.segmentation = DT::kScoreNonTeSegmentationWeight * segmentation_support;
    h1.te = DT::kScoreNonTeModelWeight * std::max(te_model_support, 0.0);
    h1.total = h1.existence + h1.segmentation + h1.te;
    h1.reason = "NON_TE_INSERTION";

    JointHypothesisScore h2 = make_hypothesis(FinalHypothesisKind::kTeUnknown);
    h2.existence = existence_support;
    h2.segmentation = DT::kScoreTeSegmentationWeight * segmentation_support;
    h2.te = te_model_support;
    h2.boundary = DT::kScoreTeBoundaryWeight * te_boundary_support;
    h2.total = h2.existence + h2.segmentation + h2.te + h2.boundary;
    h2.reason = "TE_UNKNOWN";
    h2.hard_veto = compute_te_unknown_hard_veto(
        existence, segmentation, te_alignment, boundary);

    JointHypothesisScore h3 = make_hypothesis(FinalHypothesisKind::kTeResolved);
    h3.existence = existence_support;
    h3.segmentation = DT::kScoreTeSegmentationWeight * segmentation_support;
    h3.te = te_model_support;
    h3.boundary = DT::kScoreTeBoundaryWeight * te_boundary_support;
    h3.total = h3.existence + h3.segmentation + h3.te + h3.boundary;
    h3.reason = "TE_RESOLVED";
    h3.hard_veto = compute_te_resolved_hard_veto(
        existence, segmentation, te_alignment, boundary);

    std::array<JointHypothesisScore, 4> all = {h0, h1, h2, h3};
    std::sort(all.begin(), all.end(), [](const JointHypothesisScore& lhs, const JointHypothesisScore& rhs) {
        const double lhs_total = lhs.hard_veto ? -1e9 : lhs.total;
        const double rhs_total = rhs.hard_veto ? -1e9 : rhs.total;
        if (lhs_total != rhs_total) {
            return lhs_total > rhs_total;
        }
        return static_cast<int>(lhs.kind) > static_cast<int>(rhs.kind);
    });

    result.best = all[0];
    result.runner_up = all[1];

    const ExplanationDecision explanation_decision =
        evaluate_event_explanations_for_test(
            existence,
            segmentation,
            te_alignment,
            boundary,
            clip_insert_concordance);

    result.explanation_decision = explanation_decision;
    result.best_explanation = explanation_kind_name(explanation_decision.best.kind);
    result.explanation_residual = serialize_residual(explanation_decision.best.residual);
    result.explanation_path = serialize_explanation_path(explanation_decision.best);
    // The explanation engine is an advisory hypothesis classifier: it proposes
    // whether the event is TE-like, non-TE, or artifact. It does not gate
    // emission.
    result.emit_te_call = explanation_decision.emit_te_call;
    result.emit_unknown_te = explanation_decision.emit_unknown_te;
    result.emit_evidence_te_call = explanation_decision.emit_evidence_te_call;
    result.final_qc = explanation_decision.final_qc;

    populate_joint_diagnostics(
        result,
        posterior,
        latent_lfdr,
        mechanistic_certificate,
        robust_mechanistic,
        te_structure);

    // PRIMARY EMISSION GATE: a final TE call is emitted only when the calibrated
    // robust worst-case local FDR is at or below the single target risk q, and
    // the event definitionally can carry TE sequence. There are no read-count /
    // insert-length / genotype-quality / segmentation-score ladders. A TE-like
    // event that fails the gate is not rejected as reference/artifact; it is
    // retained as TE-evidence (abstention), matching the precision-first design.
    const double q = DT::kTargetQ;
    const bool robust_te_lfdr_pass =
        result.robust_mechanistic_qc == "PASS_TE_LFDR" &&
        result.robust_mechanistic_worst_case_lfdr <= q;
    const bool structural_sanity =
        segmentation.has_insert_seq && has_te_sequence_evidence(te_alignment);
    if (result.emit_te_call && !(robust_te_lfdr_pass && structural_sanity)) {
        result.emit_te_call = false;
        result.emit_unknown_te = false;
        result.emit_evidence_te_call = has_te_sequence_evidence(te_alignment);
        result.final_qc = result.emit_evidence_te_call
            ? "TE_AMBIGUOUS"
            : "REFERENCE_OR_ARTIFACT";
    }

    // Non-TE structural insertion (h1): a separate closed-structural log-evidence
    // gate. Considered only when no confident TE call was emitted; it asserts
    // "there is a real non-reference insertion here" on decisive event-existence
    // evidence with a low artifact posterior, even when the TE identity is
    // ambiguous.
    if (!result.emit_te_call &&
        should_emit_structural_event_call(
            existence,
            segmentation,
            boundary,
            posterior,
            mechanistic_certificate,
            explanation_decision.best.kind == ExplanationKind::kInsertionNonTe)) {
        result.emit_structural_event_call = true;
        result.emit_unknown_te = true;
        result.emit_evidence_te_call = true;
        result.final_qc = "PASS_STRUCTURAL_INSERTION";
        result.best = h1;
        result.best.hard_veto = false;
        result.best.reason = "STRUCTURAL_INSERTION";
    }

    if (result.emit_te_call) {
        if (result.emit_unknown_te) {
            result.best = h2;
            result.best.hard_veto = false;
            result.best.reason = "TE_UNKNOWN";
        } else {
            result.best = h3;
            result.best.hard_veto = false;
            result.best.reason = "TE_RESOLVED";
        }
    }
    return result;
}

EventSegmentationEvidence analyze_event_segmentation_for_test(
    bool has_consensus,
    const EventSegmentation& segmentation) {
    EventSegmentationEvidence evidence;
    evidence.has_consensus = has_consensus;
    evidence.has_left_flank = segmentation.left_flank_align_len > 0;
    evidence.has_right_flank = segmentation.right_flank_align_len > 0;
    evidence.has_insert_seq = !segmentation.insert_seq.empty();
    evidence.pair_valid = segmentation.pass;
    evidence.left_align_len = segmentation.left_flank_align_len;
    evidence.right_align_len = segmentation.right_flank_align_len;
    evidence.left_identity = segmentation.left_flank_identity;
    evidence.right_identity = segmentation.right_flank_identity;
    evidence.insert_len = static_cast<int32_t>(segmentation.insert_seq.size());
    evidence.qc = segmentation.qc_reason;

    int32_t scoring_flank_len = std::min(evidence.left_align_len, evidence.right_align_len);
    double scoring_flank_identity = std::min(evidence.left_identity, evidence.right_identity);
    const bool one_sided_pass =
        evidence.pair_valid &&
        evidence.has_insert_seq &&
        (evidence.has_left_flank != evidence.has_right_flank);
    if (one_sided_pass) {
        scoring_flank_len =
            evidence.has_left_flank ? evidence.left_align_len : evidence.right_align_len;
        scoring_flank_identity =
            evidence.has_left_flank ? evidence.left_identity : evidence.right_identity;
    }
    evidence.score =
        0.5 * clamp_score((static_cast<double>(scoring_flank_len) - 50.0) / 25.0, -2.0, 2.0) +
        0.5 * clamp_score((scoring_flank_identity - 0.90) / 0.05, -2.0, 2.0);
    if (one_sided_pass) {
        evidence.score -= 0.5;
    } else if (!evidence.has_left_flank || !evidence.has_right_flank) {
        evidence.score -= 1.0;
    }
    if (!evidence.pair_valid) {
        evidence.score -= 1.5;
    }
    return evidence;
}

// Biological boundary-structure prior. Instead of a flat score for any boundary
// that clears the discrete range test, grade the evidence by how well the
// junction geometry matches a genuine TE-insertion signature: a target-site
// duplication has a characteristic length (broadly log-normal around ~12 bp for
// common families, so a 10 bp TSD is strong and a 55 bp "TSD" is not), a blunt
// join is clean, and a target-site deletion is credible but decays with size.
// Returns a bounded log-LR versus a featureless artifact junction; the discrete
// boundary_type label is preserved for reporting.
double boundary_structure_log_lr(const std::string& type, int32_t len) {
    if (type == "BLUNT") {
        return 1.0;
    }
    if (type == "TSD") {
        const double x = std::max(1, len);
        const double z = (std::log(x) - std::log(12.0)) / 0.7;
        return std::clamp(1.2 - (0.5 * z * z), -2.0, 1.2);
    }
    if (type == "SMALL_DEL") {
        return std::clamp(0.8 - (static_cast<double>(std::max(0, len)) / 80.0), 0.2, 0.8);
    }
    if (type == "NONCANONICAL") {
        return 0.25;
    }
    return -2.0;
}

BoundaryEvidence evaluate_boundary_evidence(
    const FinalBoundaryInput& input,
    int32_t breakpoint_envelope_width) {
    BoundaryEvidence evidence;
    const FinalBoundaryDecision canonical = check_boundary_consistency(input);
    evidence.geometry_defined =
        input.left_ref_start >= 0 && input.left_ref_end >= 0 &&
        input.right_ref_start >= 0 && input.right_ref_end >= 0 &&
        input.left_ref_start < input.left_ref_end &&
        input.right_ref_start < input.right_ref_end;
    evidence.canonical_pass = canonical.pass;
    evidence.boundary_type = canonical.boundary_type;
    evidence.boundary_len = canonical.boundary_len;
    evidence.qc = canonical.qc;

    if (!evidence.geometry_defined) {
        evidence.score = -2.0;
        return evidence;
    }

    if (canonical.pass) {
        evidence.evidence_consistent = true;
        evidence.score = boundary_structure_log_lr(
            canonical.boundary_type, canonical.boundary_len);
        return evidence;
    }

    const int32_t link_delta = input.right_ref_start - input.left_ref_end;
    const int32_t noncanonical_span = std::abs(link_delta);
    evidence.evidence_consistent =
        breakpoint_envelope_width > 0 && noncanonical_span <= breakpoint_envelope_width;
    if (evidence.evidence_consistent) {
        evidence.boundary_type = "NONCANONICAL";
        evidence.boundary_len = noncanonical_span;
        evidence.qc = "PASS_BOUNDARY_NONCANONICAL_CONSISTENT";
        evidence.score = boundary_structure_log_lr("NONCANONICAL", noncanonical_span);
    } else {
        evidence.score = -2.0;
    }
    return evidence;
}

}  // namespace placer
