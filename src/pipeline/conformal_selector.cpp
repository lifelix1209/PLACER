#include "conformal_selector.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <utility>

namespace placer {
namespace {

bool finite_feature_vector(const ConformalFeatureVector& feature) {
    if (!std::isfinite(feature.ref_span_reads)) {
        return false;
    }
    for (double value : feature.pro_te) {
        if (!std::isfinite(value)) {
            return false;
        }
    }
    return true;
}

bool null_dominates_candidate(
    const ConformalFeatureVector& null_feature,
    const ConformalFeatureVector& candidate) {
    if (null_feature.pro_te.empty() ||
        null_feature.pro_te.size() != candidate.pro_te.size() ||
        !finite_feature_vector(null_feature) ||
        !finite_feature_vector(candidate)) {
        return false;
    }
    for (size_t i = 0; i < candidate.pro_te.size(); ++i) {
        if (null_feature.pro_te[i] < candidate.pro_te[i]) {
            return false;
        }
    }
    return null_feature.ref_span_reads <= candidate.ref_span_reads;
}

double harmonic_number(size_t n) {
    double value = 0.0;
    for (size_t i = 1; i <= n; ++i) {
        value += 1.0 / static_cast<double>(i);
    }
    return value;
}

// L4: the null controls relevant to a candidate are those in the same local
// context. When a context has too few controls to calibrate against, fall back
// to the full pooled null so power is not lost. Deterministic.
constexpr size_t kMinContextNulls = 12;

std::vector<const ConformalFeatureVector*> context_relevant_nulls(
    const std::vector<ConformalFeatureVector>& null_controls,
    const ConformalFeatureVector& candidate) {
    std::vector<const ConformalFeatureVector*> same_context;
    for (const auto& null_feature : null_controls) {
        if (null_feature.context == candidate.context) {
            same_context.push_back(&null_feature);
        }
    }
    if (same_context.size() >= kMinContextNulls) {
        return same_context;
    }
    std::vector<const ConformalFeatureVector*> pooled;
    pooled.reserve(null_controls.size());
    for (const auto& null_feature : null_controls) {
        pooled.push_back(&null_feature);
    }
    return pooled;
}

double structural_support_p_value(
    const std::vector<const ConformalFeatureVector*>& null_controls,
    const ConformalFeatureVector& candidate) {
    if (null_controls.empty() ||
        candidate.pro_te.empty() ||
        !finite_feature_vector(candidate)) {
        return 1.0;
    }
    size_t ge = 0;
    for (const auto* null_feature : null_controls) {
        if (null_feature->pro_te.empty() || !finite_feature_vector(*null_feature)) {
            continue;
        }
        if (null_feature->pro_te.front() >= candidate.pro_te.front()) {
            ++ge;
        }
    }
    return static_cast<double>(1 + ge) /
           static_cast<double>(1 + null_controls.size());
}

}  // namespace

void ConformalNullSelector::add_null_control(const ConformalFeatureVector& feature) {
    if (!feature.pro_te.empty() && finite_feature_vector(feature)) {
        null_controls_.push_back(feature);
    }
}

size_t ConformalNullSelector::null_count() const {
    return null_controls_.size();
}

size_t ConformalNullSelector::dominated_null_count(
    const ConformalFeatureVector& candidate) const {
    const auto relevant = context_relevant_nulls(null_controls_, candidate);
    size_t count = 0;
    for (const auto* null_feature : relevant) {
        if (null_dominates_candidate(*null_feature, candidate)) {
            ++count;
        }
    }
    return count;
}

double ConformalNullSelector::dominance_p_value(
    const ConformalFeatureVector& candidate) const {
    if (null_controls_.empty() ||
        candidate.pro_te.empty() ||
        !finite_feature_vector(candidate)) {
        return 1.0;
    }
    const auto relevant = context_relevant_nulls(null_controls_, candidate);
    const size_t dominated = dominated_null_count(candidate);
    return static_cast<double>(1 + dominated) /
           static_cast<double>(1 + relevant.size());
}

std::vector<ConformalSelectionResult> ConformalNullSelector::select(
    const std::vector<ConformalFeatureVector>& candidates,
    double target_fdr) const {
    std::vector<ConformalSelectionResult> results;
    results.reserve(candidates.size());
    for (const auto& candidate : candidates) {
        const auto relevant = context_relevant_nulls(null_controls_, candidate);
        ConformalSelectionResult result;
        result.id = candidate.id;
        result.null_count = relevant.size();
        result.dominated_null_count = dominated_null_count(candidate);
        result.conformal_p = std::max(
            dominance_p_value(candidate),
            structural_support_p_value(relevant, candidate));
        result.qc = null_controls_.empty()
            ? "CONFORMAL_NULL_INSUFFICIENT"
            : "CONFORMAL_FDR_REJECT";
        results.push_back(std::move(result));
    }

    if (results.empty() || null_controls_.empty()) {
        return results;
    }

    std::vector<size_t> order(results.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&results](size_t lhs, size_t rhs) {
        if (results[lhs].conformal_p != results[rhs].conformal_p) {
            return results[lhs].conformal_p < results[rhs].conformal_p;
        }
        return results[lhs].id < results[rhs].id;
    });

    const double q = std::clamp(target_fdr, 0.0, 1.0);
    const double harmonic = std::max(1.0, harmonic_number(results.size()));
    size_t selected_prefix = 0;
    double selected_threshold = 0.0;
    for (size_t rank = 1; rank <= order.size(); ++rank) {
        const double threshold =
            (static_cast<double>(rank) * q) /
            (static_cast<double>(order.size()) * harmonic);
        const size_t index = order[rank - 1];
        if (results[index].conformal_p <= threshold) {
            selected_prefix = rank;
            selected_threshold = threshold;
        }
    }

    for (size_t rank = 1; rank <= order.size(); ++rank) {
        const size_t index = order[rank - 1];
        const double threshold =
            (static_cast<double>(rank) * q) /
            (static_cast<double>(order.size()) * harmonic);
        results[index].by_threshold = threshold;
        if (selected_prefix > 0 && rank <= selected_prefix) {
            results[index].pass = true;
            results[index].by_threshold = selected_threshold;
            results[index].qc = "PASS_CONFORMAL_FDR";
        }
    }
    return results;
}

}  // namespace placer
