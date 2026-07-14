#ifdef NDEBUG
#undef NDEBUG
#endif

#include "conformal_selector.h"

#include <cassert>
#include <vector>

namespace {

placer::ConformalFeatureVector make_feature(
    const char* id,
    std::vector<double> pro_te,
    double ref_span_reads) {
    placer::ConformalFeatureVector feature;
    feature.id = id;
    feature.pro_te = std::move(pro_te);
    feature.ref_span_reads = ref_span_reads;
    return feature;
}

void dominance_p_value_is_weight_free_and_monotone() {
    placer::ConformalNullSelector selector;
    selector.add_null_control(make_feature("n1", {4.0, 0.70, 0.60, 0.02}, 8.0));
    selector.add_null_control(make_feature("n2", {6.0, 0.82, 0.70, 0.05}, 4.0));
    selector.add_null_control(make_feature("n3", {9.0, 0.90, 0.82, 0.10}, 1.0));

    const auto weak = make_feature("weak", {5.0, 0.75, 0.65, 0.04}, 6.0);
    const auto strong = make_feature("strong", {12.0, 0.96, 0.93, 0.18}, 0.0);

    const double weak_p = selector.dominance_p_value(weak);
    const double strong_p = selector.dominance_p_value(strong);
    assert(strong_p < weak_p);
    assert(strong_p > 0.0);
}

void by_selector_keeps_only_sample_null_significant_candidates() {
    placer::ConformalNullSelector selector;
    for (int i = 0; i < 20; ++i) {
        selector.add_null_control(make_feature(
            ("n" + std::to_string(i)).c_str(),
            {static_cast<double>(i % 5), 0.60, 0.50, 0.01},
            5.0));
    }

    std::vector<placer::ConformalFeatureVector> candidates;
    candidates.push_back(make_feature("strong_a", {20.0, 0.98, 0.95, 0.20}, 0.0));
    candidates.push_back(make_feature("strong_b", {18.0, 0.97, 0.94, 0.19}, 0.0));
    candidates.push_back(make_feature("weak", {3.0, 0.60, 0.50, 0.01}, 6.0));

    const auto results = selector.select(candidates, 0.20);
    assert(results.size() == 3);
    assert(results[0].id == "strong_a");
    assert(results[0].qc == "PASS_CONFORMAL_FDR");
    assert(results[1].qc == "PASS_CONFORMAL_FDR");
    assert(results[2].qc == "CONFORMAL_FDR_REJECT");
    assert(results[0].by_threshold > 0.0);
}

void selection_requires_structural_support_to_be_null_unusual() {
    placer::ConformalNullSelector selector;
    for (int i = 0; i < 30; ++i) {
        selector.add_null_control(make_feature(
            ("n" + std::to_string(i)).c_str(),
            {10.0, 0.40, 0.40, 0.01},
            5.0));
    }

    std::vector<placer::ConformalFeatureVector> candidates;
    candidates.push_back(make_feature("sequence_only", {4.0, 0.99, 0.99, 0.80}, 0.0));
    candidates.push_back(make_feature("structural_and_sequence", {20.0, 0.99, 0.99, 0.80}, 0.0));

    const auto results = selector.select(candidates, 0.20);
    assert(results.size() == 2);
    assert(results[0].id == "sequence_only");
    assert(results[0].qc == "CONFORMAL_FDR_REJECT");
    assert(results[1].qc == "PASS_CONFORMAL_FDR");
    assert(results[0].conformal_p > results[1].conformal_p);
}

void empty_null_controls_reject_all_candidates() {
    placer::ConformalNullSelector selector;
    std::vector<placer::ConformalFeatureVector> candidates;
    candidates.push_back(make_feature("candidate", {20.0, 0.98, 0.95, 0.20}, 0.0));

    const auto results = selector.select(candidates, 0.10);
    assert(results.size() == 1);
    assert(results[0].conformal_p == 1.0);
    assert(results[0].qc == "CONFORMAL_NULL_INSUFFICIENT");
    assert(!results[0].pass);
}

}  // namespace

int main() {
    dominance_p_value_is_weight_free_and_monotone();
    by_selector_keeps_only_sample_null_significant_candidates();
    selection_requires_structural_support_to_be_null_unusual();
    empty_null_controls_reject_all_candidates();
    return 0;
}
