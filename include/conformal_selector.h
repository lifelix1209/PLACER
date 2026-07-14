#ifndef PLACER_CONFORMAL_SELECTOR_H
#define PLACER_CONFORMAL_SELECTOR_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace placer {

struct ConformalFeatureVector {
    std::string id = "NA";
    std::vector<double> pro_te;
    double ref_span_reads = 0.0;
    // L4: discrete local-context bucket (e.g. clean / low-complexity-tandem /
    // one-sided). A candidate is compared against null controls in the same
    // context, so the empirical null is context-conditioned; the selector falls
    // back to the pooled null when a context has too few controls.
    int32_t context = 0;
};

struct ConformalSelectionResult {
    std::string id = "NA";
    double conformal_p = 1.0;
    double by_threshold = 0.0;
    size_t dominated_null_count = 0;
    size_t null_count = 0;
    bool pass = false;
    std::string qc = "CONFORMAL_NULL_INSUFFICIENT";
};

class ConformalNullSelector {
public:
    void add_null_control(const ConformalFeatureVector& feature);
    size_t null_count() const;

    double dominance_p_value(const ConformalFeatureVector& candidate) const;
    size_t dominated_null_count(const ConformalFeatureVector& candidate) const;

    std::vector<ConformalSelectionResult> select(
        const std::vector<ConformalFeatureVector>& candidates,
        double target_fdr) const;

private:
    std::vector<ConformalFeatureVector> null_controls_;
};

}  // namespace placer

#endif  // PLACER_CONFORMAL_SELECTOR_H
