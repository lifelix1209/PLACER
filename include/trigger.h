#ifndef PLACER_TRIGGER_H
#define PLACER_TRIGGER_H

#include "window_buffer.h"
#include <string>
#include <functional>
#include <atomic>
#include <algorithm>

namespace placer {

/**
 * TriggerResult: Result of a trigger evaluation
 */
struct TriggerResult {
    std::string window_id;
    bool triggered = false;
    double score = 0.0;

    // Debug info: which conditions triggered
    bool clip_trig = false;
    bool sa_trig = false;
    bool ins_trig = false;
};

/**
 * TriggerConfig: Configuration for trigger thresholds
 */
struct TriggerConfig {
    // Relative threshold multipliers (adaptive thresholding)
    double clip_threshold_factor = 2.0;
    double sa_threshold_factor = 2.0;
    double ins_threshold_factor = 2.0;

    // Absolute baseline thresholds (prevent division by zero and false positives)
    // Even if Q99 is 0, require biologically meaningful signals
    double min_baseline_clip = 50.0;      // At least 50bp of clip
    double min_baseline_sa = 1.0;         // At least 2 SA reads (> 1.0)
    double min_baseline_ins = 2.0;        // At least 2 insertion events

    // Score weights
    double weight_clip = 0.4;
    double weight_sa = 0.4;
    double weight_ins = 0.2;
};

/**
 * Trigger: Threshold-based triggering logic
 * Evaluates window statistics against quantile thresholds with robust division handling
 */
class Trigger {
public:
    explicit Trigger(TriggerConfig config = TriggerConfig());

    /**
     * Evaluate if a window should be triggered
     * Single call handles both condition checking and score calculation
     */
    TriggerResult evaluate(const std::string& window_id, const WindowStats& stats) const;

    TriggerConfig get_config() const { return config_; }

private:
    TriggerConfig config_;
};

/**
 * TriggerCallback: Callback for triggered windows
 */
using TriggerCallback = std::function<void(const TriggerResult&)>;

/**
 * TriggerEngine: Orchestrates triggering across windows
 */
class TriggerEngine {
public:
    explicit TriggerEngine(Trigger trigger, TriggerCallback callback);

    /**
     * Process a window buffer
     * Only calls evaluate once - no redundant condition checks
     */
    void process_window(WindowID window_id, const WindowBuffer* buffer);

    /**
     * Set the callback
     */
    void set_callback(TriggerCallback callback);

    /**
     * Get statistics
     */
    int64_t get_total_processed() const { return total_processed_.load(); }
    int64_t get_total_triggered() const { return total_triggered_.load(); }

private:
    Trigger trigger_;
    TriggerCallback callback_;
    std::atomic<int64_t> total_triggered_{0};
    std::atomic<int64_t> total_processed_{0};
};

}  // namespace placer

#endif  // PLACER_TRIGGER_H
