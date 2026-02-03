#include "trigger.h"
#include <cmath>
#include <algorithm>

namespace placer {

Trigger::Trigger(TriggerConfig config) : config_(std::move(config)) {}

TriggerResult Trigger::evaluate(const std::string& window_id, const WindowStats& stats) const {
    TriggerResult result;
    result.window_id = window_id;

    // 1. Calculate effective thresholds (max of adaptive and baseline)
    // This solves:
    //   a) Division by zero (baseline guarantees denominator > 0)
    //   b) Cold start / low-noise regions (prevents triggering on 1bp clip)
    double thresh_clip = std::max(stats.q99_clip * config_.clip_threshold_factor,
                                  config_.min_baseline_clip);
    double thresh_sa = std::max(stats.q99_sa * config_.sa_threshold_factor,
                                 config_.min_baseline_sa);
    double thresh_ins = std::max(stats.q99_ins * config_.ins_threshold_factor,
                                  config_.min_baseline_ins);

    // 2. Check trigger conditions
    if (static_cast<double>(stats.clip_bp) > thresh_clip) result.clip_trig = true;
    if (static_cast<double>(stats.sa_reads) > thresh_sa) result.sa_trig = true;
    if (static_cast<double>(stats.ins_events) > thresh_ins) result.ins_trig = true;

    result.triggered = result.clip_trig || result.sa_trig || result.ins_trig;

    // 3. Calculate normalized score [0, 1]
    // Using effective thresholds as normalizers ensures:
    //   - Safe division (thresh > 0 guaranteed)
    // - Score reflects relative intensity above baseline
    double score_clip = std::min(1.0, static_cast<double>(stats.clip_bp) / thresh_clip);
    double score_sa = std::min(1.0, static_cast<double>(stats.sa_reads) / thresh_sa);
    double score_ins = std::min(1.0, static_cast<double>(stats.ins_events) / thresh_ins);

    result.score = config_.weight_clip * score_clip +
                   config_.weight_sa * score_sa +
                   config_.weight_ins * score_ins;

    return result;
}

TriggerEngine::TriggerEngine(Trigger trigger, TriggerCallback callback)
    : trigger_(std::move(trigger)), callback_(std::move(callback)) {}

void TriggerEngine::process_window(WindowID window_id, const WindowBuffer* buffer) {
    total_processed_++;

    const WindowStats& stats = buffer->get_stats(window_id);

    // Single evaluate call - no redundant condition checks
    TriggerResult result = trigger_.evaluate("", stats);

    if (result.triggered) {
        total_triggered_++;
        if (callback_) {
            callback_(result);
        }
    }
}

void TriggerEngine::set_callback(TriggerCallback callback) {
    callback_ = std::move(callback);
}

}  // namespace placer
