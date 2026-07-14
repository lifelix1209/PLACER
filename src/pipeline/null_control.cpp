#include "null_control.h"

#include <algorithm>

namespace placer {

std::vector<BreakpointShiftControl> make_breakpoint_shift_controls(
    int32_t bp_left,
    int32_t bp_right,
    int32_t window_start,
    int32_t window_end,
    int32_t shift_step,
    int32_t max_controls) {
    std::vector<BreakpointShiftControl> controls;
    if (bp_left < 0 || bp_right < 0 || window_start > window_end ||
        shift_step <= 0 || max_controls <= 0) {
        return controls;
    }

    const int32_t width = std::max(0, bp_right - bp_left);
    for (int32_t multiplier = 1;
         static_cast<int32_t>(controls.size()) < max_controls;
         ++multiplier) {
        bool added = false;
        for (int32_t direction : {-1, 1}) {
            if (static_cast<int32_t>(controls.size()) >= max_controls) {
                break;
            }
            const int32_t shifted_left = bp_left + (direction * multiplier * shift_step);
            const int32_t shifted_right = shifted_left + width;
            if (shifted_left < window_start || shifted_right > window_end) {
                continue;
            }
            if (shifted_left == bp_left && shifted_right == bp_right) {
                continue;
            }
            BreakpointShiftControl control;
            control.bp_left = shifted_left;
            control.bp_right = shifted_right;
            controls.push_back(control);
            added = true;
        }
        if (!added &&
            bp_left - (multiplier * shift_step) < window_start &&
            bp_left + (multiplier * shift_step) + width > window_end) {
            break;
        }
    }
    return controls;
}

void EmpiricalNullTail::add(double value) {
    values_.push_back(value);
}

double EmpiricalNullTail::upper_tail_p(double observed) const {
    size_t ge = 0;
    for (double value : values_) {
        if (value >= observed) {
            ++ge;
        }
    }
    return static_cast<double>(1 + ge) /
           static_cast<double>(2 + values_.size());
}

size_t EmpiricalNullTail::size() const {
    return values_.size();
}

}  // namespace placer
