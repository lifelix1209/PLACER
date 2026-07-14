#ifndef PLACER_NULL_CONTROL_H
#define PLACER_NULL_CONTROL_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace placer {

struct BreakpointShiftControl {
    int32_t bp_left = -1;
    int32_t bp_right = -1;
    std::string kind = "breakpoint_shift";
};

std::vector<BreakpointShiftControl> make_breakpoint_shift_controls(
    int32_t bp_left,
    int32_t bp_right,
    int32_t window_start,
    int32_t window_end,
    int32_t shift_step,
    int32_t max_controls);

class EmpiricalNullTail {
public:
    void add(double value);
    double upper_tail_p(double observed) const;
    size_t size() const;

private:
    std::vector<double> values_;
};

}  // namespace placer

#endif  // PLACER_NULL_CONTROL_H
