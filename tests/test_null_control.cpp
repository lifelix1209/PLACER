#ifdef NDEBUG
#undef NDEBUG
#endif

#include "null_control.h"

#include <cassert>

int main() {
    const auto shifts = placer::make_breakpoint_shift_controls(
        1000, 1012, 800, 1300, 50, 4);
    assert(shifts.size() == 4);
    for (const auto& control : shifts) {
        assert(control.bp_left >= 800);
        assert(control.bp_right <= 1300);
        assert(!(control.bp_left == 1000 && control.bp_right == 1012));
    }

    placer::EmpiricalNullTail tail;
    tail.add(0.1);
    tail.add(0.2);
    tail.add(0.8);
    assert(tail.upper_tail_p(0.9) < tail.upper_tail_p(0.15));
    assert(tail.upper_tail_p(0.9) > 0.0);
    assert(tail.size() == 3);
    return 0;
}
