#include <cassert>
#include <string>
#include <vector>

#include "local_interval_cache.h"

int main() {
    using namespace placer;

    std::vector<LocalIntervalRequest> requests = {
        LocalIntervalRequest{"chr1", 1000, 2100, 0},
        LocalIntervalRequest{"chr1", 1500, 2600, 1},
        LocalIntervalRequest{"chr1", 8000, 9000, 2},
    };

    const std::vector<CanonicalLocalInterval> intervals =
        build_canonical_local_intervals(requests, 128);

    assert(intervals.size() == 2);
    assert(intervals[0].chrom == "chr1");
    assert(intervals[0].start == 1000);
    assert(intervals[0].end == 2600);
    assert((intervals[0].request_ids == std::vector<size_t>{0, 1}));
    assert(intervals[1].start == 8000);
    assert(intervals[1].end == 9000);
    assert((intervals[1].request_ids == std::vector<size_t>{2}));
    return 0;
}
