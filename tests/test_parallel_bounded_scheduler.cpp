#include <cassert>
#include <cstdint>
#include <string>
#include <vector>

#include "parallel_executor_internal.h"

namespace {

placer::ExactBinTask make_task(
    int32_t tid,
    int32_t bin_index,
    size_t merged_input_records,
    size_t current_bin_records,
    size_t context_bin_records) {
    placer::ExactBinTask task;
    task.tid = tid;
    task.bin_index = bin_index;
    task.merged_input_records = merged_input_records;
    task.current_bin_records = current_bin_records;
    task.context_bin_records = context_bin_records;
    return task;
}

}  // namespace

int main() {
    using namespace placer;

    BoundedTaskQueue<int> snapshot_queue(2);
    assert(snapshot_queue.capacity() == 2);
    assert(snapshot_queue.try_push(1));
    assert(snapshot_queue.try_push(2));
    assert(!snapshot_queue.try_push(3));

    int value = 0;
    assert(snapshot_queue.pop(value));
    assert(value == 1);
    assert(snapshot_queue.try_push(3));

    CostAwareTaskQueue ready_queue(3);
    assert(ready_queue.try_push(make_task(0, 0, 4, 2, 2)));
    assert(ready_queue.try_push(make_task(0, 1, 11, 6, 5)));
    assert(ready_queue.try_push(make_task(0, 2, 7, 4, 3)));

    ExactBinTask popped;
    assert(ready_queue.pop(popped));
    assert(popped.bin_index == 1);
    assert(ready_queue.pop(popped));
    assert(popped.bin_index == 2);
    assert(ready_queue.pop(popped));
    assert(popped.bin_index == 0);

    DeterministicBinReducer reducer;
    reducer.expect(0, 0);
    reducer.expect(0, 1);
    reducer.expect(0, 2);

    reducer.push_ready(0, 2, ExactBinTaskResult{0, 2, placer::PipelineResult{}, 0.30});
    reducer.push_ready(0, 0, ExactBinTaskResult{0, 0, placer::PipelineResult{}, 0.10});
    reducer.push_ready(0, 1, ExactBinTaskResult{0, 1, placer::PipelineResult{}, 0.20});

    std::vector<int> committed;
    reducer.drain([&](int32_t tid, int32_t bin_index, const ExactBinTaskResult& payload) {
        assert(tid == 0);
        committed.push_back(bin_index);
        assert(payload.tid == 0);
        assert(payload.bin_index == bin_index);
    });

    assert((committed == std::vector<int>{0, 1, 2}));

    ParallelProgressSnapshot progress;
    progress.reads_scanned = 120000;
    progress.raw_bins_discovered = 3;
    progress.snapshot_window_size = 2;
    progress.tasks_running = 2;
    progress.tasks_completed = 7;
    progress.reducer_committed = 5;
    progress.ready_queue_depth = 4;
    progress.result_queue_depth = 1;

    const std::string text = render_parallel_progress(progress);
    assert(text.find("reads_scanned=120000") != std::string::npos);
    assert(text.find("raw_bins_discovered=3") != std::string::npos);
    assert(text.find("tasks_completed=7") != std::string::npos);
    assert(text.find("snapshot_window_size=2") != std::string::npos);
    assert(text.find("ready_queue_depth=4") != std::string::npos);
    return 0;
}
