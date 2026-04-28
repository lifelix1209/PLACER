#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <string>
#include <vector>

#include "parallel_executor_internal.h"

namespace {

void require(bool condition, const char* message) {
    if (!condition) {
        std::fprintf(
            stderr,
            "test_parallel_bounded_scheduler failed: %s\n",
            message);
        std::fflush(stderr);
        std::exit(1);
    }
}

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

    BoundedTaskQueue<int> close_queue(1);
    const bool close_pushed = close_queue.push(7);
    require(close_pushed, "closed queue first push failed");
    close_queue.close();
    const bool close_second_push = close_queue.push(8);
    require(!close_second_push, "closed queue accepted push");
    require(close_queue.size() == 1, "closed queue wrong retained size");
    require(close_queue.pop(value), "closed queue failed to pop retained value");
    require(value == 7, "closed queue retained wrong value");
    require(!close_queue.pop(value), "closed queue popped after drain");

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

    CostAwareTaskQueue closed_ready_queue(1);
    const bool closed_ready_first_push =
        closed_ready_queue.push(make_task(0, 3, 1, 1, 0));
    require(closed_ready_first_push, "closed ready queue first push failed");
    closed_ready_queue.close();
    const bool closed_ready_second_push =
        closed_ready_queue.push(make_task(0, 4, 1, 1, 0));
    require(!closed_ready_second_push, "closed ready queue accepted push");

    DeterministicBinReducer reducer;
    reducer.expect(0, 0);
    reducer.expect(0, 1);
    reducer.expect(0, 2);

    reducer.push_ready(0, 2, ExactBinTaskResult{0, 2, placer::PipelineResult{}, 0.30});
    reducer.push_ready(0, 0, ExactBinTaskResult{0, 0, placer::PipelineResult{}, 0.10});
    reducer.push_ready(0, 1, ExactBinTaskResult{0, 1, placer::PipelineResult{}, 0.20});

    std::vector<int> committed;
    reducer.drain([&](int32_t tid, int32_t bin_index, ExactBinTaskResult payload) {
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

    {
        DeterministicBinReducer move_reducer;
        move_reducer.expect(1, 10);
        move_reducer.expect(1, 11);

        PipelineResult first;
        first.final_pass_calls = 3;
        PipelineResult second;
        second.final_pass_calls = 5;
        move_reducer.push_ready(1, 11, ExactBinTaskResult{1, 11, std::move(second), 0.20});
        move_reducer.push_ready(1, 10, ExactBinTaskResult{1, 10, std::move(first), 0.10});

        std::vector<int32_t> moved_bins;
        int64_t moved_final_calls = 0;
        move_reducer.drain([&](int32_t tid, int32_t bin_index, ExactBinTaskResult payload) {
            assert(tid == 1);
            assert(payload.tid == 1);
            assert(payload.bin_index == bin_index);
            moved_bins.push_back(bin_index);
            moved_final_calls += payload.partial.final_pass_calls;
        });

        assert((moved_bins == std::vector<int32_t>{10, 11}));
        assert(moved_final_calls == 8);
        require(move_reducer.expected_count() == 2, "move reducer wrong expected count");
        require(move_reducer.committed_count() == 2, "move reducer wrong committed count");
        require(!move_reducer.has_pending_ready(), "move reducer kept pending ready payload");
    }

    {
        DeterministicBinReducer delayed_reducer;
        delayed_reducer.expect(2, 0);
        delayed_reducer.expect(2, 1);
        delayed_reducer.expect(2, 2);

        std::vector<int32_t> committed_bins;
        delayed_reducer.push_ready(
            2,
            2,
            ExactBinTaskResult{2, 2, PipelineResult{}, 0.30});
        delayed_reducer.drain([&](int32_t, int32_t bin_index, ExactBinTaskResult) {
            committed_bins.push_back(bin_index);
        });
        require(committed_bins.empty(), "delayed reducer committed past missing earlier bin");

        delayed_reducer.push_ready(
            2,
            0,
            ExactBinTaskResult{2, 0, PipelineResult{}, 0.10});
        delayed_reducer.drain([&](int32_t, int32_t bin_index, ExactBinTaskResult) {
            committed_bins.push_back(bin_index);
        });
        require(
            committed_bins == std::vector<int32_t>{0},
            "delayed reducer committed non-contiguous result");

        delayed_reducer.push_ready(
            2,
            1,
            ExactBinTaskResult{2, 1, PipelineResult{}, 0.20});
        delayed_reducer.drain([&](int32_t, int32_t bin_index, ExactBinTaskResult) {
            committed_bins.push_back(bin_index);
        });
        require(
            committed_bins == std::vector<int32_t>{0, 1, 2},
            "delayed reducer did not commit contiguous ready suffix");
    }

    {
        CostAwareTaskQueue pressure_queue(1);
        require(
            pressure_queue.push(make_task(3, 0, 10, 5, 5)),
            "pressure queue initial push failed");
        require(pressure_queue.size() == 1, "pressure queue wrong size");
        require(
            !pressure_queue.try_push(make_task(3, 1, 9, 5, 4)),
            "pressure queue accepted task past capacity");
        ExactBinTask popped_task;
        require(pressure_queue.pop(popped_task), "pressure queue failed to pop");
        require(popped_task.bin_index == 0, "pressure queue popped wrong task");
        require(
            pressure_queue.push(make_task(3, 1, 9, 5, 4)),
            "pressure queue second push failed");
        pressure_queue.close();
        require(
            !pressure_queue.push(make_task(3, 2, 8, 4, 4)),
            "pressure queue accepted push after close");
    }
    return 0;
}
