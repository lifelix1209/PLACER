#include "parallel_executor_internal.h"

#include <cstdlib>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace placer {
namespace {

void append_exact_task_records(
    const std::vector<BufferedRecord>& src,
    std::vector<ExactTaskRecordRef>& dst,
    int32_t min_ref_end,
    bool from_current_bin,
    size_t& current_count,
    size_t& context_count) {
    for (const auto& buffered : src) {
        if (!buffered.record) {
            continue;
        }
        if (buffered.ref_end < min_ref_end) {
            continue;
        }
        dst.push_back(ExactTaskRecordRef{buffered.record.get(), from_current_bin});
        if (from_current_bin) {
            ++current_count;
        } else {
            ++context_count;
        }
    }
}

void prune_consumed_exact_bin_snapshots(
    std::deque<ExactBinSnapshot>& retained_snapshots,
    size_t& next_owner_offset,
    int32_t cross_bin_context_bins) {
    while (!retained_snapshots.empty() && next_owner_offset > 0) {
        bool can_drop_front = false;
        if (next_owner_offset < retained_snapshots.size()) {
            const int32_t next_owner_bin = retained_snapshots[next_owner_offset].bin_index;
            can_drop_front =
                retained_snapshots.front().bin_index <
                (next_owner_bin - cross_bin_context_bins);
        } else {
            can_drop_front =
                retained_snapshots.size() > static_cast<size_t>(cross_bin_context_bins);
        }

        if (!can_drop_front) {
            break;
        }

        retained_snapshots.pop_front();
        --next_owner_offset;
    }
}

void snapshot_current_exact_bin(ExactBinScanState& state) {
    if (state.current_tid < 0 || state.current_bin_index < 0 ||
        state.current_bin_records.empty()) {
        return;
    }

    SharedRecordBatch batch = std::make_shared<const std::vector<BufferedRecord>>(
        std::move(state.current_bin_records));
    state.current_bin_records.clear();
    state.recent_snapshots.push_back(ExactBinSnapshot{
        state.current_tid,
        state.current_bin_index,
        std::move(batch),
    });
}

std::vector<ExactBinTask> materialize_scan_ready_tasks(
    ExactBinScanState& state,
    bool flush_all,
    int32_t latest_completed_bin_index,
    int32_t bin_size,
    int32_t cross_bin_context_bins,
    int32_t window_bin_slack_bp) {
    snapshot_current_exact_bin(state);
    if (state.recent_snapshots.empty()) {
        if (flush_all) {
            state.current_tid = -1;
            state.current_bin_index = -1;
            state.next_owner_offset = 0;
        }
        return {};
    }

    const int32_t latest_completed =
        latest_completed_bin_index >= 0
            ? latest_completed_bin_index
            : state.recent_snapshots.back().bin_index;
    std::vector<ExactBinTask> ready = materialize_ready_exact_bin_tasks(
        state.recent_snapshots,
        state.next_owner_offset,
        latest_completed,
        flush_all,
        bin_size,
        cross_bin_context_bins,
        window_bin_slack_bp);

    if (flush_all) {
        state.recent_snapshots.clear();
        state.next_owner_offset = 0;
        state.current_tid = -1;
        state.current_bin_index = -1;
    }
    return ready;
}

}  // namespace

ExactBinTaskCost estimate_exact_bin_task_cost(const ExactBinTask& task) {
    ExactBinTaskCost cost;
    cost.tid = task.tid;
    cost.bin_index = task.bin_index;
    cost.merged_input_records = task.merged_input_records;
    cost.current_bin_records = task.current_bin_records;
    cost.context_bin_records = task.context_bin_records;
    return cost;
}

ExactBinTask build_exact_bin_task(
    int32_t tid,
    int32_t bin_index,
    const std::vector<ExactBinSnapshot>& snapshots,
    int32_t bin_size,
    int32_t cross_bin_context_bins,
    int32_t window_bin_slack_bp) {
    ExactBinTask task;
    task.tid = tid;
    task.bin_index = bin_index;
    task.bin_start = bin_index * bin_size;
    task.bin_end = task.bin_start + bin_size;

    const int32_t snapshot_min_ref_end = task.bin_start - window_bin_slack_bp;
    for (const auto& snapshot : snapshots) {
        if (snapshot.tid != tid || !snapshot.records) {
            continue;
        }
        if (std::abs(snapshot.bin_index - bin_index) > cross_bin_context_bins) {
            continue;
        }

        const bool is_current_bin = snapshot.bin_index == bin_index;
        task.owners.push_back(snapshot.records);
        append_exact_task_records(
            *snapshot.records,
            task.records,
            is_current_bin ? std::numeric_limits<int32_t>::min() : snapshot_min_ref_end,
            is_current_bin,
            task.current_bin_records,
            task.context_bin_records);
    }

    task.merged_input_records = task.records.size();
    if (task.current_bin_records == 0) {
        throw std::runtime_error("ExactBinTask materialized without current-bin records");
    }
    return task;
}

std::vector<ExactBinTask> materialize_ready_exact_bin_tasks(
    std::deque<ExactBinSnapshot>& retained_snapshots,
    size_t& next_owner_offset,
    int32_t latest_completed_bin_index,
    bool flush_all,
    int32_t bin_size,
    int32_t cross_bin_context_bins,
    int32_t window_bin_slack_bp) {
    std::vector<ExactBinTask> tasks;

    while (next_owner_offset < retained_snapshots.size()) {
        const ExactBinSnapshot& owner = retained_snapshots[next_owner_offset];
        if (!flush_all &&
            latest_completed_bin_index < (owner.bin_index + cross_bin_context_bins)) {
            break;
        }

        std::vector<ExactBinSnapshot> context_snapshots;
        context_snapshots.reserve(retained_snapshots.size());
        for (const auto& snapshot : retained_snapshots) {
            if (snapshot.tid != owner.tid || !snapshot.records) {
                continue;
            }
            if (std::abs(snapshot.bin_index - owner.bin_index) > cross_bin_context_bins) {
                continue;
            }
            context_snapshots.push_back(snapshot);
        }

        tasks.push_back(build_exact_bin_task(
            owner.tid,
            owner.bin_index,
            context_snapshots,
            bin_size,
            cross_bin_context_bins,
            window_bin_slack_bp));
        ++next_owner_offset;
        prune_consumed_exact_bin_snapshots(
            retained_snapshots,
            next_owner_offset,
            cross_bin_context_bins);
    }

    return tasks;
}

std::vector<ExactBinTask> append_record_to_exact_bin_scan(
    ExactBinScanState& state,
    BufferedRecord&& record,
    int32_t tid,
    int32_t bin_index,
    int32_t bin_size,
    int32_t cross_bin_context_bins,
    int32_t window_bin_slack_bp) {
    if (tid < 0 || bin_index < 0) {
        return {};
    }

    if (state.current_tid < 0) {
        state.current_tid = tid;
        state.current_bin_index = bin_index;
        state.current_bin_records.push_back(std::move(record));
        return {};
    }

    std::vector<ExactBinTask> ready;
    if (tid != state.current_tid) {
        ready = materialize_scan_ready_tasks(
            state,
            true,
            state.current_bin_index,
            bin_size,
            cross_bin_context_bins,
            window_bin_slack_bp);
        state.current_tid = tid;
        state.current_bin_index = bin_index;
    } else if (bin_index != state.current_bin_index) {
        const int32_t latest_completed_bin_index =
            bin_index > state.current_bin_index
                ? bin_index - 1
                : state.current_bin_index;
        ready = materialize_scan_ready_tasks(
            state,
            false,
            latest_completed_bin_index,
            bin_size,
            cross_bin_context_bins,
            window_bin_slack_bp);
        state.current_tid = tid;
        state.current_bin_index = bin_index;
    }

    state.current_bin_records.push_back(std::move(record));
    return ready;
}

std::vector<ExactBinTask> flush_exact_bin_scan(
    ExactBinScanState& state,
    int32_t bin_size,
    int32_t cross_bin_context_bins,
    int32_t window_bin_slack_bp) {
    return materialize_scan_ready_tasks(
        state,
        true,
        state.current_bin_index,
        bin_size,
        cross_bin_context_bins,
        window_bin_slack_bp);
}

CostAwareTaskQueue::CostAwareTaskQueue(size_t capacity)
    : capacity_(capacity) {}

bool CostAwareTaskQueue::ReadyEntryLowerPriority::operator()(
    const ReadyEntry& lhs,
    const ReadyEntry& rhs) const {
    if (lhs.cost.merged_input_records != rhs.cost.merged_input_records) {
        return lhs.cost.merged_input_records < rhs.cost.merged_input_records;
    }
    if (lhs.cost.current_bin_records != rhs.cost.current_bin_records) {
        return lhs.cost.current_bin_records < rhs.cost.current_bin_records;
    }
    if (lhs.cost.context_bin_records != rhs.cost.context_bin_records) {
        return lhs.cost.context_bin_records < rhs.cost.context_bin_records;
    }
    return lhs.fifo_rank > rhs.fifo_rank;
}

bool CostAwareTaskQueue::try_push(ExactBinTask task) {
    std::lock_guard<std::mutex> lock(mu_);
    if (closed_) {
        return false;
    }
    if (capacity_ > 0 && heap_.size() >= capacity_) {
        return false;
    }
    ReadyEntry entry;
    entry.cost = estimate_exact_bin_task_cost(task);
    entry.task = std::move(task);
    entry.fifo_rank = next_fifo_rank_++;
    heap_.push(std::move(entry));
    cv_not_empty_.notify_one();
    return true;
}

bool CostAwareTaskQueue::push(ExactBinTask task) {
    std::unique_lock<std::mutex> lock(mu_);
    cv_not_full_.wait(lock, [this]() {
        return closed_ || capacity_ == 0 || heap_.size() < capacity_;
    });
    if (closed_) {
        return false;
    }
    ReadyEntry entry;
    entry.cost = estimate_exact_bin_task_cost(task);
    entry.task = std::move(task);
    entry.fifo_rank = next_fifo_rank_++;
    heap_.push(std::move(entry));
    lock.unlock();
    cv_not_empty_.notify_one();
    return true;
}

bool CostAwareTaskQueue::pop(ExactBinTask& out) {
    std::unique_lock<std::mutex> lock(mu_);
    cv_not_empty_.wait(lock, [this]() {
        return closed_ || !heap_.empty();
    });
    if (heap_.empty()) {
        return false;
    }
    ReadyEntry entry = heap_.top();
    heap_.pop();
    out = std::move(entry.task);
    lock.unlock();
    cv_not_full_.notify_one();
    return true;
}

void CostAwareTaskQueue::close() {
    {
        std::lock_guard<std::mutex> lock(mu_);
        closed_ = true;
    }
    cv_not_empty_.notify_all();
    cv_not_full_.notify_all();
}

size_t CostAwareTaskQueue::capacity() const {
    return capacity_;
}

size_t CostAwareTaskQueue::size() const {
    std::lock_guard<std::mutex> lock(mu_);
    return heap_.size();
}

void DeterministicBinReducer::expect(int32_t tid, int32_t bin_index) {
    order_.push_back({tid, bin_index});
}

void DeterministicBinReducer::push_ready(
    int32_t tid,
    int32_t bin_index,
    ExactBinTaskResult result) {
    ready_[{tid, bin_index}] = std::move(result);
}

std::string render_parallel_progress(const ParallelProgressSnapshot& snapshot) {
    std::ostringstream oss;
    oss << "[Pipeline][parallel]"
        << " reads_scanned=" << snapshot.reads_scanned
        << " raw_bins_discovered=" << snapshot.raw_bins_discovered
        << " tasks_materialized=" << snapshot.tasks_materialized
        << " tasks_running=" << snapshot.tasks_running
        << " tasks_completed=" << snapshot.tasks_completed
        << " reducer_committed=" << snapshot.reducer_committed
        << " snapshot_window_size=" << snapshot.snapshot_window_size
        << " ready_queue_depth=" << snapshot.ready_queue_depth
        << " result_queue_depth=" << snapshot.result_queue_depth;
    return oss.str();
}

}  // namespace placer
