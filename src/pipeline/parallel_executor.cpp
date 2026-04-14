#include "parallel_executor_internal.h"

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
        if ((bin_index - snapshot.bin_index) > cross_bin_context_bins) {
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

void CostAwareTaskQueue::push(ExactBinTask task) {
    std::unique_lock<std::mutex> lock(mu_);
    cv_not_full_.wait(lock, [this]() {
        return closed_ || capacity_ == 0 || heap_.size() < capacity_;
    });
    if (closed_) {
        return;
    }
    ReadyEntry entry;
    entry.cost = estimate_exact_bin_task_cost(task);
    entry.task = std::move(task);
    entry.fifo_rank = next_fifo_rank_++;
    heap_.push(std::move(entry));
    lock.unlock();
    cv_not_empty_.notify_one();
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
