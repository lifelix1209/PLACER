#ifndef PLACER_PARALLEL_EXECUTOR_INTERNAL_H
#define PLACER_PARALLEL_EXECUTOR_INTERNAL_H

#include "pipeline.h"

#include <cstddef>
#include <condition_variable>
#include <cstdint>
#include <map>
#include <memory>
#include <mutex>
#include <queue>
#include <string>
#include <utility>
#include <vector>

namespace placer {

using SharedRecordBatch = std::shared_ptr<const std::vector<BufferedRecord>>;

struct ExactBinSnapshot {
    int32_t tid = -1;
    int32_t bin_index = -1;
    SharedRecordBatch records;
};

struct ExactTaskRecordRef {
    const bam1_t* record = nullptr;
    bool from_current_bin = false;
};

struct ExactBinTask {
    int32_t tid = -1;
    int32_t bin_index = -1;
    int32_t bin_start = -1;
    int32_t bin_end = -1;
    size_t current_bin_records = 0;
    size_t context_bin_records = 0;
    size_t merged_input_records = 0;
    std::vector<ExactTaskRecordRef> records;
    std::vector<SharedRecordBatch> owners;
};

ExactBinTask build_exact_bin_task(
    int32_t tid,
    int32_t bin_index,
    const std::vector<ExactBinSnapshot>& snapshots,
    int32_t bin_size,
    int32_t cross_bin_context_bins,
    int32_t window_bin_slack_bp);

template <typename T>
class BoundedTaskQueue {
public:
    explicit BoundedTaskQueue(size_t capacity) : capacity_(capacity) {}

    bool try_push(T value) {
        std::lock_guard<std::mutex> lock(mu_);
        if (closed_) {
            return false;
        }
        if (capacity_ > 0 && queue_.size() >= capacity_) {
            return false;
        }
        queue_.push(std::move(value));
        cv_not_empty_.notify_one();
        return true;
    }

    void push(T value) {
        std::unique_lock<std::mutex> lock(mu_);
        cv_not_full_.wait(lock, [this]() {
            return closed_ || capacity_ == 0 || queue_.size() < capacity_;
        });
        if (closed_) {
            return;
        }
        queue_.push(std::move(value));
        lock.unlock();
        cv_not_empty_.notify_one();
    }

    bool try_pop(T& out) {
        std::lock_guard<std::mutex> lock(mu_);
        if (queue_.empty()) {
            return false;
        }
        out = std::move(queue_.front());
        queue_.pop();
        cv_not_full_.notify_one();
        return true;
    }

    bool pop(T& out) {
        std::unique_lock<std::mutex> lock(mu_);
        cv_not_empty_.wait(lock, [this]() {
            return closed_ || !queue_.empty();
        });
        if (queue_.empty()) {
            return false;
        }
        out = std::move(queue_.front());
        queue_.pop();
        lock.unlock();
        cv_not_full_.notify_one();
        return true;
    }

    void close() {
        {
            std::lock_guard<std::mutex> lock(mu_);
            closed_ = true;
        }
        cv_not_empty_.notify_all();
        cv_not_full_.notify_all();
    }

    size_t capacity() const {
        return capacity_;
    }

    size_t size() const {
        std::lock_guard<std::mutex> lock(mu_);
        return queue_.size();
    }

private:
    std::queue<T> queue_;
    mutable std::mutex mu_;
    std::condition_variable cv_not_empty_;
    std::condition_variable cv_not_full_;
    size_t capacity_ = 0;
    bool closed_ = false;
};

struct ExactBinTaskCost {
    int32_t tid = -1;
    int32_t bin_index = -1;
    size_t merged_input_records = 0;
    size_t current_bin_records = 0;
    size_t context_bin_records = 0;
};

ExactBinTaskCost estimate_exact_bin_task_cost(const ExactBinTask& task);

class CostAwareTaskQueue {
public:
    explicit CostAwareTaskQueue(size_t capacity);

    bool try_push(ExactBinTask task);
    void push(ExactBinTask task);
    bool pop(ExactBinTask& out);
    void close();
    size_t capacity() const;
    size_t size() const;

private:
    struct ReadyEntry {
        ExactBinTask task;
        ExactBinTaskCost cost;
        uint64_t fifo_rank = 0;
    };

    struct ReadyEntryLowerPriority {
        bool operator()(const ReadyEntry& lhs, const ReadyEntry& rhs) const;
    };

    std::priority_queue<
        ReadyEntry,
        std::vector<ReadyEntry>,
        ReadyEntryLowerPriority> heap_;
    mutable std::mutex mu_;
    std::condition_variable cv_not_empty_;
    std::condition_variable cv_not_full_;
    size_t capacity_ = 0;
    uint64_t next_fifo_rank_ = 0;
    bool closed_ = false;
};

struct ExactBinTaskResult {
    int32_t tid = -1;
    int32_t bin_index = -1;
    PipelineResult partial;
    double elapsed_s = 0.0;
};

class DeterministicBinReducer {
public:
    void expect(int32_t tid, int32_t bin_index);
    void push_ready(int32_t tid, int32_t bin_index, ExactBinTaskResult result);

    template <typename CommitFn>
    void drain(CommitFn&& commit) {
        while (next_index_ < order_.size()) {
            const auto& key = order_[next_index_];
            auto it = ready_.find(key);
            if (it == ready_.end()) {
                break;
            }
            commit(it->second.tid, it->second.bin_index, it->second);
            ready_.erase(it);
            ++next_index_;
        }
    }

private:
    std::vector<std::pair<int32_t, int32_t>> order_;
    std::map<std::pair<int32_t, int32_t>, ExactBinTaskResult> ready_;
    size_t next_index_ = 0;
};

struct ParallelProgressSnapshot {
    int64_t reads_scanned = 0;
    int64_t raw_bins_discovered = 0;
    int64_t tasks_materialized = 0;
    int64_t tasks_running = 0;
    int64_t tasks_completed = 0;
    int64_t reducer_committed = 0;
    size_t snapshot_window_size = 0;
    size_t ready_queue_depth = 0;
    size_t result_queue_depth = 0;
};

std::string render_parallel_progress(const ParallelProgressSnapshot& snapshot);

}  // namespace placer

#endif
