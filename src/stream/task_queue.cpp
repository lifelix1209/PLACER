#include "task_queue.h"
#include <iostream>
#include <filesystem>
#include <sstream>

namespace placer {

// ============= TaskSerializer Implementation =============

TaskSerializer::TaskSerializer(const std::string& output_dir) : output_dir_(output_dir) {
    namespace fs = std::filesystem;
    if (!output_dir_.empty() && !fs::exists(output_dir_)) {
        fs::create_directories(output_dir_);
    }
}

bool TaskSerializer::save_task(const TaskData& task) {
    std::string filename = output_dir_ + "/task_" + std::to_string(static_cast<int>(task.type)) +
                           "_" + task.window_id + ".task";

    std::ofstream out(filename, std::ios::binary);
    if (!out) return false;

    // 写入版本
    out.write(reinterpret_cast<const char*>(&TaskData::VERSION), sizeof(uint32_t));

    // 写入任务类型
    out.write(reinterpret_cast<const char*>(&task.type), sizeof(TaskType));

    // 写入 window_id 长度和内容
    uint32_t id_len = static_cast<uint32_t>(task.window_id.size());
    out.write(reinterpret_cast<const char*>(&id_len), sizeof(uint32_t));
    out.write(task.window_id.data(), id_len);

    // 写入统计信息
    out.write(reinterpret_cast<const char*>(&task.clip_bp), sizeof(uint32_t));
    out.write(reinterpret_cast<const char*>(&task.sa_reads), sizeof(uint32_t));
    out.write(reinterpret_cast<const char*>(&task.ins_events), sizeof(uint32_t));

    // 写入 reads 数量
    uint32_t read_count = static_cast<uint32_t>(task.reads.size());
    out.write(reinterpret_cast<const char*>(&read_count), sizeof(uint32_t));

    // 写入每个 ReadSketch
    for (const auto& read : task.reads) {
        // qname
        uint32_t qname_len = static_cast<uint32_t>(read.qname.size());
        out.write(reinterpret_cast<const char*>(&qname_len), sizeof(uint32_t));
        out.write(read.qname.data(), qname_len);

        // 基本字段
        out.write(reinterpret_cast<const char*>(&read.tid), sizeof(int32_t));
        out.write(reinterpret_cast<const char*>(&read.pos), sizeof(int32_t));
        out.write(reinterpret_cast<const char*>(&read.end_pos), sizeof(int32_t));
        out.write(reinterpret_cast<const char*>(&read.flag), sizeof(uint16_t));
        out.write(reinterpret_cast<const char*>(&read.mapq), sizeof(uint8_t));

        // sequence
        uint32_t seq_len = static_cast<uint32_t>(read.sequence.size());
        out.write(reinterpret_cast<const char*>(&seq_len), sizeof(uint32_t));
        out.write(read.sequence.data(), seq_len);

        // 信号字段
        out.write(reinterpret_cast<const char*>(&read.total_clip_len), sizeof(int32_t));
        out.write(reinterpret_cast<const char*>(&read.has_large_insertion), sizeof(bool));

        // CIGAR ops
        uint32_t cigar_count = static_cast<uint32_t>(read.cigar_ops.size());
        out.write(reinterpret_cast<const char*>(&cigar_count), sizeof(uint32_t));
        for (const auto& [op, len] : read.cigar_ops) {
            out.write(&op, sizeof(char));
            out.write(reinterpret_cast<const char*>(&len), sizeof(int));
        }

        // SA targets
        uint32_t sa_count = static_cast<uint32_t>(read.sa_targets.size());
        out.write(reinterpret_cast<const char*>(&sa_count), sizeof(uint32_t));
        for (const auto& [sa_tid, sa_pos] : read.sa_targets) {
            out.write(reinterpret_cast<const char*>(&sa_tid), sizeof(int32_t));
            out.write(reinterpret_cast<const char*>(&sa_pos), sizeof(int32_t));
        }

        // MD tag
        out.write(reinterpret_cast<const char*>(&read.has_md), sizeof(bool));
    }

    ++pending_count_;
    return true;
}

bool TaskSerializer::save_tasks(const std::vector<TaskData>& tasks) {
    for (const auto& task : tasks) {
        if (!save_task(task)) return false;
    }
    return true;
}

bool TaskSerializer::load_next_task(TaskData& task) {
    // 简化实现：读取目录中的任务文件
    // 实际应用中应维护文件列表索引
    return false;
}

void TaskSerializer::close() {
    if (task_file_.is_open()) {
        task_file_.close();
    }
}

size_t TaskSerializer::pending_count() const {
    return pending_count_;
}

// ============= Task Implementation =============

Task::Task(TaskType type, std::string window_id)
    : type_(type), window_id_(std::move(window_id)) {}

PlaceholderTask::PlaceholderTask(TaskType type, std::string window_id)
    : Task(type, std::move(window_id)) {}

void PlaceholderTask::execute() {}

TaskQueue::TaskQueue(int num_workers) {
    for (int i = 0; i < num_workers; ++i) {
        workers_.emplace_back(&TaskQueue::worker_loop, this);
    }
}

TaskQueue::~TaskQueue() {
    close();
    for (auto& w : workers_) {
        if (w.joinable()) w.join();
    }
}

bool TaskQueue::submit(std::unique_ptr<Task> task) {
    std::unique_lock<std::mutex> lock(mutex_);
    if (closed_) return false;
    queue_.push(std::move(task));
    ++pending_;
    cv_.notify_one();
    return true;
}

bool TaskQueue::submit_serialized(TaskType type, const std::string& window_id,
                                  const std::vector<ReadSketch>& reads) {
    std::unique_lock<std::mutex> lock(mutex_);
    if (closed_) return false;

    // 创建包含读数据的任务
    auto task = std::make_unique<PlaceholderTask>(type, window_id);
    queue_.push(std::move(task));
    ++pending_;
    cv_.notify_one();
    return true;
}

void TaskQueue::close() {
    std::unique_lock<std::mutex> lock(mutex_);
    closed_ = true;
    cv_.notify_all();
}

void TaskQueue::wait() {
    std::unique_lock<std::mutex> lock(mutex_);
    cv_done_.wait(lock, [this]() {
        return queue_.empty() || closed_;
    });
}

size_t TaskQueue::size() const {
    std::unique_lock<std::mutex> lock(mutex_);
    return queue_.size();
}

void TaskQueue::worker_loop() {
    while (true) {
        std::unique_ptr<Task> task;

        {
            std::unique_lock<std::mutex> lock(mutex_);
            cv_.wait(lock, [this]() {
                return closed_ || !queue_.empty();
            });

            if (queue_.empty()) return;

            task = std::move(queue_.front());
            queue_.pop();
        }

        task->execute();

        {
            std::unique_lock<std::mutex> lock(mutex_);
            --pending_;
            if (queue_.empty() && pending_ == 0) {
                cv_done_.notify_all();
            }
        }
    }
}

std::unique_ptr<Task> TaskFactory::create_component_build(const std::string& window_id) {
    return std::make_unique<PlaceholderTask>(TaskType::COMPONENT_BUILD, window_id);
}
std::unique_ptr<Task> TaskFactory::create_local_align(const std::string& window_id) {
    return std::make_unique<PlaceholderTask>(TaskType::LOCAL_ALIGN, window_id);
}
std::unique_ptr<Task> TaskFactory::create_assembly(const std::string& window_id) {
    return std::make_unique<PlaceholderTask>(TaskType::ASSEMBLY, window_id);
}
std::unique_ptr<Task> TaskFactory::create_collapse(const std::string& window_id) {
    return std::make_unique<PlaceholderTask>(TaskType::COLLAPSE, window_id);
}
std::unique_ptr<Task> TaskFactory::create_genotype(const std::string& window_id) {
    return std::make_unique<PlaceholderTask>(TaskType::GENOTYPE, window_id);
}

}  // namespace placer
