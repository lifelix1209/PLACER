#include "task_queue.h"
#include <iostream>

namespace placer {

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
