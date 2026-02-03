#ifndef PLACER_TASK_QUEUE_H
#define PLACER_TASK_QUEUE_H

#include <queue>
#include <string>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>

namespace placer {

enum class TaskType {
    COMPONENT_BUILD,
    LOCAL_ALIGN,
    ASSEMBLY,
    COLLAPSE,
    GENOTYPE
};

class Task {
public:
    explicit Task(TaskType type, std::string window_id);
    virtual ~Task() = default;

    TaskType get_type() const { return type_; }
    const std::string& get_window_id() const { return window_id_; }
    virtual void execute() = 0;

private:
    TaskType type_;
    std::string window_id_;
};

class PlaceholderTask : public Task {
public:
    PlaceholderTask(TaskType type, std::string window_id);
    void execute() override;
};

/**
 * TaskQueue: 简单的同步任务队列
 * - 单生产者（Stream Layer）提交任务
 * - 多消费者（Task Layer workers）执行任务
 */
class TaskQueue {
public:
    explicit TaskQueue(int num_workers = 4);

    ~TaskQueue();

    // 不可复制
    TaskQueue(const TaskQueue&) = delete;
    TaskQueue& operator=(const TaskQueue&) = delete;

    // 提交任务，成功返回 true，队列已关闭返回 false
    bool submit(std::unique_ptr<Task> task);

    // 关闭队列，停止接受新任务
    void close();

    // 等待所有待处理任务完成
    void wait();

    size_t size() const;

private:
    void worker_loop();

    std::queue<std::unique_ptr<Task>> queue_;
    std::vector<std::thread> workers_;
    mutable std::mutex mutex_;
    std::condition_variable cv_;
    std::condition_variable cv_done_;
    bool closed_ = false;
    size_t pending_ = 0;  // 待处理任务数
};

class TaskFactory {
public:
    static std::unique_ptr<Task> create_component_build(const std::string& window_id);
    static std::unique_ptr<Task> create_local_align(const std::string& window_id);
    static std::unique_ptr<Task> create_assembly(const std::string& window_id);
    static std::unique_ptr<Task> create_collapse(const std::string& window_id);
    static std::unique_ptr<Task> create_genotype(const std::string& window_id);
};

}  // namespace placer

#endif  // PLACER_TASK_QUEUE_H
