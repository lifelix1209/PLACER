#ifndef PLACER_TASK_QUEUE_H
#define PLACER_TASK_QUEUE_H

#include "bam_reader.h"
#include <queue>
#include <string>
#include <vector>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <fstream>

namespace placer {

enum class TaskType {
    COMPONENT_BUILD,
    LOCAL_ALIGN,
    ASSEMBLY,
    COLLAPSE,
    GENOTYPE
};

/**
 * TaskData: 序列化任务数据的通用结构
 */
struct TaskData {
    TaskType type;
    std::string window_id;

    // 可序列化的读片段数据
    std::vector<ReadSketch> reads;

    // 窗口统计信息
    uint32_t clip_bp = 0;
    uint32_t sa_reads = 0;
    uint32_t ins_events = 0;

    // 序列化版本
    static constexpr uint32_t VERSION = 1;
};

class TaskSerializer {
public:
    explicit TaskSerializer(const std::string& output_dir);

    // 保存任务到磁盘
    bool save_task(const TaskData& task);

    // 批量保存任务
    bool save_tasks(const std::vector<TaskData>& tasks);

    // 从磁盘加载任务
    bool load_next_task(TaskData& task);

    // 关闭序列化器
    void close();

    // 获取待处理任务数量
    size_t pending_count() const;

private:
    std::string output_dir_;
    std::ofstream task_file_;
    size_t pending_count_ = 0;
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

    // 提交序列化的任务数据
    bool submit_serialized(TaskType type, const std::string& window_id,
                          const std::vector<ReadSketch>& reads);

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
