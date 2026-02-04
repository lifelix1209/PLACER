#include "task_queue.h"
#include <iostream>
#include <filesystem>
#include <sstream>
#include <cstdint>
#include <arpa/inet.h>  // for htonl, ntohl

namespace placer {

// ============= Endian-Safe Serialization Helpers =============

namespace endian {

// 检查是否为小端序
inline bool is_little_endian() {
    static const uint16_t test = 0x0001;
    return reinterpret_cast<const uint8_t*>(&test)[0] == 0x01;
}

// 字节序转换：逐字节反转
inline uint32_t swap32(uint32_t v) {
    return ((v & 0x000000FFU) << 24) |
           ((v & 0x0000FF00U) << 8) |
           ((v & 0x00FF0000U) >> 8) |
           ((v & 0xFF000000U) >> 24);
}

inline uint16_t swap16(uint16_t v) {
    return ((v & 0x00FFU) << 8) |
           ((v & 0xFF00U) >> 8);
}

inline int32_t swap32(int32_t v) {
    return static_cast<int32_t>(swap32(static_cast<uint32_t>(v)));
}

// 网络字节序转换（小端转大端）
inline uint32_t to_network_u32(uint32_t v) {
    return is_little_endian() ? swap32(v) : v;
}

inline int32_t to_network_i32(int32_t v) {
    return is_little_endian() ? swap32(v) : v;
}

inline uint16_t to_network_u16(uint16_t v) {
    return is_little_endian() ? swap16(v) : v;
}

// 从网络字节序转换（大端转本地）
inline uint32_t from_network_u32(uint32_t v) {
    return to_network_u32(v);  // 对称操作
}

inline int32_t from_network_i32(int32_t v) {
    return to_network_i32(v);  // 对称操作
}

inline uint16_t from_network_u16(uint16_t v) {
    return to_network_u16(v);  // 对称操作
}

}  // namespace endian

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

    // 写入版本（网络字节序）
    uint32_t version_net = endian::to_network_u32(TaskData::VERSION);
    out.write(reinterpret_cast<const char*>(&version_net), sizeof(uint32_t));

    // 写入任务类型（作为 uint32_t 存储）
    uint32_t type_net = endian::to_network_u32(static_cast<uint32_t>(task.type));
    out.write(reinterpret_cast<const char*>(&type_net), sizeof(uint32_t));

    // 写入 window_id 长度和内容
    uint32_t id_len_net = endian::to_network_u32(static_cast<uint32_t>(task.window_id.size()));
    out.write(reinterpret_cast<const char*>(&id_len_net), sizeof(uint32_t));
    out.write(task.window_id.data(), static_cast<std::streamsize>(task.window_id.size()));

    // 写入统计信息（网络字节序）
    uint32_t clip_bp_net = endian::to_network_u32(task.clip_bp);
    out.write(reinterpret_cast<const char*>(&clip_bp_net), sizeof(uint32_t));

    uint32_t sa_reads_net = endian::to_network_u32(task.sa_reads);
    out.write(reinterpret_cast<const char*>(&sa_reads_net), sizeof(uint32_t));

    uint32_t ins_events_net = endian::to_network_u32(task.ins_events);
    out.write(reinterpret_cast<const char*>(&ins_events_net), sizeof(uint32_t));

    // 写入 reads 数量（网络字节序）
    uint32_t read_count_net = endian::to_network_u32(static_cast<uint32_t>(task.reads.size()));
    out.write(reinterpret_cast<const char*>(&read_count_net), sizeof(uint32_t));

    // 写入每个 ReadSketch
    for (const auto& read : task.reads) {
        // qname
        uint32_t qname_len_net = endian::to_network_u32(static_cast<uint32_t>(read.qname.size()));
        out.write(reinterpret_cast<const char*>(&qname_len_net), sizeof(uint32_t));
        out.write(read.qname.data(), static_cast<std::streamsize>(read.qname.size()));

        // 基本字段（网络字节序）
        int32_t tid_net = endian::to_network_i32(read.tid);
        out.write(reinterpret_cast<const char*>(&tid_net), sizeof(int32_t));

        int32_t pos_net = endian::to_network_i32(read.pos);
        out.write(reinterpret_cast<const char*>(&pos_net), sizeof(int32_t));

        int32_t end_pos_net = endian::to_network_i32(read.end_pos);
        out.write(reinterpret_cast<const char*>(&end_pos_net), sizeof(int32_t));

        uint16_t flag_net = endian::to_network_u16(read.flag);
        out.write(reinterpret_cast<const char*>(&flag_net), sizeof(uint16_t));

        uint8_t mapq_val = read.mapq;
        out.write(reinterpret_cast<const char*>(&mapq_val), sizeof(uint8_t));

        // sequence
        uint32_t seq_len_net = endian::to_network_u32(static_cast<uint32_t>(read.sequence.size()));
        out.write(reinterpret_cast<const char*>(&seq_len_net), sizeof(uint32_t));
        out.write(read.sequence.data(), static_cast<std::streamsize>(read.sequence.size()));

        // 信号字段（网络字节序）
        int32_t total_clip_len_net = endian::to_network_i32(read.total_clip_len);
        out.write(reinterpret_cast<const char*>(&total_clip_len_net), sizeof(int32_t));

        uint8_t has_large = read.has_large_insertion ? 1 : 0;
        out.write(reinterpret_cast<const char*>(&has_large), sizeof(uint8_t));

        // CIGAR ops
        uint32_t cigar_count_net = endian::to_network_u32(static_cast<uint32_t>(read.cigar_ops.size()));
        out.write(reinterpret_cast<const char*>(&cigar_count_net), sizeof(uint32_t));
        for (const auto& [op, len] : read.cigar_ops) {
            out.write(&op, sizeof(char));
            int32_t len_net = endian::to_network_i32(len);
            out.write(reinterpret_cast<const char*>(&len_net), sizeof(int32_t));
        }

        // SA targets
        uint32_t sa_count_net = endian::to_network_u32(static_cast<uint32_t>(read.sa_targets.size()));
        out.write(reinterpret_cast<const char*>(&sa_count_net), sizeof(uint32_t));
        for (const auto& [sa_tid, sa_pos] : read.sa_targets) {
            int32_t sa_tid_net = endian::to_network_i32(sa_tid);
            out.write(reinterpret_cast<const char*>(&sa_tid_net), sizeof(int32_t));

            int32_t sa_pos_net = endian::to_network_i32(sa_pos);
            out.write(reinterpret_cast<const char*>(&sa_pos_net), sizeof(int32_t));
        }

        // MD tag
        uint8_t has_md = read.has_md ? 1 : 0;
        out.write(reinterpret_cast<const char*>(&has_md), sizeof(uint8_t));
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

bool TaskSerializer::load_task(const std::string& filename, TaskData& task) {
    std::ifstream in(filename, std::ios::binary);
    if (!in) return false;

    // 读取版本
    uint32_t version_net;
    in.read(reinterpret_cast<char*>(&version_net), sizeof(uint32_t));
    uint32_t version = endian::from_network_u32(version_net);
    if (version != TaskData::VERSION) {
        return false;  // 版本不匹配
    }

    // 读取任务类型
    uint32_t type_net;
    in.read(reinterpret_cast<char*>(&type_net), sizeof(uint32_t));
    task.type = static_cast<TaskType>(endian::from_network_u32(type_net));

    // 读取 window_id
    uint32_t id_len_net;
    in.read(reinterpret_cast<char*>(&id_len_net), sizeof(uint32_t));
    uint32_t id_len = endian::from_network_u32(id_len_net);
    task.window_id.resize(id_len);
    in.read(&task.window_id[0], static_cast<std::streamsize>(id_len));

    // 读取统计信息
    uint32_t clip_bp_net, sa_reads_net, ins_events_net;
    in.read(reinterpret_cast<char*>(&clip_bp_net), sizeof(uint32_t));
    in.read(reinterpret_cast<char*>(&sa_reads_net), sizeof(uint32_t));
    in.read(reinterpret_cast<char*>(&ins_events_net), sizeof(uint32_t));
    task.clip_bp = endian::from_network_u32(clip_bp_net);
    task.sa_reads = endian::from_network_u32(sa_reads_net);
    task.ins_events = endian::from_network_u32(ins_events_net);

    // 读取 reads 数量
    uint32_t read_count_net;
    in.read(reinterpret_cast<char*>(&read_count_net), sizeof(uint32_t));
    uint32_t read_count = endian::from_network_u32(read_count_net);

    // 读取每个 ReadSketch
    task.reads.clear();
    task.reads.reserve(read_count);

    for (uint32_t i = 0; i < read_count; ++i) {
        ReadSketch read;

        // qname
        uint32_t qname_len_net;
        in.read(reinterpret_cast<char*>(&qname_len_net), sizeof(uint32_t));
        uint32_t qname_len = endian::from_network_u32(qname_len_net);
        read.qname.resize(qname_len);
        in.read(&read.qname[0], static_cast<std::streamsize>(qname_len));

        // 基本字段
        int32_t tid_net, pos_net, end_pos_net;
        in.read(reinterpret_cast<char*>(&tid_net), sizeof(int32_t));
        in.read(reinterpret_cast<char*>(&pos_net), sizeof(int32_t));
        in.read(reinterpret_cast<char*>(&end_pos_net), sizeof(int32_t));
        read.tid = endian::from_network_i32(tid_net);
        read.pos = endian::from_network_i32(pos_net);
        read.end_pos = endian::from_network_i32(end_pos_net);

        uint16_t flag_net;
        in.read(reinterpret_cast<char*>(&flag_net), sizeof(uint16_t));
        read.flag = endian::from_network_u16(flag_net);

        uint8_t mapq_val;
        in.read(reinterpret_cast<char*>(&mapq_val), sizeof(uint8_t));
        read.mapq = mapq_val;

        // sequence
        uint32_t seq_len_net;
        in.read(reinterpret_cast<char*>(&seq_len_net), sizeof(uint32_t));
        uint32_t seq_len = endian::from_network_u32(seq_len_net);
        read.sequence.resize(seq_len);
        in.read(&read.sequence[0], static_cast<std::streamsize>(seq_len));

        // 信号字段
        int32_t total_clip_len_net;
        in.read(reinterpret_cast<char*>(&total_clip_len_net), sizeof(int32_t));
        read.total_clip_len = endian::from_network_i32(total_clip_len_net);

        uint8_t has_large;
        in.read(reinterpret_cast<char*>(&has_large), sizeof(uint8_t));
        read.has_large_insertion = (has_large != 0);

        // CIGAR ops
        uint32_t cigar_count_net;
        in.read(reinterpret_cast<char*>(&cigar_count_net), sizeof(uint32_t));
        uint32_t cigar_count = endian::from_network_u32(cigar_count_net);
        read.cigar_ops.clear();
        for (uint32_t j = 0; j < cigar_count; ++j) {
            char op;
            int32_t len_net;
            in.read(&op, sizeof(char));
            in.read(reinterpret_cast<char*>(&len_net), sizeof(int32_t));
            read.cigar_ops.push_back({op, endian::from_network_i32(len_net)});
        }

        // SA targets
        uint32_t sa_count_net;
        in.read(reinterpret_cast<char*>(&sa_count_net), sizeof(uint32_t));
        uint32_t sa_count = endian::from_network_u32(sa_count_net);
        read.sa_targets.clear();
        for (uint32_t j = 0; j < sa_count; ++j) {
            int32_t sa_tid_net, sa_pos_net;
            in.read(reinterpret_cast<char*>(&sa_tid_net), sizeof(int32_t));
            in.read(reinterpret_cast<char*>(&sa_pos_net), sizeof(int32_t));
            read.sa_targets.push_back({
                endian::from_network_i32(sa_tid_net),
                endian::from_network_i32(sa_pos_net)
            });
        }

        // MD tag
        uint8_t has_md;
        in.read(reinterpret_cast<char*>(&has_md), sizeof(uint8_t));
        read.has_md = (has_md != 0);

        task.reads.push_back(std::move(read));
    }

    return true;
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

SerializedTask::SerializedTask(TaskType type, std::string window_id, std::vector<ReadSketch> reads)
    : Task(type, std::move(window_id)), reads_(std::move(reads)) {}

void SerializedTask::execute() {
    // TODO: 实现真正的任务执行逻辑
    // 此时可通过 get_reads() 访问读数据
}

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
                                  std::vector<ReadSketch> reads) {
    std::unique_lock<std::mutex> lock(mutex_);
    if (closed_) return false;

    // 通过值传递，reads 被移动进任务，避免拷贝
    auto task = std::make_unique<SerializedTask>(type, window_id, std::move(reads));
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
