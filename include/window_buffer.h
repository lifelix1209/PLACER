#ifndef PLACER_WINDOW_BUFFER_H
#define PLACER_WINDOW_BUFFER_H

#include "bam_reader.h"
#include <deque>
#include <unordered_map>
#include <vector>
#include <memory>
#include <cstdint>

namespace placer {

/**
 * WindowID: 避免字符串拼接的性能开销
 * 格式: (int64_t)chrom_tid << 32 | window_start_idx
 */
using WindowID = uint64_t;

/**
 * 窗口统计：归一化后的值，用于触发判断
 */
struct WindowStats {
    uint32_t total_reads = 0;
    uint32_t clip_bp = 0;
    uint32_t sa_reads = 0;
    uint32_t ins_events = 0;
    uint32_t bases_covered = 0;  // 归一化分母

    // Q99 thresholds (populated from WindowStatsCollector per chromosome)
    double q99_clip = 0.0;
    double q99_sa = 0.0;
    double q99_ins = 0.0;
};

/**
 * 窗口：包含分层 buffer
 */
struct Window {
    int32_t chrom_tid;
    int32_t start;
    int32_t end;

    WindowStats stats;
    bool triggered = false;

    // 分层 Buffer：存储副本以避免悬空指针
    std::vector<ReadSketch> priority_reads;  // SA/Split/Large-Indel
    std::vector<ReadSketch> normal_reads;   // 背景 reads
};

/**
 * WindowBuffer: 滑动窗口 + 分层采样 + 内存回收
 */
class WindowBuffer {
public:
    struct Config {
        int window_size = 10000;
        int window_step = 5000;
        int max_priority_reads = 50;
        int max_normal_reads = 200;
        int min_clip_bp = 50;
        int min_sa_reads = 3;
    };

    explicit WindowBuffer(Config config);

    // 添加一条 read，它会贡献给它覆盖的所有窗口
    void add_read(const ReadSketch& read);

    // 批量添加
    void add_reads(const std::vector<ReadSketch>& reads);

    // 真正的内存回收：返回被触发且完成封存的窗口，清理落后于 safe_pos 的所有窗口
    // 注意：此方法不再处理染色体切换，染色体切换由 add_read 自动处理
    std::vector<std::unique_ptr<Window>> seal_and_flush(int32_t safe_pos);

    // 染色体切换时的清理：刷新前一个染色体的所有窗口
    std::vector<std::unique_ptr<Window>> flush_current_chromosome();

    // 获取窗口统计
    const Window* get_window(WindowID id) const;
    const WindowStats& get_stats(WindowID id) const;
    std::vector<WindowID> list_windows() const;

    Config get_config() const { return config_; }

private:
    static WindowID make_id(int32_t chrom_tid, int32_t start) {
        return (static_cast<uint64_t>(static_cast<uint32_t>(chrom_tid)) << 32) |
               static_cast<uint32_t>(start);
    }

    static int32_t get_start(WindowID id) {
        return static_cast<int32_t>(id & 0xFFFFFFFF);
    }

    static int32_t get_chrom_tid(WindowID id) {
        return static_cast<int32_t>(id >> 32);
    }

    Window* get_or_create(int32_t chrom_tid, int32_t window_start);

    bool is_evidence_read(const ReadSketch& read) const;

    void add_to_stratified_buffer(Window* win, const ReadSketch& read);

    bool check_trigger(const WindowStats& stats, uint32_t bases_covered) const;

    Config config_;

    // 当前活跃染色体（用于检测染色体切换）
    int32_t current_chrom_tid_ = -1;

    // 活跃窗口
    std::unordered_map<WindowID, std::unique_ptr<Window>> active_windows_;

    // 维护顺序用于清理：(chrom_tid, window_start)
    std::deque<WindowID> window_order_;
};

}  // namespace placer

#endif  // PLACER_WINDOW_BUFFER_H
