#include "window_buffer.h"
#include <random>
#include <algorithm>

namespace placer {

WindowBuffer::WindowBuffer(Config config) : config_(std::move(config)) {}

void WindowBuffer::add_read(const ReadSketch& read) {
    if (read.pos >= read.end_pos) return;

    // 检测染色体切换：当 read 的 tid 与当前活跃染色体不同时
    // 需要先刷新前一个染色体的所有窗口
    if (current_chrom_tid_ != read.tid) {
        if (current_chrom_tid_ != -1) {
            // 切换前，先刷新前一个染色体的所有窗口（无论位置）
            flush_current_chromosome();
        }
        current_chrom_tid_ = read.tid;
    }

    // 空间逻辑修正：read 贡献给它覆盖的所有窗口
    int32_t start_idx = read.pos / config_.window_size;
    int32_t end_idx = (read.end_pos - 1) / config_.window_size;  // -1 确保在区间内

    for (int32_t idx = start_idx; idx <= end_idx; ++idx) {
        int32_t window_start = idx * config_.window_size;
        Window* win = get_or_create(read.tid, window_start);

        // 累加该 read 在此窗口内的覆盖度
        int32_t win_end = std::min(read.end_pos, win->end);
        int32_t overlap = win_end - std::max(read.pos, win->start);
        if (overlap > 0) {
            win->stats.bases_covered += overlap;
        }

        // 累加统计信号（简化：每条 read 的信号只记一次，不管跨多少窗口）
        // 更精细的做法是只统计落在窗口内的 event
        if (read.has_sa) {
            win->stats.sa_reads++;
        }

        // 计算 clip_bp 和 ins_events（简化处理）
        uint32_t clip_bp = 0;
        uint32_t ins_events = 0;
        for (const auto& [op, len] : read.cigar_ops) {
            if (op == 'S') clip_bp += len;
            if (op == 'I') ins_events += len;
        }
        if (clip_bp > 0) {
            win->stats.clip_bp += clip_bp;
        }
        if (ins_events > 0) {
            win->stats.ins_events += ins_events;
        }

        // 分层采样
        add_to_stratified_buffer(win, read);

        // 触发检查
        if (!win->triggered && check_trigger(win->stats, win->stats.bases_covered)) {
            win->triggered = true;
        }
    }
}

void WindowBuffer::add_reads(const std::vector<ReadSketch>& reads) {
    for (const auto& read : reads) {
        add_read(read);
    }
}

std::vector<std::unique_ptr<Window>> WindowBuffer::seal_and_flush(int32_t safe_pos) {
    std::vector<std::unique_ptr<Window>> sealed;

    // 真正的内存回收：清理所有 end < safe_pos 的窗口
    // 注意：染色体切换已在 add_read 中处理，这里只处理位置推进
    while (!window_order_.empty()) {
        WindowID old_id = window_order_.front();
        int32_t old_start = get_start(old_id);

        // 检查是否超出安全线
        if (old_start >= safe_pos) {
            break;  // 还在活跃区，停止
        }

        auto it = active_windows_.find(old_id);
        if (it != active_windows_.end()) {
            if (it->second->triggered) {
                // 移动 ownership 到输出
                sealed.push_back(std::move(it->second));
            }
            // 无论是否触发，都从内存中移除
            active_windows_.erase(it);
        }
        window_order_.pop_front();
    }

    return sealed;
}

std::vector<std::unique_ptr<Window>> WindowBuffer::flush_current_chromosome() {
    std::vector<std::unique_ptr<Window>> flushed;

    // 清理当前染色体的所有窗口
    // 这在染色体切换时被调用，确保旧染色体不会占用内存
    while (!window_order_.empty()) {
        WindowID old_id = window_order_.front();
        int32_t old_tid = get_chrom_tid(old_id);

        // 如果是当前染色体，继续清理
        if (old_tid != current_chrom_tid_) {
            break;  // 遇到其他染色体，停止
        }

        auto it = active_windows_.find(old_id);
        if (it != active_windows_.end()) {
            // 只保留触发的窗口
            if (it->second->triggered) {
                flushed.push_back(std::move(it->second));
            }
            active_windows_.erase(it);
        }
        window_order_.pop_front();
    }

    return flushed;
}

const Window* WindowBuffer::get_window(WindowID id) const {
    auto it = active_windows_.find(id);
    if (it != active_windows_.end()) {
        return it->second.get();
    }
    return nullptr;
}

const WindowStats& WindowBuffer::get_stats(WindowID id) const {
    static const WindowStats empty_stats{};
    auto it = active_windows_.find(id);
    if (it != active_windows_.end()) {
        return it->second->stats;
    }
    return empty_stats;
}

std::vector<WindowID> WindowBuffer::list_windows() const {
    std::vector<WindowID> ids;
    for (const auto& [id, _] : active_windows_) {
        ids.push_back(id);
    }
    return ids;
}

Window* WindowBuffer::get_or_create(int32_t chrom_tid, int32_t window_start) {
    WindowID id = make_id(chrom_tid, window_start);
    auto it = active_windows_.find(id);
    if (it != active_windows_.end()) {
        return it->second.get();
    }

    auto win = std::make_unique<Window>();
    win->chrom_tid = chrom_tid;
    win->start = window_start;
    win->end = window_start + config_.window_size;

    Window* ptr = win.get();
    active_windows_[id] = std::move(win);
    window_order_.push_back(id);

    return ptr;
}

bool WindowBuffer::is_evidence_read(const ReadSketch& read) const {
    // 定义什么是"证据 read"
    return read.has_sa ||
           read.mapq < 20 ||
           [&read]() {
               for (const auto& [op, len] : read.cigar_ops) {
                   if ((op == 'S' || op == 'I' || op == 'D') && len >= 10) return true;
               }
               return false;
           }();
}

void WindowBuffer::add_to_stratified_buffer(Window* win, const ReadSketch& read) {
    static std::random_device rd;
    static std::mt19937 gen(rd());

    bool is_priority = is_evidence_read(read);

    if (is_priority) {
        // 证据 reads：优先保留
        if (win->priority_reads.size() < static_cast<size_t>(config_.max_priority_reads)) {
            win->priority_reads.push_back(read);
        } else {
            // 随机替换，保持多样性
            std::uniform_int_distribution<> dis(0, win->priority_reads.size() - 1);
            int idx = dis(gen);
            win->priority_reads[idx] = read;
        }
    } else {
        // 背景 reads：蓄水池采样
        if (win->normal_reads.size() < static_cast<size_t>(config_.max_normal_reads)) {
            win->normal_reads.push_back(read);
        } else {
            // 蓄水池：替换概率 k/n
            std::uniform_int_distribution<> dis(0, win->normal_reads.size());
            int idx = dis(gen);
            if (idx < static_cast<int>(win->normal_reads.size())) {
                win->normal_reads[idx] = read;
            }
        }
    }
}

bool WindowBuffer::check_trigger(const WindowStats& stats, uint32_t bases_covered) const {
    // 绝对阈值
    if (stats.sa_reads >= static_cast<uint32_t>(config_.min_sa_reads)) return true;
    if (stats.clip_bp >= static_cast<uint32_t>(config_.min_clip_bp)) return true;

    // 相对阈值（归一化）：防止高覆盖度区域的假阳性
    if (bases_covered > 0) {
        double clip_ratio = static_cast<double>(stats.clip_bp) / bases_covered;
        if (clip_ratio > 0.1) return true;  // >10% 的 clip 比例异常
    }

    return false;
}

}  // namespace placer
