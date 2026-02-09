#include "component_builder.h"
#include <algorithm>
#include <cmath>
#include <execution>

namespace placer {

ComponentBuilder::ComponentBuilder(ComponentBuilderConfig config)
    : config_(std::move(config)) {}

std::vector<Component> ComponentBuilder::build(const TaskData& task_data) {
    return build(task_data.reads, 0);  // TODO: chrom_tid from task_data if available
}

std::vector<Component> ComponentBuilder::build(const std::vector<ReadSketch>& reads, int32_t chrom_tid) {
    total_anchors_ = 0;
    total_components_ = 0;
    breakers_triggered_ = 0;

    // ===== Step 1: 多重锚点提取 =====
    std::vector<Anchor> all_anchors;
    all_anchors.reserve(reads.size() * 4);  // 预分配：考虑HardClip后平均每条read约4个anchor

    for (size_t i = 0; i < reads.size(); ++i) {
        // MapQ 过滤：跳过低质量比对
        if (reads[i].mapq < config_.min_mapq) {
            continue;
        }

        auto anchors = extract_anchors(reads[i], chrom_tid);

        // 修正 read_idx
        for (auto& anchor : anchors) {
            anchor.read_idx = i;
            anchor.mapq = reads[i].mapq;  // 记录 MapQ
            // 从 BAM flag 提取链信息：0x10 (16) = reverse strand
            anchor.strand = !(reads[i].flag & 0x10);  // true=forward, false=reverse
        }

        if (!anchors.empty()) {
            merge_close_anchors(anchors);  // 同一 read 内去重
            all_anchors.insert(all_anchors.end(), anchors.begin(), anchors.end());
        }
    }

    total_anchors_ = all_anchors.size();
    if (total_anchors_ == 0) {
        return {};  // 无锚点，无组件
    }

    // 全局排序：先按 tid，再按 pos（使用并行排序优化）
    std::sort(std::execution::par, all_anchors.begin(), all_anchors.end());

    // ===== Step 2: 密度聚类 =====
    AnchorSpan anchor_span(all_anchors.data(), all_anchors.size());
    auto cluster_ranges = density_clustering(anchor_span);

    // ===== Step 3: 递归式 Cluster Breaker & 构建 Components =====
    std::vector<Component> components;
    components.reserve(cluster_ranges.size() * 2);

    int32_t component_id = 0;
    for (const auto& range : cluster_ranges) {
        if (range.empty()) continue;
        recursive_cluster_breaker(anchor_span, range, components, chrom_tid, component_id, 0);
    }

    total_components_ = components.size();
    return components;
}

std::vector<Anchor> ComponentBuilder::extract_anchors(const ReadSketch& read, int32_t chrom_tid) {
    std::vector<Anchor> anchors;

    // ===== 0. 完整 CIGAR 解析 =====
    // 计算参考消耗：M/D/N/EQ/X 消耗参考坐标
    int32_t current_pos = read.pos;
    int32_t read_consumed = 0;

    // 记录第一个和最后一个 S/H 操作
    int32_t first_clip_len = 0;
    int32_t last_clip_len = 0;
    char first_clip_op = 0;
    char last_clip_op = 0;

    size_t op_idx = 0;
    for (const auto& [op, len] : read.cigar_ops) {
        if (op == 'S' || op == 'H') {
            if (op_idx == 0) {
                first_clip_len = len;
                first_clip_op = op;
            }
            last_clip_len = len;
            last_clip_op = op;
        }

        // 累加参考消耗以追踪位置
        if (op == 'M' || op == 'D' || op == 'N' || op == '=' || op == 'X') {
            current_pos += len;
        }
        // I/P 不消耗参考坐标

        op_idx++;
    }

    // ===== 1. Clip 信号处理 (同时支持 SoftClip 和 HardClip) =====
    // 左端 Clip -> Right Break (5' 端断点)
    if (first_clip_len >= config_.min_clip_len) {
        anchors.push_back({
            chrom_tid,
            read.pos,           // Clip 开始位置 = 断点位置
            1,                  // RightBreak
            1,                  // CLIP
            0,                  // read_idx - 临时占位
            0,                  // mapq - 临时占位
            false               // strand - 临时占位
        });
    }

    // 右端 Clip -> Left Break (3' 端断点)
    if (last_clip_len >= config_.min_clip_len) {
        // 计算右端位置：CIGAR 最后一个操作结束后的位置
        // 右端断点位置 = end_pos - 1 (0-based)，或者 pos + clip_len
        int32_t right_break_pos = read.end_pos - 1;
        anchors.push_back({
            chrom_tid,
            right_break_pos,    // Clip 结束位置 = 断点位置
            0,                  // LeftBreak
            1,                  // CLIP
            0,                  // read_idx - 临时占位
            0,                  // mapq - 临时占位
            false               // strand - 临时占位
        });
    }

    // ===== 2. Large Insertion 信号 =====
    current_pos = read.pos;
    read_consumed = 0;
    for (const auto& [op, len] : read.cigar_ops) {
        if (op == 'I' && len >= config_.min_ins_len) {
            // 插入信号：当前位置 + 上游消耗的 read 长度 = 断点位置
            anchors.push_back({
                chrom_tid,
                current_pos,      // 当前参考位置
                3,               // Unknown orientation (需要根据上下文推断)
                2,               // INS
                0,               // read_idx - 临时占位
                0,               // mapq - 临时占位
                false            // strand - 临时占位
            });
        }
        // 累加参考消耗以追踪位置
        if (op == 'M' || op == 'D' || op == 'N' || op == '=' || op == 'X') {
            current_pos += len;
            read_consumed += len;
        } else if (op == 'I' || op == 'S' || op == 'H' || op == '=') {
            read_consumed += len;
        }
    }

    // ===== 3. SA Tag 信号 =====
    // 遍历所有 SA split，每个都生成 Anchor
    if (read.has_sa) {
        for (const auto& [sa_tid, sa_pos] : read.sa_targets) {
            // 根据 SA 的位置和 Read 方向综合判断断点类型
            int orient = 3;  // Unknown

            // 根据 SA 位置与 primary alignment 的相对关系判断
            if (sa_pos < read.pos) {
                orient = 0;  // Left Break (比对到上游)
            } else if (sa_pos > read.end_pos) {
                orient = 1;  // Right Break (比对到下游)
            } else {
                // SA 位置在 primary 范围内，可能是复杂重排
                orient = 2;  // Mixed
            }

            anchors.push_back({
                sa_tid,       // SA 可能指向不同染色体（易位场景）
                sa_pos,
                orient,       // 根据位置判断的方向
                0,            // SA
                0,            // read_idx - 临时占位
                0,            // mapq - 临时占位
                false         // strand - 临时占位
            });
        }
    }

    return anchors;
}

void ComponentBuilder::merge_close_anchors(std::vector<Anchor>& anchors) const {
    if (anchors.size() <= 1) return;

    // 按 pos 排序
    std::sort(anchors.begin(), anchors.end(),
        [](const Anchor& a, const Anchor& b) { return a.pos < b.pos; });

    // 合并距离过近的 anchors
    std::vector<Anchor> merged;
    merged.reserve(anchors.size());

    for (const auto& curr : anchors) {
        if (merged.empty()) {
            merged.push_back(curr);
            continue;
        }

        Anchor& last = merged.back();

        // 距离检查
        if (curr.pos - last.pos <= config_.anchor_merge_distance) {
            // 冲突解决：优先保留 SA (source_type=0)，因为精度最高
            // 如果精度相同，保留 MapQ 较高的
            if (curr.source_type == 0 && last.source_type != 0) {
                last = curr;
            } else if (curr.source_type == last.source_type && curr.mapq > last.mapq) {
                last = curr;
            }
            // 否则保留原来的（通常是 CLIP，位置更稳定）
        } else {
            merged.push_back(curr);
        }
    }

    anchors.swap(merged);
}

std::vector<ClusterRange> ComponentBuilder::density_clustering(AnchorSpan anchors) const {
    std::vector<ClusterRange> ranges;

    if (anchors.empty()) return ranges;

    ClusterRange current_range{0, 1};
    int32_t last_pos = anchors[0].pos;
    int32_t last_tid = anchors[0].chrom_tid;

    for (size_t i = 1; i < anchors.size(); ++i) {
        // 检查染色体切换
        if (anchors[i].chrom_tid != last_tid) {
            // 完成当前 range
            ranges.push_back(current_range);
            // 开始新 range
            current_range = {i, i + 1};
        }
        // 检查 gap
        else if (anchors[i].pos - last_pos <= config_.cluster_gap) {
            // 扩展当前 range
            current_range.end_idx = i + 1;
        } else {
            // Gap 过大，完成当前 range，开始新的
            ranges.push_back(current_range);
            current_range = {i, i + 1};
        }
        last_pos = anchors[i].pos;
        last_tid = anchors[i].chrom_tid;
    }

    // 添加最后一个 range
    ranges.push_back(current_range);

    return ranges;
}

Component ComponentBuilder::build_component(
    AnchorSpan anchors,
    const ClusterRange& range,
    int32_t chrom_tid,
    int32_t id) const {

    Component comp;
    comp.id = id;
    comp.chrom_tid = chrom_tid;

    if (range.empty()) {
        return comp;
    }

    const Anchor& first = anchors[range.start_idx];
    const Anchor& last = anchors[range.end_idx - 1];

    comp.start = first.pos;
    comp.end = last.pos;
    comp.centroid = (comp.start + comp.end) / 2;

    // 收集唯一的 read indices
    std::vector<size_t> unique_reads;
    unique_reads.reserve(range.size());

    // 密度计算
    int32_t span = comp.end - comp.start + 1;
    uint32_t anchor_count = range.size();

    // 跨链统计
    int forward_count = 0;
    int reverse_count = 0;

    for (size_t i = range.start_idx; i < range.end_idx; ++i) {
        const auto& anchor = anchors[i];
        comp.read_indices.push_back(anchor.read_idx);

        if (anchor.read_idx >= unique_reads.size() ||
            std::find(unique_reads.begin(), unique_reads.end(), anchor.read_idx) == unique_reads.end()) {
            unique_reads.push_back(anchor.read_idx);
        }

        // 跨链统计
        if (anchor.strand) forward_count++;
        else reverse_count++;
    }

    // 去重
    std::sort(comp.read_indices.begin(), comp.read_indices.end());
    comp.read_indices.erase(
        std::unique(comp.read_indices.begin(), comp.read_indices.end()),
        comp.read_indices.end());

    comp.anchor_count = anchor_count;
    comp.read_count = comp.read_indices.size();

    // 密度统计
    comp.density = span > 0 ? static_cast<double>(anchor_count) / span : 0.0;
    comp.unique_strand_count = (forward_count > 0 && reverse_count > 0) ? 2 : 1;

    return comp;
}

void ComponentBuilder::recursive_cluster_breaker(
    AnchorSpan anchors,
    const ClusterRange& range,
    std::vector<Component>& components,
    int32_t chrom_tid,
    int32_t& component_id,
    int depth) const {

    // 边界条件检查
    if (range.empty() || depth >= config_.max_recursive_depth) {
        Component comp = build_component(anchors, range, chrom_tid, component_id++);
        components.push_back(std::move(comp));
        return;
    }

    // 计算范围信息
    int32_t cluster_start = anchors[range.start_idx].pos;
    int32_t cluster_end = anchors[range.end_idx - 1].pos;
    int32_t span = cluster_end - cluster_start;

    // 密度过滤：如果密度低于阈值，跳过此组件
    double density = static_cast<double>(range.size()) / (span + 1);
    if (density < config_.min_density) {
        return;  // 低密度噪声，直接丢弃
    }

    // 如果跨度正常，直接构建
    if (span <= config_.max_cluster_span) {
        Component comp = build_component(anchors, range, chrom_tid, component_id++);
        components.push_back(std::move(comp));
        return;
    }

    // 触发 Breaker：尝试递归切分
    breakers_triggered_++;

    // 寻找最佳切点：遍历所有可能的位置，找最大 gap
    int32_t best_gap = 0;
    size_t best_split = range.start_idx;

    for (size_t i = range.start_idx + 1; i < range.end_idx; ++i) {
        int32_t gap = anchors[i].pos - anchors[i-1].pos;
        // 只在合理位置切：gap < span / 2 且 gap 足够大
        if (gap > best_gap && gap < span / 2) {
            best_gap = gap;
            best_split = i;
        }
    }

    if (best_gap > 0 && best_split > range.start_idx && best_split < range.end_idx) {
        // 计算左右子簇的大小比例
        size_t left_size = best_split - range.start_idx;
        size_t right_size = range.end_idx - best_split;

        double left_ratio = static_cast<double>(left_size) / range.size();
        double right_ratio = static_cast<double>(right_size) / range.size();

        // 只有当两个子簇都足够大时才切分
        if (left_ratio >= config_.min_split_ratio && right_ratio >= config_.min_split_ratio) {
            // 递归切分左右子簇
            ClusterRange left_range{range.start_idx, best_split};
            ClusterRange right_range{best_split, range.end_idx};

            recursive_cluster_breaker(anchors, left_range, components, chrom_tid, component_id, depth + 1);
            recursive_cluster_breaker(anchors, right_range, components, chrom_tid, component_id, depth + 1);
            return;
        }
    }

    // 无法有效切分，作为散乱组件保留
    Component comp = build_component(anchors, range, chrom_tid, component_id++);
    components.push_back(std::move(comp));
}

}  // namespace placer
