#include "component_builder.h"
#include <algorithm>
#include <cmath>
#include <execution>
#include <numeric>
#include <unordered_set>
#include <cassert>

namespace placer {

ComponentBuilder::ComponentBuilder(ComponentBuilderConfig config)
    : config_(std::move(config)) {}

std::vector<Component> ComponentBuilder::build(const TaskData& task_data) {
    // 从 reads 中提取 chrom_tid（假设所有 reads 在同一染色体）
    int32_t chrom_tid = 0;
    if (!task_data.reads.empty()) {
        chrom_tid = task_data.reads[0].tid;
    }
    return build(task_data.reads, chrom_tid);
}

std::vector<Component> ComponentBuilder::build(
    const std::vector<ReadSketch>& reads, int32_t chrom_tid) {

    total_anchors_ = 0;
    total_components_ = 0;
    breakers_triggered_ = 0;

    // ===== Step 1: 多重锚点提取 =====
    std::vector<Anchor> all_anchors;
    all_anchors.reserve(reads.size() * 4);

    for (size_t i = 0; i < reads.size(); ++i) {
        // MapQ 过滤
        if (reads[i].mapq < config_.min_mapq) {
            continue;
        }

        auto anchors = extract_anchors(reads[i], chrom_tid);

        for (auto& anchor : anchors) {
            anchor.read_idx = i;
            anchor.mapq = reads[i].mapq;
            // [修正] 从 BAM flag 提取链信息：0x10 = reverse strand
            anchor.strand = !(reads[i].flag & 0x10);
        }

        if (!anchors.empty()) {
            merge_close_anchors(anchors);
            all_anchors.insert(all_anchors.end(), anchors.begin(), anchors.end());
        }
    }

    total_anchors_ = all_anchors.size();
    if (total_anchors_ == 0) {
        return {};
    }

    // 全局排序：先按 tid，再按 pos
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

        // [修正] 最小锚点数过滤：太少的 cluster 直接跳过
        if (range.size() < static_cast<size_t>(config_.min_anchors_per_cluster)) {
            continue;
        }

        recursive_cluster_breaker(
            anchor_span, range, components, chrom_tid, component_id, 0);
    }

    total_components_ = components.size();

    // [修正] 按基因组位置排序输出
    std::sort(components.begin(), components.end(),
        [](const Component& a, const Component& b) {
            if (a.chrom_tid != b.chrom_tid) return a.chrom_tid < b.chrom_tid;
            return a.start < b.start;
        });

    return components;
}

// ============================================================================
// Extract Anchors（从单条 read 提取所有锚点）
// [修正] CIGAR 解析中 read_consumed 的累加逻辑
// [修正] 右端 clip 位置计算
// ============================================================================

std::vector<Anchor> ComponentBuilder::extract_anchors(
    const ReadSketch& read, int32_t chrom_tid) {

    std::vector<Anchor> anchors;

    if (read.cigar_ops.empty()) {
        // 没有 CIGAR 信息，只能用位置信息
        if (read.pos >= 0) {
            anchors.push_back({
                chrom_tid, read.pos, 3, 0, 0, 0, false
            });
        }
        return anchors;
    }

    // ===== 0. 完整 CIGAR 解析 =====
    int32_t ref_pos = read.pos;       // 当前参考坐标
    int32_t query_pos = 0;            // 当前 query 坐标

    // 记录第一个和最后一个 clip 操作
    int32_t first_clip_len = 0;
    int32_t last_clip_len = 0;
    char first_clip_op = 0;
    char last_clip_op = 0;

    // [修正] 先扫描一遍确定 clip 信息
    if (!read.cigar_ops.empty()) {
        auto [op0, len0] = read.cigar_ops.front();
        if (op0 == 'S' || op0 == 'H') {
            first_clip_len = len0;
            first_clip_op = op0;
        }

        if (read.cigar_ops.size() > 1) {
            auto [opN, lenN] = read.cigar_ops.back();
            if (opN == 'S' || opN == 'H') {
                last_clip_len = lenN;
                last_clip_op = opN;
            }
        }
    }

    // [修正] 计算参考消耗的总长度（用于右端 clip 位置）
    int32_t total_ref_consumed = 0;
    for (const auto& [op, len] : read.cigar_ops) {
        if (op == 'M' || op == 'D' || op == 'N' || op == '=' || op == 'X') {
            total_ref_consumed += len;
        }
    }
    int32_t alignment_end = read.pos + total_ref_consumed;

    // ===== 1. Clip 信号处理 =====
    // 左端 Clip → RightBreak (5' 端断点)
    if (first_clip_len >= config_.min_clip_len) {
        anchors.push_back({
            chrom_tid,
            read.pos,           // 比对起始位置 = 断点位置
            1,                  // RightBreak
            1,                  // CLIP source
            0, 0, false
        });
    }

    // 右端 Clip → LeftBreak (3' 端断点)
    if (last_clip_len >= config_.min_clip_len) {
        // [修正] 右端断点位置 = 比对结束位置（参考坐标）
        anchors.push_back({
            chrom_tid,
            alignment_end,      // 比对结束位置
            0,                  // LeftBreak
            1,                  // CLIP source
            0, 0, false
        });
    }

    // ===== 2. Large Insertion / Deletion 信号 =====
    ref_pos = read.pos;
    query_pos = 0;

    // [修正] 跳过开头的 clip
    size_t cigar_start = 0;
    size_t cigar_end = read.cigar_ops.size();
    if (first_clip_op == 'S' || first_clip_op == 'H') cigar_start = 1;
    if (last_clip_op == 'S' || last_clip_op == 'H') cigar_end--;

    for (size_t ci = cigar_start; ci < cigar_end; ++ci) {
        auto [op, len] = read.cigar_ops[ci];

        if (op == 'I' && len >= config_.min_ins_len) {
            // 大插入信号
            anchors.push_back({
                chrom_tid,
                ref_pos,          // 插入发生在当前参考位置
                3,                // Unknown orientation
                2,                // INS source
                0, 0, false
            });
        }

        if (op == 'D' && len >= config_.min_del_len) {
            // [新增] 大缺失信号
            anchors.push_back({
                chrom_tid,
                ref_pos,          // 缺失起始位置
                3,                // Unknown orientation
                3,                // DEL source
                0, 0, false
            });

            // 缺失的另一端也产生锚点
            if (len > config_.min_del_len * 2) {
                anchors.push_back({
                    chrom_tid,
                    ref_pos + len,    // 缺失结束位置
                    3,
                    3,                // DEL source
                    0, 0, false
                });
            }
        }

        // [修正] 正确累加坐标
        switch (op) {
            case 'M': case '=': case 'X':
                ref_pos += len;
                query_pos += len;
                break;
            case 'D': case 'N':
                ref_pos += len;
                // query 不动
                break;
            case 'I':
                query_pos += len;
                // ref 不动
                break;
            case 'S':
                query_pos += len;
                break;
            case 'H':
                // 不消耗任何坐标
                break;
            case 'P':
                // padding，不消耗
                break;
            default:
                break;
        }
    }

    // ===== 3. SA Tag 信号 =====
    if (read.has_sa) {
        for (const auto& [sa_tid, sa_pos] : read.sa_targets) {
            int orient = 3;  // Unknown

            if (sa_tid == chrom_tid) {
                // 同染色体：根据相对位置判断方向
                if (sa_pos < read.pos) {
                    orient = 0;  // Left Break
                } else if (sa_pos > alignment_end) {
                    orient = 1;  // Right Break
                } else {
                    orient = 2;  // Mixed (重叠区域)
                }
            }
            // 不同染色体：保持 Unknown (3)，可能是易位

            anchors.push_back({
                sa_tid,
                sa_pos,
                orient,
                0,            // SA source
                0, 0, false
            });

            // [新增] 如果 SA 指向不同染色体，在当前染色体也放一个锚点
            // 标记为易位候选
            if (sa_tid != chrom_tid) {
                // 在 primary alignment 的 clip 端放锚点
                int32_t local_pos = (first_clip_len >= last_clip_len) ?
                                    read.pos : alignment_end;
                anchors.push_back({
                    chrom_tid,
                    local_pos,
                    2,            // Mixed (易位)
                    0,            // SA source
                    0, 0, false
                });
            }
        }
    }

    return anchors;
}

// ============================================================================
// Merge Close Anchors（同一 read 内去重合并）
// [修正] 合并逻辑：考虑 source_type 优先级 + 位置加权
// ============================================================================

void ComponentBuilder::merge_close_anchors(std::vector<Anchor>& anchors) const {
    if (anchors.size() <= 1) return;

    // 按 (chrom_tid, pos) 排序
    std::sort(anchors.begin(), anchors.end(),
        [](const Anchor& a, const Anchor& b) {
            if (a.chrom_tid != b.chrom_tid) return a.chrom_tid < b.chrom_tid;
            return a.pos < b.pos;
        });

    std::vector<Anchor> merged;
    merged.reserve(anchors.size());

    for (const auto& curr : anchors) {
        if (merged.empty()) {
            merged.push_back(curr);
            continue;
        }

        Anchor& last = merged.back();

        // 不同染色体不合并
        if (curr.chrom_tid != last.chrom_tid) {
            merged.push_back(curr);
            continue;
        }

        // 距离检查
        if (curr.pos - last.pos <= config_.anchor_merge_distance) {
            // [修正] 合并策略：
            // source_type 优先级: SA(0) > CLIP(1) > INS(2) > DEL(3)
            // 同优先级时取 MapQ 更高的
            // 位置取加权平均

            bool curr_better = false;

            if (curr.source_type < last.source_type) {
                curr_better = true;
            } else if (curr.source_type == last.source_type) {
                curr_better = (curr.mapq > last.mapq);
            }

            if (curr_better) {
                // [修正] 保留更好的锚点，但位置取加权平均
                int32_t avg_pos = (last.pos + curr.pos) / 2;
                last = curr;
                last.pos = avg_pos;
            } else {
                // 保留原来的，但更新位置为加权平均
                last.pos = (last.pos + curr.pos) / 2;
            }

            // [修正] 合并 mapq：取最大值
            last.mapq = std::max(last.mapq, curr.mapq);
        } else {
            merged.push_back(curr);
        }
    }

    anchors.swap(merged);
}

// ============================================================================
// Density Clustering（基于 gap 的线性聚类）
// [修正] 增加最小 cluster 大小过滤
// ============================================================================

std::vector<ClusterRange> ComponentBuilder::density_clustering(
    AnchorSpan anchors) const {

    std::vector<ClusterRange> ranges;

    if (anchors.empty()) return ranges;

    ClusterRange current_range{0, 1};
    int32_t last_pos = anchors[0].pos;
    int32_t last_tid = anchors[0].chrom_tid;

    for (size_t i = 1; i < anchors.size(); ++i) {
        bool should_break = false;

        // 染色体切换：必须断开
        if (anchors[i].chrom_tid != last_tid) {
            should_break = true;
        }
        // Gap 过大：断开
        else if (anchors[i].pos - last_pos > config_.cluster_gap) {
            should_break = true;
        }

        if (should_break) {
            // [修正] 只保留足够大的 cluster
            if (current_range.size() >= static_cast<size_t>(
                    config_.min_anchors_per_cluster)) {
                ranges.push_back(current_range);
            }
            current_range = {i, i + 1};
        } else {
            current_range.end_idx = i + 1;
        }

        last_pos = anchors[i].pos;
        last_tid = anchors[i].chrom_tid;
    }

    // 添加最后一个 range
    if (current_range.size() >= static_cast<size_t>(
            config_.min_anchors_per_cluster)) {
        ranges.push_back(current_range);
    }

    return ranges;
}

// ============================================================================
// Build Component（从 anchor range 构建 Component）
// [修正] 正确的 unique read 统计 + 链平衡计算
// ============================================================================

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

    // [修正] centroid 用加权平均（高 mapq 的锚点权重更大）
    double weighted_sum = 0.0;
    double weight_total = 0.0;

    // 收集唯一 read indices（用 unordered_set 去重）
    std::unordered_set<size_t> unique_read_set;

    // 统计信息
    uint32_t anchor_count = 0;
    int forward_count = 0;
    int reverse_count = 0;

    // source type 统计
    int sa_count = 0;
    int clip_count = 0;
    int ins_count = 0;
    int del_count = 0;

    for (size_t i = range.start_idx; i < range.end_idx; ++i) {
        const auto& anchor = anchors[i];

        // 加权 centroid
        double w = std::max(1.0, static_cast<double>(anchor.mapq));
        weighted_sum += anchor.pos * w;
        weight_total += w;

        // 唯一 reads
        unique_read_set.insert(anchor.read_idx);

        // 链统计
        if (anchor.strand) forward_count++;
        else reverse_count++;

        // source type 统计
        switch (anchor.source_type) {
            case 0: sa_count++; break;
            case 1: clip_count++; break;
            case 2: ins_count++; break;
            case 3: del_count++; break;
            default: break;
        }

        anchor_count++;
    }

    // [修正] 加权 centroid
    comp.centroid = (weight_total > 0) ?
        static_cast<int32_t>(weighted_sum / weight_total) :
        (comp.start + comp.end) / 2;

    // 填充 read_indices（排序去重）
    comp.read_indices.assign(unique_read_set.begin(), unique_read_set.end());
    std::sort(comp.read_indices.begin(), comp.read_indices.end());

    comp.anchor_count = anchor_count;
    comp.read_count = comp.read_indices.size();

    // 密度统计
    int32_t span = comp.end - comp.start + 1;
    comp.density = span > 0 ?
        static_cast<double>(anchor_count) / span : 0.0;

    // [修正] 链平衡：用 min/max 比值
    comp.unique_strand_count = 0;
    if (forward_count > 0) comp.unique_strand_count++;
    if (reverse_count > 0) comp.unique_strand_count++;

    // [新增] 额外统计信息
    comp.sa_anchor_count = sa_count;
    comp.clip_anchor_count = clip_count;
    comp.ins_anchor_count = ins_count;

    // [新增] 证据强度评估
    // SA 证据最强（精确断点），CLIP 次之，INS/DEL 最弱
    comp.evidence_score = sa_count * 3.0 + clip_count * 2.0 +
                          ins_count * 1.0 + del_count * 1.0;

    return comp;
}

// ============================================================================
// Recursive Cluster Breaker
// [修正] 切分策略：使用密度谷检测而非简单最大 gap
// ============================================================================

void ComponentBuilder::recursive_cluster_breaker(
    AnchorSpan anchors,
    const ClusterRange& range,
    std::vector<Component>& components,
    int32_t chrom_tid,
    int32_t& component_id,
    int depth) const {

    // 边界条件
    if (range.empty()) return;

    if (depth >= config_.max_recursive_depth) {
        Component comp = build_component(
            anchors, range, chrom_tid, component_id++);
        if (comp.read_count >= static_cast<size_t>(config_.min_reads_per_component)) {
            components.push_back(std::move(comp));
        }
        return;
    }

    // 计算范围信息
    int32_t cluster_start = anchors[range.start_idx].pos;
    int32_t cluster_end = anchors[range.end_idx - 1].pos;
    int32_t span = cluster_end - cluster_start;

    // 密度过滤
    double density = static_cast<double>(range.size()) / std::max(span, int32_t(1));
    if (density < config_.min_density && range.size() < 5) {
        return;  // 低密度噪声，丢弃
    }

    // 跨度正常，直接构建
    if (span <= config_.max_cluster_span) {
        Component comp = build_component(
            anchors, range, chrom_tid, component_id++);
        if (comp.read_count >= static_cast<size_t>(config_.min_reads_per_component)) {
            components.push_back(std::move(comp));
        }
        return;
    }

    // ===== 触发 Breaker =====
    breakers_triggered_++;

    // [修正] 使用"密度谷"检测而非简单最大 gap
    // 策略：滑动窗口计算局部密度，找密度最低的位置作为切点
    size_t n_anchors = range.size();

    // 方法1：最大 gap（快速）
    int32_t best_gap = 0;
    size_t best_split_gap = range.start_idx;

    for (size_t i = range.start_idx + 1; i < range.end_idx; ++i) {
        int32_t gap = anchors[i].pos - anchors[i - 1].pos;
        if (gap > best_gap) {
            best_gap = gap;
            best_split_gap = i;
        }
    }

    // 方法2：密度谷（更精确）
    size_t best_split_density = range.start_idx;
    double min_local_density = std::numeric_limits<double>::max();

    if (n_anchors >= 10) {
        // 滑动窗口大小 = 总锚点数的 10%，至少 3
        size_t window = std::max(size_t(3), n_anchors / 10);

        for (size_t i = range.start_idx + window;
             i + window < range.end_idx; ++i) {

            // 计算以 i 为中心的局部密度
            int32_t local_span = anchors[i + window - 1].pos -
                                 anchors[i - window].pos;
            if (local_span <= 0) continue;

            double local_density = static_cast<double>(2 * window) / local_span;

            if (local_density < min_local_density) {
                min_local_density = local_density;
                best_split_density = i;
            }
        }
    }

    // [修正] 选择切点：优先密度谷，回退到最大 gap
    size_t best_split = best_split_gap;

    if (n_anchors >= 10 && min_local_density < density * 0.3) {
        // 密度谷足够深，使用密度谷切点
        best_split = best_split_density;
    }

    // 验证切点有效性
    if (best_split <= range.start_idx || best_split >= range.end_idx) {
        // 无法切分
        Component comp = build_component(
            anchors, range, chrom_tid, component_id++);
        if (comp.read_count >= static_cast<size_t>(config_.min_reads_per_component)) {
            components.push_back(std::move(comp));
        }
        return;
    }

    // 计算左右子簇大小比例
    size_t left_size = best_split - range.start_idx;
    size_t right_size = range.end_idx - best_split;

    double left_ratio = static_cast<double>(left_size) / n_anchors;
    double right_ratio = static_cast<double>(right_size) / n_anchors;

    // [修正] 切分条件：两个子簇都要足够大
    bool can_split = (left_ratio >= config_.min_split_ratio &&
                      right_ratio >= config_.min_split_ratio);

    // 额外条件：gap 要足够大（至少是平均间距的 3 倍）
    if (can_split && best_split == best_split_gap) {
        double avg_gap = static_cast<double>(span) / std::max(n_anchors - 1, size_t(1));
        can_split = (best_gap >= avg_gap * 3.0);
    }

    if (can_split) {
        ClusterRange left_range{range.start_idx, best_split};
        ClusterRange right_range{best_split, range.end_idx};

        recursive_cluster_breaker(
            anchors, left_range, components, chrom_tid,
            component_id, depth + 1);
        recursive_cluster_breaker(
            anchors, right_range, components, chrom_tid,
            component_id, depth + 1);
    } else {
        // 无法有效切分
        Component comp = build_component(
            anchors, range, chrom_tid, component_id++);
        if (comp.read_count >= static_cast<size_t>(config_.min_reads_per_component)) {
            components.push_back(std::move(comp));
        }
    }
}

}  // namespace placer
