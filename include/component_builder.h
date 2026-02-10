#ifndef PLACER_COMPONENT_BUILDER_H
#define PLACER_COMPONENT_BUILDER_H

#include "bam_reader.h"
#include "task_queue.h"
#include <cstdint>
#include <vector>
#include <unordered_map>
#include <numeric>

namespace placer {

/**
 * Anchor: 断点证据的坐标抽象
 *
 * 物理意义：代表 Read 上的一个潜在断点位置
 * 核心原则：一条 Read 可能产生多个 Anchor（倒置/易位场景）
 */
struct Anchor {
    int32_t chrom_tid;      // 染色体 ID
    int32_t pos;            // 基因组坐标（断点位置）
    int source_type;        // 0=SA, 1=CLIP, 2=INS, 3=TE_HIT
    int orientation;        // 0=LeftBreak(3'端), 1=RightBreak(5'端), 2=Mixed, 3=Unknown
    size_t read_idx;       // 指向 TaskData.reads 的索引
    uint8_t mapq;          // 比对质量 (用于过滤低质量锚点)
    bool strand;           // 链信息：true=forward, false=reverse

    // 排序仅基于 (chrom_tid, pos)
    bool operator<(const Anchor& other) const {
        if (chrom_tid != other.chrom_tid) return chrom_tid < other.chrom_tid;
        return pos < other.pos;
    }

    bool operator==(const Anchor& other) const {
        return chrom_tid == other.chrom_tid && pos == other.pos;
    }
};

/**
 * LocusCandidate: 候选落点
 *
 * 物理意义：Component 内的候选插入位置
 */
struct LocusCandidate {
    int32_t chrom_tid;      // 染色体 ID
    int32_t pos;             // 候选位置 (断点中心)
    double score = 0.0;      // Placeability score
    uint32_t support_reads = 0;  // 支持的 reads 数量
    uint32_t evidence_mask = 0;  // 证据类型掩码 (bit0=SA, bit1=CLIP, bit2=INS)

    bool operator<(const LocusCandidate& other) const {
        return score > other.score;  // 降序排列
    }
};

/**
 * Component: 单断点簇
 *
 * 物理意义：代表一个物理上的断点事件
 * 核心原则：一个 Component 内的所有 Anchors 应该代表同一个断点
 */
struct Component {
    int32_t id;             // 组件唯一标识
    int32_t chrom_tid;      // 染色体 ID
    int32_t start;          // 聚类范围 start
    int32_t end;            // 聚类范围 end
    int32_t centroid;       // 重心坐标 (用于后续对齐)

    // 属于该断点的 reads (一个 read 可能出现在多个 component)
    std::vector<size_t> read_indices;

    // 候选落点集合 (用于 Phase 4 受限对齐)
    std::vector<LocusCandidate> locus_set;

    // 统计信息
    uint32_t anchor_count = 0;
    uint32_t read_count = 0;

    // 密度统计 (用于置信度评估)
    double density = 0.0;   // anchors per bp
    uint32_t unique_strand_count = 0;  // 跨链支持

    // 锚点来源统计
    uint32_t sa_anchor_count = 0;
    uint32_t clip_anchor_count = 0;
    uint32_t ins_anchor_count = 0;

    // 证据强度
    double evidence_score = 0.0;
};

/**
 * AnchorSpan: 轻量级范围视图 (C++17 compatible, 类似 std::span)
 */
class AnchorSpan {
public:
    AnchorSpan() : data_(nullptr), size_(0) {}
    AnchorSpan(Anchor* data, size_t size) : data_(data), size_(size) {}

    Anchor* data() const { return data_; }
    size_t size() const { return size_; }

    bool empty() const { return size_ == 0; }

    Anchor& operator[](size_t idx) const { return data_[idx]; }
    Anchor* begin() const { return data_; }
    Anchor* end() const { return data_ + size_; }

private:
    Anchor* data_;
    size_t size_;
};

/**
 * ClusterRange: 避免拷贝的轻量级范围表示
 */
struct ClusterRange {
    size_t start_idx;  // 在 all_anchors 中的起始索引
    size_t end_idx;    // 结束索引 (exclusive)

    bool empty() const { return start_idx >= end_idx; }
    size_t size() const { return end_idx - start_idx; }
};

/**
 * ComponentBuilderConfig: Component 构建配置
 */
struct ComponentBuilderConfig {
    // 锚点提取阈值
    int min_clip_len = 20;           // 最小 Clip 长度 (S/H)
    int min_ins_len = 50;           // 最小 Insertion 长度
    int min_del_len = 50;           // 最小 Deletion 长度

    // 聚类参数
    int cluster_gap = 50;            // 聚类间隔阈值 (bp)
    int max_cluster_span = 200;     // 最大聚类跨度 (bp) - 超过则触发 Breaker

    // Read 内去重参数
    int anchor_merge_distance = 20;  // 同一 read 内 Anchor 合并距离 (bp)

    // 质量过滤
    uint8_t min_mapq = 20;          // 最小 MapQ (0-255)
    double min_density = 0.01;      // 最小密度阈值 (anchors/bp)

    // SA tag 处理
    bool parse_all_sa_splits = true; // 解析所有 SA split 还是仅 primary

    // 递归拆分参数
    int max_recursive_depth = 4;    // 最大递归深度
    double min_split_ratio = 0.3;   // 最小切分比例 (子簇 / 原始簇)

    // 过滤参数
    int min_anchors_per_cluster = 3;  // 每个 cluster 最少锚点数
    int min_reads_per_component = 2;  // 每个 component 最少 reads 数
};

/**
 * ComponentBuilder: 基于锚点聚类的组件构建器
 *
 * 数据流：
 *   TaskData.reads (触发窗口内的所有 reads)
 *       ↓
 *   extract_anchors() - 多源锚点提取
 *       ↓
 *   density_clustering() - 基于 Gap 的聚类
 *       ↓
 *   cluster_breaker() - 跨度超限强制切分
 *       ↓
 *   std::vector<Component> - 构建好的组件列表
 */
class ComponentBuilder {
public:
    explicit ComponentBuilder(ComponentBuilderConfig config = ComponentBuilderConfig());

    /**
     * 从 TaskData 构建 Components
     *
     * @param task_data 任务数据（包含 reads 向量）
     * @return 构建的 Component 列表
     */
    std::vector<Component> build(const TaskData& task_data);

    /**
     * 从独立的 reads 向量构建 Components（用于测试）
     */
    std::vector<Component> build(const std::vector<ReadSketch>& reads, int32_t chrom_tid = 0);

    // 统计访问
    uint64_t get_total_anchors() const { return total_anchors_; }
    uint64_t get_total_components() const { return total_components_; }
    uint64_t get_breakers_triggered() const { return breakers_triggered_; }

    const ComponentBuilderConfig& config() const { return config_; }

private:
    /**
     * 从单条 Read 提取所有 Anchors
     * 原则：遍历所有信号源，生成可能的断点证据
     */
    std::vector<Anchor> extract_anchors(const ReadSketch& read, int32_t chrom_tid);

    /**
     * 合并同一 Read 内距离过近的 Anchors
     * 原则：同位置的多源信号合并，优先保留 SA（精度最高）
     */
    void merge_close_anchors(std::vector<Anchor>& anchors) const;

    /**
     * 基于 Gap 的密度聚类
     * 原则：Gap < threshold 归为一类
     */
    std::vector<ClusterRange> density_clustering(AnchorSpan anchors) const;

    /**
     * Cluster Breaker：跨度超限时强制切分（递归版）
     * 原则：span > max_span 时寻找低谷切开，支持多次递归切分
     */
    void recursive_cluster_breaker(
        AnchorSpan anchors,
        const ClusterRange& range,
        std::vector<Component>& components,
        int32_t chrom_tid,
        int32_t& component_id,
        int depth) const;

    /**
     * 根据单个范围构建 Component
     */
    Component build_component(
        AnchorSpan anchors,
        const ClusterRange& range,
        int32_t chrom_tid,
        int32_t id) const;

    /**
     * 计算断点方向
     * 根据 CIGAR 和 SA 信息判断 Left/Right break
     */
    static int determine_orientation(const ReadSketch& read, int anchor_pos);

    ComponentBuilderConfig config_;

    // 统计
    mutable uint64_t total_anchors_ = 0;
    mutable uint64_t total_components_ = 0;
    mutable uint64_t breakers_triggered_ = 0;
};

}  // namespace placer

#endif  // PLACER_COMPONENT_BUILDER_H
