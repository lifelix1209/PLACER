#ifndef PLACER_ASSEMBLY_H
#define PLACER_ASSEMBLY_H

#include "component_builder.h"
#include "local_realign.h"
#include <vector>
#include <string>
#include <cstdint>
#include <string_view>
#include <optional>
#include <array>
#include <limits>
#include <algorithm>
#include <cstring>
#include <queue>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <memory>

namespace placer {

// ============================================================================
// Forward Declarations
// ============================================================================

class POAArena;

// ============================================================================
// Configuration
// ============================================================================

struct AssemblyConfig {
    // POA 参数
    int min_reads_for_poa = 3;
    int max_reads_for_poa = 50;
    float min_consensus_identity = 0.8f;

    // 对齐参数 (Affine Gap Penalty)
    int8_t match = 2;
    int8_t mismatch = -3;
    int8_t gap_open = -5;
    int8_t gap_extend = -1;

    // 分段参数
    int flank_min_length = 100;
    int ins_min_length = 50;

    // 结构合并参数
    int max_components_per_rep = 10;
    float min_fingerprint_similarity = 0.9f;

    // POA 输出
    int max_output_paths = 2;

    // 采样
    bool use_stratified_sampling = true;
    float high_quality_ratio = 0.7f;

    // 内存管理
    size_t max_arena_nodes = 1000000;  // 最大节点数
    size_t max_path_length = 10000;   // 最大路径长度
};

// ============================================================================
// SmallVector-style Predecessor Storage（内联前驱，避免64位限制）
// ============================================================================

/**
 * 工业级前驱存储：
 * - 内联存储 4 个前驱（大多数节点入度 <= 2）
 * - 溢出时使用动态分配
 */
struct PredList {
    static constexpr int INLINE_CAP = 4;
    uint32_t inline_preds[INLINE_CAP] = {UINT32_MAX, UINT32_MAX, UINT32_MAX, UINT32_MAX};
    uint8_t inline_count = 0;
    std::vector<uint32_t>* overflow = nullptr;  // 溢出时使用

    ~PredList() {
        delete overflow;
    }

    void add(uint32_t pred) {
        // 检查是否已存在
        for (uint8_t i = 0; i < inline_count; ++i) {
            if (inline_preds[i] == pred) return;
        }
        if (overflow) {
            for (uint32_t p : *overflow) {
                if (p == pred) return;
            }
            overflow->push_back(pred);
            return;
        }

        if (inline_count < INLINE_CAP) {
            inline_preds[inline_count++] = pred;
        } else {
            overflow = new std::vector<uint32_t>();
            overflow->push_back(pred);
        }
    }

    bool contains(uint32_t pred) const {
        for (uint8_t i = 0; i < inline_count; ++i) {
            if (inline_preds[i] == pred) return true;
        }
        if (overflow) {
            for (uint32_t p : *overflow) {
                if (p == pred) return true;
            }
        }
        return false;
    }

    void clear() {
        inline_count = 0;
        delete overflow;
        overflow = nullptr;
    }

    bool empty() const {
        return inline_count == 0 && (overflow == nullptr || overflow->empty());
    }

    size_t size() const {
        return inline_count + (overflow ? overflow->size() : 0);
    }
};

// ============================================================================
// Graph Node with Predecessors
// ============================================================================

/**
 * 工业级 POA 节点
 * - SmallVector 前驱存储（无64位限制）
 * - 32字节对齐
 */
struct alignas(32) POANode {
    char base = 'N';
    uint32_t count = 0;
    float quality_sum = 0.0f;

    // 路径信息
    uint32_t first_out = UINT32_MAX;
    uint32_t next_sibling = UINT32_MAX;

    // 前驱信息（SmallVector 模式）
    PredList predecessors;

    // 循环检测标记（避免回溯死循环）
    uint32_t visit_token = 0;
    static constexpr uint32_t VISIT_TOKEN_MASK = 0xFFFFFFFE;

    // 拓扑信息
    uint32_t topo_order = UINT32_MAX;
    uint32_t indegree = 0;

    // 路径权重
    uint32_t path_weight = 0;
};

// ============================================================================
// Edge Storage（边列表）
// ============================================================================

struct Edge {
    uint32_t from;
    uint32_t to;
    uint32_t weight;

    bool operator<(const Edge& other) const {
        if (from != other.from) return from < other.from;
        return to < other.to;
    }
};

// ============================================================================
// Arena Allocator for POA Nodes
// ============================================================================

class POAArena {
public:
    explicit POAArena(size_t estimated_nodes = 4096);
    ~POAArena();
    POAArena(const POAArena&) = delete;
    POAArena& operator=(const POAArena&) = delete;
    POAArena(POAArena&&) noexcept;
    POAArena& operator=(POAArena&&) noexcept;

    uint32_t allocate_node();
    uint32_t allocate_node(char base);

    POANode* get_node(uint32_t idx);
    const POANode* get_node(uint32_t idx) const;

    void reset();

    size_t size() const { return node_count_; }
    size_t capacity() const { return nodes_.size(); }

    // 添加边
    bool add_edge(uint32_t from, uint32_t to);

    // 获取所有边
    const std::set<Edge>& get_edges() const { return edges_; }

    // 获取前 k 个拓扑节点
    std::vector<uint32_t> get_topo_order() const;

    // 获取节点的所有后继
    std::vector<uint32_t> get_successors(uint32_t node_idx) const;

    // 拓扑排序（ Kahn 算法）
    void topological_sort();

private:
    static constexpr size_t NODE_SIZE = sizeof(POANode);
    std::vector<POANode> nodes_;
    std::set<Edge> edges_;  // 边集合（用于去重）
    size_t node_count_ = 0;
};

// ============================================================================
// Banded DP Matrix（带状动态规划）- O(N*W) 内存复杂度
// ============================================================================

/**
 * 带状动态规划缓冲区
 *
 * 工业级核心优化：
 * - 内存复杂度 O(N×W) 而非 O(N×M)，W = bandwidth（典型值 16-64）
 * - 适用于高相似度序列（<5% 差异）的长读段比对
 * - 自动适应超长读段（100kb+）而不会崩溃
 *
 * 原理：在高质量比对中，query 和 graph 的位置偏移（offset）保持在小范围内
 * 只计算对角线附近的 band 区域即可
 */
class BandedDPBuffer {
public:
    explicit BandedDPBuffer(int bandwidth = 32) : bandwidth_(bandwidth) {
        // 确保带宽为偶数，便于 SIMD
        if (bandwidth_ % 2 != 0) bandwidth_++;
        band_size_ = bandwidth_ + 1;  // +1 用于 offset=0
    }

    // 重置矩阵
    void reset(int query_len, int graph_len) {
        query_len_ = query_len;
        graph_len_ = graph_len;

        // 动态分配（按需扩展，而非固定 4096²）
        size_t needed = static_cast<size_t>(query_len_) * band_size_;
        if (M_.size() < needed) {
            M_.resize(needed, INT_MIN / 4);
            X_.resize(needed, INT_MIN / 4);
            Y_.resize(needed, INT_MIN / 4);
        }
        std::fill(M_.begin(), M_.end(), INT_MIN / 4);
        std::fill(X_.begin(), X_.end(), INT_MIN / 4);
        std::fill(Y_.begin(), Y_.end(), INT_MIN / 4);
    }

    // 将 (i, j) 映射到 band 索引
    // offset = j - i, 范围 [-bandwidth/2, +bandwidth/2]
    inline int band_index(int i, int j) const {
        int offset = j - i + bandwidth_ / 2;
        if (offset < 0 || offset >= band_size_) return -1;
        return i * band_size_ + offset;
    }

    // 带边界检查的访问（简化版：确保调用者在 band 内）
    inline int& M(int i, int j) {
        int idx = band_index(i, j);
        return M_[idx >= 0 ? idx : 0];
    }
    inline int& X(int i, int j) {
        int idx = band_index(i, j);
        return X_[idx >= 0 ? idx : 0];
    }
    inline int& Y(int i, int j) {
        int idx = band_index(i, j);
        return Y_[idx >= 0 ? idx : 0];
    }

    inline const int& M(int i, int j) const {
        int idx = band_index(i, j);
        return M_[idx >= 0 ? idx : 0];
    }
    inline const int& X(int i, int j) const {
        int idx = band_index(i, j);
        return X_[idx >= 0 ? idx : 0];
    }
    inline const int& Y(int i, int j) const {
        int idx = band_index(i, j);
        return Y_[idx >= 0 ? idx : 0];
    }

    // 检查 (i, j) 是否在 band 内
    inline bool in_band(int i, int j) const {
        int offset = j - i;
        return offset >= -bandwidth_ / 2 && offset <= bandwidth_ / 2;
    }

    int bandwidth() const { return bandwidth_; }
    int band_size() const { return band_size_; }

    // 估算序列是否适合带状比对
    static bool is_suitable_for_banded(int qlen, int glen, float max_divergence = 0.05f) {
        int64_t diff = static_cast<int64_t>(qlen) - static_cast<int64_t>(glen);
        float divergence = std::abs(static_cast<float>(diff)) / std::max(qlen, glen);
        return divergence <= max_divergence * 2;  // 放宽到 10%
    }

private:
    int bandwidth_;
    int band_size_;
    int query_len_ = 0;
    int graph_len_ = 0;
    std::vector<int> M_, X_, Y_;
    static constexpr int invalid_cell_ = INT_MIN / 4;
};

// ============================================================================
// DSU (Disjoint Set Union)
// ============================================================================

struct DSU {
    std::vector<int> parent;
    std::vector<int> rank;

    explicit DSU(int n) : parent(n), rank(n, 0) {
        std::iota(parent.begin(), parent.end(), 0);
    }

    int find(int x) {
        if (parent[x] != x) {
            parent[x] = find(parent[x]);
        }
        return parent[x];
    }

    void unite(int x, int y) {
        int px = find(x);
        int py = find(y);
        if (px == py) return;
        if (rank[px] < rank[py]) {
            std::swap(px, py);
        }
        parent[py] = px;
        if (rank[px] == rank[py]) {
            rank[px]++;
        }
    }
};

// ============================================================================
// 2D R-Tree Spatial Index（真正的空间索引）
// ============================================================================

/**
 * 2D R-Tree 实现
 *
 * 用于 SV 聚类的工业级空间索引：
 * - O(log N) 范围查询，而非线性扫描
 * - 支持高密度区域（着丝粒附近）的高效查询
 * - 完整的插入、删除、查询接口
 */
struct RTreeNode {
    float min_x, min_y;  // 左下角
    float max_x, max_y;  // 右上角
    int data;            // 关联的 contig 索引，-1 表示内部节点
    RTreeNode* parent = nullptr;  // 父节点指针
    std::vector<RTreeNode*> children;

    RTreeNode(float min_x, float min_y, float max_x, float max_y, int data = -1)
        : min_x(min_x), min_y(min_y), max_x(max_x), max_y(max_y), data(data) {}

    float area() const {
        return (max_x - min_x) * (max_y - min_y);
    }

    bool is_leaf() const {
        return data >= 0;
    }
};

class RTree {
public:
    explicit RTree(int max_capacity = 16, int min_capacity = 4);
    ~RTree();

    // 禁止拷贝
    RTree(const RTree&) = delete;
    RTree& operator=(const RTree&) = delete;
    RTree(RTree&&) noexcept;
    RTree& operator=(RTree&&) noexcept;

    // 插入 (x1, y1) 到 (x2, y2) 的矩形
    void insert(float x1, float y1, float x2, float y2, int data);

    // 范围查询：返回与 [q_x1, q_y1] x [q_x2, q_y2] 重叠的所有数据点
    std::vector<int> range_query(float q_x1, float q_y1, float q_x2, float q_y2) const;

    bool empty() const { return root_ == nullptr || (root_->data == -1 && root_->children.empty()); }
    size_t size() const { return size_; }

    void clear();

private:
    static constexpr int MAX_CAPACITY = 16;
    static constexpr int MIN_CAPACITY = 4;

    RTreeNode* root_ = nullptr;
    size_t size_ = 0;

    // 辅助函数
    RTreeNode* choose_leaf(RTreeNode* node, float x1, float y1, float x2, float y2);
    void split_node(RTreeNode* node, RTreeNode*& new_node);
    void adjust_tree(RTreeNode* node, RTreeNode* sibling = nullptr);
    RTreeNode* find_leaf(RTreeNode* node, float x1, float y1, float x2, float y2, int data);
    void condense_tree(RTreeNode* node, std::vector<RTreeNode*>& reinsert_list);
    void query_recursive(RTreeNode* node, float q_x1, float q_y1, float q_x2, float q_y2,
                         std::vector<int>& results) const;
    void destroy_node(RTreeNode* node);
    static void get_mbr(RTreeNode* node, float& x1, float& y1, float& x2, float& y2);
    static float enlarged_area(RTreeNode* node, float x1, float y1, float x2, float y2);
    static int overlap_enlarged(RTreeNode* node, RTreeNode* child, float x1, float y1, float x2, float y2);
    static void recalc_mbr(RTreeNode* node);
};

// ============================================================================
// Translocation-aware Clustering（易位感知聚类）
// ============================================================================

/**
 * 跨染色体事件索引
 *
 * 用于处理易位（Translocation）等结构变异：
 * - 键为 (tid1, tid2) 的有序对
 * - 值为该染色体对的 2D R-Tree
 */

// 前向声明
struct StructuralFingerprint;

class TranslocationIndex {
public:
    // 获取或创建染色体对的索引
    RTree& get_or_create(int32_t tid1, int32_t tid2);

    // 插入 SV 候选
    void insert(const StructuralFingerprint& fp, int contig_idx);

    // 查询与指定 SV 重叠的所有候选
    std::vector<int> query(const StructuralFingerprint& fp) const;

    bool empty() const { return index_.empty(); }
    size_t size() const { return index_.size(); }

    void clear();

private:
    // (tid1, tid2) -> RTree，tid1 <= tid2
    std::map<std::pair<int32_t, int32_t>, std::unique_ptr<RTree>> index_;
};

// ============================================================================
// Heaviest Bundle Algorithm（最重束共识提取）
// ============================================================================

/**
 * Heaviest Bundle 结果
 *
 * 在 DAG 上找到最优路径：
 * - 边的权重 = min(count_u, count_v)
 * - 节点权重 = count
 * - 使用动态规划找最重路径
 */
struct BundlePath {
    std::vector<uint32_t> nodes;       // 路径上的节点序列
    std::vector<int32_t> offsets;      // 每个节点的偏移量
    int32_t start_offset = 0;
    int32_t end_offset = 0;
    int64_t total_weight = 0;
    std::string consensus;
};

/**
 * Heaviest Bundle 路径提取器
 *
 * 工业级共识提取：
 * - 在拓扑序上做动态规划
 * - 边权重 = min(count_u, count_v)，反映序列一致性
 * - 处理分叉图，输出最优路径
 */
class HeaviestBundleExtractor {
public:
    explicit HeaviestBundleExtractor(int bandwidth = 32) : bandwidth_(bandwidth) {}

    // 提取最重束路径
    BundlePath extract(POAArena& arena, uint32_t graph_start);

    // 提取多条候选路径（用于多态性分析）
    std::vector<BundlePath> extract_top_k(POAArena& arena, uint32_t graph_start, int k = 3);

private:
    int bandwidth_;

    // 内部状态
    std::vector<int64_t> dp_weight_;
    std::vector<int> dp_prev_;
    std::vector<int32_t> dp_offset_;

    void ensure_capacity(int n) {
        if (static_cast<int>(dp_weight_.size()) < n) {
            dp_weight_.resize(n, INT64_MIN / 4);
            dp_prev_.resize(n, -1);
            dp_offset_.resize(n, 0);
        }
    }
};

// ============================================================================
// Sequence View
// ============================================================================

struct SeqSpan {
    const char* data = nullptr;
    size_t length = 0;

    SeqSpan() = default;
    SeqSpan(const char* d, size_t l) : data(d), length(l) {}
    explicit SeqSpan(std::string_view sv) : data(sv.data()), length(sv.size()) {}

    bool empty() const { return length == 0 || data == nullptr; }
    char operator[](size_t i) const { return data[i]; }
    std::string_view sv() const { return std::string_view(data, length); }
};

// ============================================================================
// Structural Fingerprint
// ============================================================================

struct StructuralFingerprint {
    int32_t tid = -1;           // 染色体 ID（必须匹配才能聚类）
    int32_t breakpoint_l = -1;  // 左断点
    int32_t breakpoint_l_end = -1;
    int32_t breakpoint_r = -1;  // 右断点
    int32_t breakpoint_r_end = -1;
    int32_t te_family_id = -1;
    int8_t orientation = 0;
    int8_t trunc_level = 0;
    int32_t ins_length_min = 0;
    int32_t ins_length_max = 0;
    bool has_inversion = false;

    static constexpr int32_t BP_TOLERANCE = 20;
    static constexpr int32_t LEN_TOLERANCE = 50;

    uint64_t hash() const;
    bool matches(const StructuralFingerprint& other) const;
    static StructuralFingerprint from_contig(
        std::string_view contig,
        int32_t left_bp,
        int32_t right_bp,
        int32_t te_family,
        int8_t orient);

    // 双端匹配：检查右端点重叠（用于发现左端点很远但右端点重叠的 SV）
    bool right_overlaps(const StructuralFingerprint& other) const {
        if (tid != other.tid) return false;
        // 右端点重叠判定
        return !(breakpoint_r_end < other.breakpoint_r ||
                 other.breakpoint_r_end < breakpoint_r);
    }

    // 双端匹配：检查左端点重叠
    bool left_overlaps(const StructuralFingerprint& other) const {
        if (tid != other.tid) return false;
        return !(breakpoint_l_end < other.breakpoint_l ||
                 other.breakpoint_l_end < breakpoint_l);
    }
};

// ============================================================================
// Contig Structure
// ============================================================================

struct Contig {
    std::string sequence;
    int32_t left_breakpoint = -1;
    int32_t right_breakpoint = -1;
    int32_t te_family_id = -1;
    int8_t orientation = 0;
    int8_t trunc_level = 0;
    int32_t support_reads = 0;
    double consensus_quality = 0.0;
    StructuralFingerprint fingerprint;

    std::string up_flank_seq;
    std::string ins_seq;
    std::string down_flank_seq;
};

// ============================================================================
// Structural Representative
// ============================================================================

struct StructuralRepresentative {
    int32_t rep_id = -1;
    std::vector<int32_t> component_ids;
    std::vector<int32_t> contig_ids;

    StructuralFingerprint fingerprint;
    std::string rep_sequence;

    struct Polymorphism {
        int32_t position = 0;
        char ref_base = 'N';
        char alt_base = 'N';
        int count = 0;
        double frequency = 0.0;
    };
    std::vector<Polymorphism> poly_summary;

    std::vector<int> polya_lengths;
    double polya_mean = 0.0;
    double polya_std = 0.0;

    int total_reads = 0;
    double avg_quality = 0.0;

    int tier = 3;
    double placeability_score = 0.0;
};

// ============================================================================
// Bi-interval Index（双端索引）- 解决左端点很远但右端点重叠的 SV 聚类问题
// ============================================================================

struct IntervalNodeBi {
    int32_t low_l;      // 左端点 low
    int32_t high_l;     // 左端点 high
    int32_t low_r;      // 右端点 low
    int32_t high_r;     // 右端点 high
    int32_t contig_idx;
};

class BiIntervalIndex {
public:
    explicit BiIntervalIndex(std::vector<IntervalNodeBi> nodes);

    // 查询与指定 SV 双端重叠的所有候选
    void query(const StructuralFingerprint& fp, std::vector<int>& candidates) const;

    bool empty() const { return nodes_.empty(); }
    size_t size() const { return nodes_.size(); }

private:
    std::vector<IntervalNodeBi> nodes_;
    std::vector<int> sorted_by_l_;
    std::vector<int> sorted_by_r_;
};

// ============================================================================
// Phred Quality Weighted Consensus（质量加权共识）
// ============================================================================

struct QualityWeightedConsensus {
    std::string sequence;
    std::vector<double> base_quality;  // Phred 质量分数

    // 提取质量加权的共识序列
    static QualityWeightedConsensus extract(
        const POAArena& arena,
        uint32_t graph_start,
        const AssemblyConfig& config);
};

// ============================================================================
// Circular DNA Detection（环状 DNA 检测）
// ============================================================================

struct CircularDNAConfig {
    static constexpr int32_t MIN_CIRCULAR_LENGTH = 1000;   // 最小环长
    static constexpr int32_t MAX_CIRCULAR_LENGTH = 100000;  // 最大环长
    static constexpr int32_t BREAKPOINT_TOLERANCE = 50;   // 断点容忍度
    static constexpr double MIN_COVERAGE_RATIO = 0.3;     // 最小覆盖比例
};

bool detect_circular_dna(
    const std::vector<ReadSketch>& reads,
    int32_t left_breakpoint,
    int32_t right_breakpoint,
    const GenomeAccessor& genome,
    const CircularDNAConfig& config = CircularDNAConfig());

// ============================================================================
// Industrial-Grade POA Engine（真图对齐）
// ============================================================================

class POAEngine {
public:
    explicit POAEngine(const AssemblyConfig& config = AssemblyConfig());

    // 构建图（接收外部 arena）
    uint32_t build_graph(POAArena& arena, const std::vector<SeqSpan>& sequences);

    // 真图对齐（Smith-Waterman on DAG）
    int align_to_graph(POAArena& arena, uint32_t graph_start, SeqSpan sequence);

    // 添加对齐后的序列到图（修改拓扑）
    void add_aligned_sequence(POAArena& arena,
                              const std::vector<uint32_t>& topo_order,
                              SeqSpan sequence,
                              const std::vector<uint8_t>& traceback,
                              int best_i,
                              int best_j);

    // 提取共识（最重路径）
    std::string extract_consensus(POAArena& arena, uint32_t graph_start);

    // 提取多路径
    std::vector<std::string> extract_paths(POAArena& arena,
                                           uint32_t graph_start,
                                           int max_paths);

    void reset();

private:
    AssemblyConfig config_;

    // 带状 DP 缓冲区（O(N*W) 内存复杂度）
    BandedDPBuffer dp_;

    // =========================================================================
    // 核心算法：真图对齐
    // =========================================================================

    /**
     * 在图上执行 Smith-Waterman
     * 关键：处理多前驱
     */
    int graph_smith_waterman(
        POAArena& arena,
        const std::vector<uint32_t>& topo_order,
        SeqSpan sequence,
        std::vector<uint8_t>& traceback,
        int& best_j);

    /**
     * 根据 traceback 更新图拓扑
     * - Match: 增加计数
     * - Insertion: 创建新节点链
     * - Deletion: 添加跳过边
     */
    void update_graph_topology(POAArena& arena,
                              const std::vector<uint32_t>& topo_order,
                              SeqSpan sequence,
                              const std::vector<uint8_t>& traceback,
                              int best_i,
                              int best_j);

    /**
     * 找到或创建匹配节点
     */
    uint32_t find_or_create_node(POAArena& arena,
                                  uint32_t current,
                                  char base,
                                  uint32_t& next_available);

    /**
     * 递归提取路径
     */
    void extract_paths_recursive(
        POAArena& arena,
        uint32_t node_idx,
        std::string& current_path,
        std::vector<std::string>& results,
        int max_paths,
        std::unordered_set<uint32_t>& visited);
};

// ============================================================================
// Thread-Safe POA Context（线程安全的 POA 计算上下文）
// ============================================================================

struct POAContext {
    POAEngine engine;
    POAArena arena;

    void reset() {
        arena.reset();
        engine.reset();
    }
};

// ============================================================================
// Assembly Engine
// ============================================================================

class AssemblyEngine {
public:
    explicit AssemblyEngine(AssemblyConfig config = AssemblyConfig());

    std::vector<Contig> assemble_component(
        const Component& component,
        const std::vector<ReadSketch>& reads,
        const GenomeAccessor& genome);

    // 线程安全：返回 POAEngine 工厂方法，每个组件使用独立的 POAEngine
    std::vector<Contig> assemble_batch(
        std::vector<Component>& components,
        const std::vector<ReadSketch>& reads,
        const GenomeAccessor& genome);

    std::vector<StructuralRepresentative> collapse_structurally(
        std::vector<Contig>& contigs);

    const AssemblyConfig& config() const { return config_; }

private:
    AssemblyConfig config_;

    // 移除共享的 poa_engine_，改为工厂方法创建线程本地上下文

    std::array<std::vector<std::string>, 3> extract_segments(
        const Component& component,
        const std::vector<ReadSketch>& reads,
        const GenomeAccessor& genome);

    std::vector<size_t> stratified_sample(
        const Component& component,
        const std::vector<ReadSketch>& reads,
        size_t max_samples);

    StructuralFingerprint build_fingerprint(const Contig& contig) const;
    StructuralRepresentative merge_contigs(
        const std::vector<int>& indices,
        std::vector<Contig>& contigs);
    int extract_polya_length(std::string_view seq);
    std::string simple_majority_consensus(
        const std::vector<std::string>& sequences);

    // 线程安全的组件组装辅助函数
    std::vector<Contig> assemble_component_thread_safe(
        const Component& component,
        const std::vector<ReadSketch>& reads,
        const GenomeAccessor& genome,
        POAContext& ctx);
};

}  // namespace placer

#endif  // PLACER_ASSEMBLY_H
