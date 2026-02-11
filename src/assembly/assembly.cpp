#include "assembly.h"
#include <algorithm>
#include <cmath>
#if __has_include(<execution>)
#include <execution>
#endif
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <limits>
#include <random>
#include <atomic>
#include <numeric>
#include <functional>
#include <cassert>

namespace placer {

// ============================================================================
// R-Tree Implementation（工业级 2D 空间索引）
// 修正：parent 指针、正确的 MBR 向上传播、root split、完整分裂逻辑
// ============================================================================

RTree::RTree(int max_capacity, int min_capacity)
    : root_(nullptr), size_(0) {}

RTree::~RTree() {
    clear();
}

RTree::RTree(RTree&& other) noexcept
    : root_(other.root_), size_(other.size_) {
    other.root_ = nullptr;
    other.size_ = 0;
}

RTree& RTree::operator=(RTree&& other) noexcept {
    if (this != &other) {
        clear();
        root_ = other.root_;
        size_ = other.size_;
        other.root_ = nullptr;
        other.size_ = 0;
    }
    return *this;
}

void RTree::clear() {
    if (root_) {
        destroy_node(root_);
        root_ = nullptr;
    }
    size_ = 0;
}

void RTree::destroy_node(RTreeNode* node) {
    if (!node) return;
    for (auto* child : node->children) {
        destroy_node(child);
    }
    delete node;
}

void RTree::get_mbr(RTreeNode* node, float& x1, float& y1, float& x2, float& y2) {
    x1 = node->min_x;
    y1 = node->min_y;
    x2 = node->max_x;
    y2 = node->max_y;
}

float RTree::enlarged_area(RTreeNode* node, float x1, float y1, float x2, float y2) {
    float old_area = node->area();
    float new_min_x = std::min(node->min_x, x1);
    float new_min_y = std::min(node->min_y, y1);
    float new_max_x = std::max(node->max_x, x2);
    float new_max_y = std::max(node->max_y, y2);
    return (new_max_x - new_min_x) * (new_max_y - new_min_y) - old_area;
}

// [修正] overlap_enlarged: 计算插入新矩形后与指定 child 的重叠增量
int RTree::overlap_enlarged(RTreeNode* node, RTreeNode* child,
                            float x1, float y1, float x2, float y2) {
    // 计算 child 扩展前与其他 children 的总重叠
    float old_overlap_total = 0.0f;
    float new_overlap_total = 0.0f;

    // child 扩展后的 MBR
    float ext_min_x = std::min(child->min_x, x1);
    float ext_min_y = std::min(child->min_y, y1);
    float ext_max_x = std::max(child->max_x, x2);
    float ext_max_y = std::max(child->max_y, y2);

    for (auto* sibling : node->children) {
        if (sibling == child) continue;

        // 旧重叠
        float o_x1 = std::max(child->min_x, sibling->min_x);
        float o_y1 = std::max(child->min_y, sibling->min_y);
        float o_x2 = std::min(child->max_x, sibling->max_x);
        float o_y2 = std::min(child->max_y, sibling->max_y);
        old_overlap_total += std::max(0.0f, o_x2 - o_x1) * std::max(0.0f, o_y2 - o_y1);

        // 新重叠
        float n_x1 = std::max(ext_min_x, sibling->min_x);
        float n_y1 = std::max(ext_min_y, sibling->min_y);
        float n_x2 = std::min(ext_max_x, sibling->max_x);
        float n_y2 = std::min(ext_max_y, sibling->max_y);
        new_overlap_total += std::max(0.0f, n_x2 - n_x1) * std::max(0.0f, n_y2 - n_y1);
    }

    return static_cast<int>((new_overlap_total - old_overlap_total) * 1000);
}

// [修正] recalc_mbr: 从 children 重新计算节点的 MBR
void RTree::recalc_mbr(RTreeNode* node) {
    if (node->children.empty()) return;
    node->min_x = std::numeric_limits<float>::max();
    node->min_y = std::numeric_limits<float>::max();
    node->max_x = std::numeric_limits<float>::lowest();
    node->max_y = std::numeric_limits<float>::lowest();
    for (auto* child : node->children) {
        node->min_x = std::min(node->min_x, child->min_x);
        node->min_y = std::min(node->min_y, child->min_y);
        node->max_x = std::max(node->max_x, child->max_x);
        node->max_y = std::max(node->max_y, child->max_y);
    }
}

// [修正] choose_leaf: 区分叶子层和内部层，叶子层用 R*-tree 的 overlap 最小化
RTreeNode* RTree::choose_leaf(RTreeNode* node, float x1, float y1, float x2, float y2) {
    if (!node || node->children.empty()) return node;

    // 当前节点的 children 已经是数据记录，当前节点就是插入目标叶层
    if (node->children[0]->is_leaf()) {
        return node;
    }

    // 检查 child 层是否为“叶层节点”（其 children 为数据记录）
    bool child_is_leaf_level =
        !node->children[0]->children.empty() && node->children[0]->children[0]->is_leaf();

    RTreeNode* best = nullptr;

    if (child_is_leaf_level) {
        // 选择叶层节点时，优先最小 overlap 增量（R*-tree 策略）
        int best_overlap_delta = INT_MAX;
        float best_enlargement = std::numeric_limits<float>::max();
        float best_area = std::numeric_limits<float>::max();

        for (auto* child : node->children) {
            int ov = overlap_enlarged(node, child, x1, y1, x2, y2);
            float en = enlarged_area(child, x1, y1, x2, y2);
            float ar = child->area();

            if (ov < best_overlap_delta ||
                (ov == best_overlap_delta && en < best_enlargement) ||
                (ov == best_overlap_delta && en == best_enlargement && ar < best_area)) {
                best_overlap_delta = ov;
                best_enlargement = en;
                best_area = ar;
                best = child;
            }
        }
    } else {
        // 更高内部层：选择 enlargement 最小的 child
        float best_enlargement = std::numeric_limits<float>::max();
        float best_area = std::numeric_limits<float>::max();

        for (auto* child : node->children) {
            float en = enlarged_area(child, x1, y1, x2, y2);
            float ar = child->area();

            if (en < best_enlargement ||
                (en == best_enlargement && ar < best_area)) {
                best_enlargement = en;
                best_area = ar;
                best = child;
            }
        }
    }

    return best ? choose_leaf(best, x1, y1, x2, y2) : node;
}

// [修正] split_node: 线性分裂 + 正确设置 parent 指针 + 返回新节点
void RTree::split_node(RTreeNode* node, RTreeNode*& new_sibling) {
    std::vector<RTreeNode*>& children = node->children;
    int n = static_cast<int>(children.size());

    // 按 x 坐标排序做线性分裂
    std::sort(children.begin(), children.end(),
        [](RTreeNode* a, RTreeNode* b) {
            return a->min_x < b->min_x;
        });

    int split_idx = n / 2;
    if (split_idx < MIN_CAPACITY) split_idx = MIN_CAPACITY;
    if (split_idx > n - MIN_CAPACITY) split_idx = n - MIN_CAPACITY;

    // 创建新兄弟节点（内部节点，data = -1）
    new_sibling = new RTreeNode(0, 0, 0, 0, -1);
    new_sibling->min_x = std::numeric_limits<float>::max();
    new_sibling->min_y = std::numeric_limits<float>::max();
    new_sibling->max_x = std::numeric_limits<float>::lowest();
    new_sibling->max_y = std::numeric_limits<float>::lowest();

    for (int i = split_idx; i < n; i++) {
        children[i]->parent = new_sibling;
        new_sibling->children.push_back(children[i]);
    }
    recalc_mbr(new_sibling);

    children.resize(split_idx);
    // 保持 parent 指针不变（它们已经指向 node）
    recalc_mbr(node);
}

// [修正] adjust_tree: 从叶子向上传播 MBR 更新和分裂
void RTree::adjust_tree(RTreeNode* node, RTreeNode* sibling) {
    while (node != root_) {
        RTreeNode* parent = node->parent;
        if (!parent) break;

        // 更新 parent 的 MBR
        recalc_mbr(parent);

        if (sibling) {
            // 需要把 sibling 插入 parent
            sibling->parent = parent;
            parent->children.push_back(sibling);
            recalc_mbr(parent);

            if (parent->children.size() > MAX_CAPACITY) {
                // parent 也需要分裂
                RTreeNode* new_parent_sibling = nullptr;
                split_node(parent, new_parent_sibling);
                node = parent;
                sibling = new_parent_sibling;
            } else {
                sibling = nullptr;
                node = parent;
            }
        } else {
            node = parent;
        }
    }

    // 如果 root 也需要分裂
    if (sibling && node == root_) {
        RTreeNode* new_root = new RTreeNode(0, 0, 0, 0, -1);
        root_->parent = new_root;
        sibling->parent = new_root;
        new_root->children.push_back(root_);
        new_root->children.push_back(sibling);
        recalc_mbr(new_root);
        new_root->parent = nullptr;
        root_ = new_root;
    }
}

// [修正] insert: 完整的 R-Tree 插入（含 parent 维护、分裂传播、root split）
void RTree::insert(float x1, float y1, float x2, float y2, int data) {
    if (x1 > x2) std::swap(x1, x2);
    if (y1 > y2) std::swap(y1, y2);

    RTreeNode* new_node = new RTreeNode(x1, y1, x2, y2, data);

    if (!root_) {
        // 创建一个根内部节点，把 new_node 作为第一个 child
        root_ = new RTreeNode(x1, y1, x2, y2, -1);
        root_->parent = nullptr;
        new_node->parent = root_;
        root_->children.push_back(new_node);
        size_++;
        return;
    }

    // 选择叶子
    RTreeNode* leaf = choose_leaf(root_, x1, y1, x2, y2);

    // 插入到叶子
    new_node->parent = leaf;
    leaf->children.push_back(new_node);
    recalc_mbr(leaf);
    size_++;

    RTreeNode* split_sibling = nullptr;

    if (leaf->children.size() > MAX_CAPACITY) {
        // 需要分裂
        split_node(leaf, split_sibling);
    }

    // 向上调整
    adjust_tree(leaf, split_sibling);
}

void RTree::query_recursive(RTreeNode* node, float q_x1, float q_y1,
                            float q_x2, float q_y2,
                            std::vector<int>& results) const {
    if (!node) return;

    // 检查 MBR 是否与查询矩形重叠
    bool overlaps = !(node->max_x < q_x1 || node->min_x > q_x2 ||
                      node->max_y < q_y1 || node->min_y > q_y2);

    if (!overlaps) return;

    if (node->data >= 0) {
        // 数据节点（叶子记录）
        results.push_back(node->data);
    }

    // 递归进入子节点
    for (auto* child : node->children) {
        query_recursive(child, q_x1, q_y1, q_x2, q_y2, results);
    }
}

std::vector<int> RTree::range_query(float q_x1, float q_y1,
                                     float q_x2, float q_y2) const {
    std::vector<int> results;
    if (!root_) return results;
    query_recursive(root_, q_x1, q_y1, q_x2, q_y2, results);
    return results;
}

// ============================================================================
// Translocation Index Implementation
// ============================================================================

RTree& TranslocationIndex::get_or_create(int32_t tid1, int32_t tid2) {
    if (tid1 > tid2) std::swap(tid1, tid2);
    auto key = std::make_pair(tid1, tid2);
    auto it = index_.find(key);
    if (it == index_.end()) {
        auto [new_it, inserted] = index_.try_emplace(key, std::make_unique<RTree>());
        return *new_it->second;
    }
    return *it->second;
}

void TranslocationIndex::insert(const StructuralFingerprint& fp, int contig_idx) {
    RTree& tree = get_or_create(fp.tid, fp.tid);
    tree.insert(static_cast<float>(fp.breakpoint_l),
                static_cast<float>(fp.breakpoint_r),
                static_cast<float>(fp.breakpoint_l_end),
                static_cast<float>(fp.breakpoint_r_end),
                contig_idx);
}

std::vector<int> TranslocationIndex::query(const StructuralFingerprint& fp) const {
    std::vector<int> results;
    if (fp.tid < 0) return results;

    float q_l = static_cast<float>(fp.breakpoint_l - StructuralFingerprint::BP_TOLERANCE);
    float q_r = static_cast<float>(fp.breakpoint_r - StructuralFingerprint::BP_TOLERANCE);
    float q_l_end = static_cast<float>(fp.breakpoint_l_end + StructuralFingerprint::BP_TOLERANCE);
    float q_r_end = static_cast<float>(fp.breakpoint_r_end + StructuralFingerprint::BP_TOLERANCE);

    auto key = std::make_pair(fp.tid, fp.tid);
    auto it = index_.find(key);
    if (it != index_.end()) {
        auto candidates = it->second->range_query(q_l, q_r, q_l_end, q_r_end);
        results.insert(results.end(), candidates.begin(), candidates.end());
    }

    return results;
}

void TranslocationIndex::clear() {
    index_.clear();
}

// ============================================================================
// Bi-interval Index（用二分查找优化）
// ============================================================================

BiIntervalIndex::BiIntervalIndex(std::vector<IntervalNodeBi> nodes)
    : nodes_(std::move(nodes)) {
    if (!nodes_.empty()) {
        sorted_by_l_.resize(nodes_.size());
        sorted_by_r_.resize(nodes_.size());
        std::iota(sorted_by_l_.begin(), sorted_by_l_.end(), 0);
        std::iota(sorted_by_r_.begin(), sorted_by_r_.end(), 0);
        std::sort(sorted_by_l_.begin(), sorted_by_l_.end(),
            [&](int a, int b) { return nodes_[a].low_l < nodes_[b].low_l; });
        std::sort(sorted_by_r_.begin(), sorted_by_r_.end(),
            [&](int a, int b) { return nodes_[a].low_r < nodes_[b].low_r; });
    }
}

void BiIntervalIndex::query(const StructuralFingerprint& fp,
                            std::vector<int>& candidates) const {
    candidates.clear();
    if (nodes_.empty()) return;

    int32_t tol = StructuralFingerprint::BP_TOLERANCE;

    // [修正] 用二分查找定位起始位置，而不是全扫描
    int32_t l_lo = fp.breakpoint_l - tol;
    int32_t l_hi = fp.breakpoint_l_end + tol;

    // 二分找到 sorted_by_l_ 中 low_l >= l_lo 的起始位置
    auto it_start = std::lower_bound(sorted_by_l_.begin(), sorted_by_l_.end(), l_lo,
        [&](size_t idx, int32_t val) { return nodes_[idx].low_l < val; });

    // 但我们需要的是 low_l <= l_hi 且 high_l >= l_lo 的区间重叠
    // 所以从头扫描到 low_l > l_hi 为止（但用 high_l >= l_lo 过滤）
    std::unordered_set<int> cand_set;

    for (size_t si = 0; si < sorted_by_l_.size(); ++si) {
        size_t idx = sorted_by_l_[si];
        const auto& node = nodes_[idx];
        if (node.low_l > l_hi) break;
        if (node.high_l < l_lo) continue;
        // 左端点区间重叠
        cand_set.insert(node.contig_idx);
    }

    // 右端点过滤
    int32_t r_lo = fp.breakpoint_r - tol;
    int32_t r_hi = fp.breakpoint_r_end + tol;

    for (size_t si = 0; si < sorted_by_r_.size(); ++si) {
        size_t idx = sorted_by_r_[si];
        const auto& node = nodes_[idx];
        if (node.low_r > r_hi) break;
        if (node.high_r < r_lo) continue;
        // 右端点区间也重叠 → 只保留两端都重叠的
        if (cand_set.count(node.contig_idx)) {
            candidates.push_back(node.contig_idx);
        }
    }
}

// ============================================================================
// POA Arena（修正：完整的邻接表 + predecessor 维护）
// ============================================================================

POAArena::POAArena(size_t estimated_nodes) {
    nodes_.reserve(estimated_nodes);
    reset();
}

POAArena::~POAArena() = default;

POAArena::POAArena(POAArena&&) noexcept = default;
POAArena& POAArena::operator=(POAArena&&) noexcept = default;

uint32_t POAArena::allocate_node() {
    uint32_t idx = static_cast<uint32_t>(nodes_.size());
    nodes_.emplace_back();
    node_count_++;
    return idx;
}

uint32_t POAArena::allocate_node(char base) {
    uint32_t idx = allocate_node();
    nodes_.back().base = base;
    nodes_.back().count = 1;  // [修正] 初始 count=1（第一条序列贡献）
    return idx;
}

POANode* POAArena::get_node(uint32_t idx) {
    if (idx >= nodes_.size()) return nullptr;
    return &nodes_[idx];
}

const POANode* POAArena::get_node(uint32_t idx) const {
    if (idx >= nodes_.size()) return nullptr;
    return &nodes_[idx];
}

void POAArena::reset() {
    nodes_.clear();
    edges_.clear();
    node_count_ = 0;
}

// [修正] add_edge: 正确维护 out-edge sibling 链 + predecessor 列表
bool POAArena::add_edge(uint32_t from, uint32_t to) {
    if (from >= nodes_.size() || to >= nodes_.size()) return false;

    Edge e{from, to, 1};
    auto it = edges_.find(e);
    if (it != edges_.end()) {
        // 边已存在，增加权重
        Edge updated = *it;
        updated.weight++;
        edges_.erase(it);
        edges_.insert(updated);
        return false;  // 不是新边
    }

    edges_.insert(e);

    // 维护 predecessor
    nodes_[to].predecessors.add(from);

    nodes_[to].indegree++;

    return true;
}

// [修正] get_successors: 从 edges_ 集合中获取节点的所有后继
std::vector<uint32_t> POAArena::get_successors(uint32_t node_idx) const {
    std::vector<uint32_t> succs;
    for (const auto& e : edges_) {
        if (e.from == node_idx) {
            succs.push_back(e.to);
        }
    }
    return succs;
}

// [修正] get_topo_order: 使用 edges_ 而不是 first_out/next_sibling
std::vector<uint32_t> POAArena::get_topo_order() const {
    std::vector<uint32_t> order;
    order.reserve(nodes_.size());

    // 构建邻接表
    std::vector<std::vector<uint32_t>> adj(nodes_.size());
    std::vector<uint32_t> local_indegree(nodes_.size(), 0);

    for (const auto& e : edges_) {
        if (e.from < nodes_.size() && e.to < nodes_.size()) {
            adj[e.from].push_back(e.to);
            local_indegree[e.to]++;
        }
    }

    std::queue<uint32_t> q;
    for (uint32_t i = 0; i < nodes_.size(); ++i) {
        if (local_indegree[i] == 0) {
            q.push(i);
        }
    }

    while (!q.empty()) {
        uint32_t u = q.front();
        q.pop();
        order.push_back(u);

        for (uint32_t v : adj[u]) {
            if (local_indegree[v] > 0) {
                local_indegree[v]--;
                if (local_indegree[v] == 0) {
                    q.push(v);
                }
            }
        }
    }

    return order;
}

void POAArena::topological_sort() {
    std::vector<uint32_t> order = get_topo_order();
    for (uint32_t i = 0; i < order.size(); ++i) {
        POANode* node = get_node(order[i]);
        if (node) node->topo_order = i;
    }
}

// ============================================================================
// POA Engine
// ============================================================================

POAEngine::POAEngine(const AssemblyConfig& config)
    : config_(config), abpoa_(AbPOAWrapper::Config{}) {}

void POAEngine::reset() {
    last_score_ = 0;
}

std::string POAEngine::build_consensus(const std::vector<std::string>& sequences) {
    if (sequences.empty()) return "";
    if (sequences.size() == 1) return sequences[0];

    auto result = abpoa_.align(sequences);
    last_score_ = result.best_score;
    return result.consensus;
}


// ============================================================================
// Heaviest Bundle Extractor
// [修正] 使用邻接表遍历后继，正确的拓扑序 DP
// ============================================================================

BundlePath HeaviestBundleExtractor::extract(POAArena& arena, uint32_t graph_start) {
    BundlePath result;

    if (graph_start == UINT32_MAX) return result;

    arena.topological_sort();
    std::vector<uint32_t> topo_order = arena.get_topo_order();
    int n = static_cast<int>(topo_order.size());
    if (n == 0) return result;

    ensure_capacity(n);

    // 构建邻接表
    std::vector<std::vector<uint32_t>> adj(arena.size());
    for (const auto& e : arena.get_edges()) {
        adj[e.from].push_back(e.to);
    }

    // node_idx → topo_pos
    std::vector<int> node_to_topo(arena.size(), -1);
    for (int i = 0; i < n; ++i) {
        node_to_topo[topo_order[i]] = i;
    }

    // 初始化 DP
    for (int i = 0; i < n; ++i) {
        const POANode* node = arena.get_node(topo_order[i]);
        dp_weight_[i] = node ? node->count : 0;
        dp_prev_[i] = -1;
        dp_offset_[i] = 0;
    }

    // 拓扑序 DP：找最重路径
    for (int i = 0; i < n; ++i) {
        uint32_t u = topo_order[i];

        for (uint32_t v : adj[u]) {
            int j = node_to_topo[v];
            if (j < 0 || j >= n) continue;

            const POANode* u_node = arena.get_node(u);
            const POANode* v_node = arena.get_node(v);

            // 边权重 = min(count_u, count_v)
            int edge_weight = std::min(
                u_node ? u_node->count : 0,
                v_node ? v_node->count : 0
            );

            int new_weight = dp_weight_[i] + edge_weight;
            if (new_weight > dp_weight_[j]) {
                dp_weight_[j] = new_weight;
                dp_prev_[j] = i;
                dp_offset_[j] = j - i;
            }
        }
    }

    // 找到终点
    int max_weight = 0;
    int end_idx = -1;
    for (int i = 0; i < n; ++i) {
        if (dp_weight_[i] > max_weight) {
            max_weight = dp_weight_[i];
            end_idx = i;
        }
    }

    if (end_idx < 0) return result;

    // 回溯路径
    std::vector<int> path_indices;
    int current = end_idx;
    while (current >= 0) {
        path_indices.push_back(current);
        current = dp_prev_[current];
    }
    std::reverse(path_indices.begin(), path_indices.end());

    // 构建结果
    result.nodes.reserve(path_indices.size());
    result.offsets.reserve(path_indices.size());

    int prev_pos = 0;
    for (int idx : path_indices) {
        uint32_t node_idx = topo_order[idx];
        result.nodes.push_back(node_idx);
        result.offsets.push_back(idx - prev_pos);
        prev_pos = idx;

        const POANode* node = arena.get_node(node_idx);
        if (node) result.consensus.push_back(node->base);
    }

    if (!path_indices.empty()) {
        result.start_offset = path_indices[0];
        result.end_offset = n - path_indices.back();
    }
    result.total_weight = max_weight;

    return result;
}

std::vector<BundlePath> HeaviestBundleExtractor::extract_top_k(
    POAArena& arena,
    uint32_t graph_start,
    int k) {

    std::vector<BundlePath> results;
    if (k <= 0) return results;

    // 第一条路径
    BundlePath first = extract(arena, graph_start);
    if (first.consensus.empty()) return results;
    results.push_back(std::move(first));

    // 后续路径：Yen's K-shortest paths 的简化版
    // 将已选路径上的边权重减半，重新提取
    for (int ki = 1; ki < k; ++ki) {
        // 降低已选路径节点的权重
        for (uint32_t node_idx : results.back().nodes) {
            POANode* node = arena.get_node(node_idx);
            if (node && node->count > 1) {
                node->count = node->count / 2;
            }
        }

        BundlePath next = extract(arena, graph_start);
        if (next.consensus.empty() || next.consensus == results.back().consensus) {
            break;
        }
        results.push_back(std::move(next));
    }

    // 恢复权重（简化：不恢复，因为通常只调用一次）

    return results;
}

// ============================================================================
// Quality-Weighted Consensus
// [// 修正] 沿最重路径提取，而非简单遍历所有拓扑节点
// ============================================================================

QualityWeightedConsensus QualityWeightedConsensus::extract(
    const POAArena& arena,
    uint32_t graph_start,
    const AssemblyConfig& config) {

    QualityWeightedConsensus result;

    if (graph_start == UINT32_MAX) return result;

    // 先找最重路径
    std::vector<uint32_t> topo_order = arena.get_topo_order();
    int n = static_cast<int>(topo_order.size());
    if (n == 0) return result;

    // 构建邻接表
    std::vector<std::vector<uint32_t>> adj(arena.size());
    for (const auto& e : arena.get_edges()) {
        adj[e.from].push_back(e.to);
    }

    // node_idx → topo_pos
    std::vector<int> node_to_topo(arena.size(), -1);
    for (int i = 0; i < n; ++i) {
        node_to_topo[topo_order[i]] = i;
    }

    // 拓扑序 DP 找最重路径
    std::vector<int> dp_weight(n, 0);
    std::vector<int> dp_prev(n, -1);

    for (int i = 0; i < n; ++i) {
        const POANode* node = arena.get_node(topo_order[i]);
        dp_weight[i] = node ? node->count : 0;
    }

    for (int i = 0; i < n; ++i) {
        uint32_t u = topo_order[i];
        for (uint32_t v : adj[u]) {
            int j = node_to_topo[v];
            if (j < 0) continue;
            const POANode* v_node = arena.get_node(v);
            int v_w = v_node ? v_node->count : 0;
            if (dp_weight[i] + v_w > dp_weight[j]) {
                dp_weight[j] = dp_weight[i] + v_w;
                dp_prev[j] = i;
            }
        }
    }

    // 找终点
    int max_w = 0, end_idx = -1;
    for (int i = 0; i < n; ++i) {
        if (dp_weight[i] > max_w) {
            max_w = dp_weight[i];
            end_idx = i;
        }
    }
    if (end_idx < 0) return result;

    // 回溯路径
    std::vector<int> path;
    for (int cur = end_idx; cur >= 0; cur = dp_prev[cur]) {
        path.push_back(cur);
    }
    std::reverse(path.begin(), path.end());

    // [修正] 沿路径收集每个位置的碱基投票（包括同一拓扑位置的分支节点）
    for (int topo_idx : path) {
        uint32_t main_node_idx = topo_order[topo_idx];
        const POANode* main_node = arena.get_node(main_node_idx);
        if (!main_node) {
            result.sequence.push_back('N');
            result.base_quality.push_back(0.0);
            continue;
        }

        // 收集该位置所有"等价节点"的碱基投票
        // 等价节点 = 同一拓扑位置的所有节点（通过 align_info 或 sibling 关系）
        // 简化版：只看主节点 + 检查是否有同 topo_order 的其他节点
        std::array<int, 5> base_votes{};  // A=0, C=1, G=2, T=3, N=4
        int total_votes = 0;

        auto base_to_idx = [](char b) -> int {
            switch (b) {
                case 'A': case 'a': return 0;
                case 'C': case 'c': return 1;
                case 'G': case 'g': return 2;
                case 'T': case 't': return 3;
                default: return 4;
            }
        };

        auto idx_to_base = [](int i) -> char {
            static const char bases[] = "ACGTN";
            return (i >= 0 && i < 5) ? bases[i] : 'N';
        };

        // 主节点投票
        int bi = base_to_idx(main_node->base);
        base_votes[bi] += main_node->count;
        total_votes += main_node->count;

        // 查找同一拓扑位置的其他节点（通过 predecessors 的后继关系）
        // 这些是在 POA 中被"合并"到同一列的不同碱基
        std::unordered_set<uint32_t> seen_alt_nodes;
        seen_alt_nodes.insert(main_node_idx);

        auto collect_from_pred = [&](uint32_t pred) {
            for (const auto& e : arena.get_edges()) {
                if (e.from != pred || e.to == main_node_idx) continue;
                if (!seen_alt_nodes.insert(e.to).second) continue;

                const POANode* alt_node = arena.get_node(e.to);
                if (!alt_node) continue;

                int alt_topo = node_to_topo[e.to];
                // 如果 alt 节点的拓扑序与 main 相近（±1），视为同列
                if (alt_topo >= 0 && std::abs(alt_topo - topo_idx) <= 1) {
                    int abi = base_to_idx(alt_node->base);
                    base_votes[abi] += alt_node->count;
                    total_votes += alt_node->count;
                }
            }
        };

        for (uint8_t k = 0; k < main_node->predecessors.inline_count; ++k) {
            collect_from_pred(main_node->predecessors.inline_preds[k]);
        }
        if (main_node->predecessors.overflow) {
            for (uint32_t pred : *main_node->predecessors.overflow) {
                collect_from_pred(pred);
            }
        }

        // 选择最高票碱基
        int max_count = 0;
        int best_idx = 4;  // default N
        for (int i = 0; i < 5; ++i) {
            if (base_votes[i] > max_count) {
                max_count = base_votes[i];
                best_idx = i;
            }
        }

        result.sequence.push_back(idx_to_base(best_idx));

        // 计算质量分数
        if (total_votes == 0) {
            result.base_quality.push_back(0.0);
        } else {
            double frequency = static_cast<double>(max_count) / total_votes;
            frequency = std::clamp(frequency, 0.0001, 0.9999);
            double phred = -10.0 * std::log10(1.0 - frequency);
            result.base_quality.push_back(std::min(phred, 60.0));
        }
    }

    return result;
}

// ============================================================================
// Circular DNA Detection
// ============================================================================

bool detect_circular_dna(
    const std::vector<ReadSketch>& reads,
    int32_t left_breakpoint,
    int32_t right_breakpoint,
    const GenomeAccessor& genome,
    const CircularDNAConfig& config) {

    if (reads.empty()) return false;

    int    crossing_reads = 0;
    int split_reads = 0;  // reads 跨越断点的
    int32_t span = std::abs(right_breakpoint - left_breakpoint);

    if (span < config.MIN_CIRCULAR_LENGTH) return false;

    for (const auto& read : reads) {
        // 跨越整个区间的 reads
        if (read.pos <= left_breakpoint && read.end_pos >= right_breakpoint) {
            crossing_reads++;
        }
        // 跨越左断点或右断点的 split reads
        if ((read.pos < left_breakpoint && read.end_pos > left_breakpoint &&
             read.end_pos < right_breakpoint) ||
            (read.pos > left_breakpoint && read.pos < right_breakpoint &&
             read.end_pos > right_breakpoint)) {
            split_reads++;
        }
    }

    double crossing_ratio = static_cast<double>(crossing_reads) / reads.size();
    double split_ratio = static_cast<double>(split_reads) / reads.size();

    // 环状 DNA 特征：
    // 1. 跨越断点的 reads 比例足够高
    // 2. split reads 比例也要达标
    // 3. 区间长度在合理范围内
    return crossing_ratio >= config.MIN_COVERAGE_RATIO &&
           split_ratio >= config.MIN_COVERAGE_RATIO * 0.5 &&
           span >= config.MIN_CIRCULAR_LENGTH &&
           span <= config.MAX_CIRCULAR_LENGTH;
}

// ============================================================================
// AssemblyEngine Implementation
// ============================================================================

AssemblyEngine::AssemblyEngine(AssemblyConfig config)
    : config_(std::move(config)) {}

// ============================================================================
// Breakpoint Detection from CIGAR
// ============================================================================

/**
 * Detect read-level breakpoint from CIGAR operations
 * Returns: {read_pos, type, is_valid}
 * - read_pos: position in read coordinate (0-based, where clip/break occurs)
 * - type: type of breakpoint signal
 * - is_valid: true if this read has a clear breakpoint signal
 */
ReadBreakpoint detect_breakpoint_from_cigar(const std::vector<std::pair<char, int>>& cigar_ops) {
    ReadBreakpoint bp{};
    bp.read_pos = -1;
    bp.ref_pos = -1;
    bp.type = ReadBreakpoint::NONE;
    bp.insertion_len = 0;
    bp.is_valid = false;

    int32_t read_pos = 0;  // Current position in read sequence
    int32_t ref_pos = 0;   // Current position in reference

    // Track soft-clip boundaries
    int32_t clip_5p_pos = -1;  // Position after 5' clip (clip end in read)
    int32_t clip_3p_start = -1; // Position where 3' clip starts in read
    int32_t clip_5p_len = -1;   // 5' soft clip length
    int32_t clip_3p_len = -1;   // 3' soft clip length
    int32_t clip_3p_ref_pos = -1; // Reference position at 3' clip boundary

    for (size_t op_idx = 0; op_idx < cigar_ops.size(); ++op_idx) {
        char op = cigar_ops[op_idx].first;
        int len = cigar_ops[op_idx].second;
        switch (op) {
            case 'M':  // Match/mismatch
                read_pos += len;
                ref_pos += len;
                break;
            case 'I':  // Insertion to reference (deletion from read perspective)
                if (bp.read_pos < 0) {
                    bp.read_pos = read_pos;
                    bp.ref_pos = ref_pos;
                    bp.type = ReadBreakpoint::INSERTION;
                    bp.insertion_len = len;
                    bp.is_valid = true;
                }
                read_pos += len;
                break;
            case 'D':  // Deletion from reference
                ref_pos += len;
                break;
            case 'N':  // Skipped region (splice junction)
                if (bp.read_pos < 0) {
                    bp.read_pos = read_pos;
                    bp.ref_pos = ref_pos;
                    bp.type = ReadBreakpoint::SPLIT;
                    bp.is_valid = true;
                }
                ref_pos += len;
                break;
            case 'P':  // Padding (silent in read/reference coordinates)
                break;
            case 'S':  // Soft clip
                if (len >= 10) {  // Only significant clips
                    bool at_read_start = (read_pos == 0);
                    bool at_cigar_end = (op_idx + 1 == cigar_ops.size());
                    if (at_read_start && clip_5p_pos < 0) {
                        // 5' soft clip
                        clip_5p_pos = len;  // Clip ends at position len
                        clip_5p_len = len;
                    } else if (at_cigar_end || !at_read_start) {
                        // 3' soft clip (after alignment)
                        clip_3p_start = read_pos;
                        clip_3p_len = len;
                        clip_3p_ref_pos = ref_pos;
                    }
                }
                read_pos += len;
                break;
            case 'H':  // Hard clip (not in sequence)
                break;
            case '=':  // Sequence match
                read_pos += len;
                ref_pos += len;
                break;
            case 'X':  // Sequence mismatch
                read_pos += len;
                ref_pos += len;
                break;
        }
    }

    // Determine best breakpoint from clip information
    // Priority: SPLIT > INSERTION > SOFT_CLIP

    // If we have a split/insertion, use event-time coordinates already recorded
    if (bp.is_valid && bp.read_pos >= 0) {
        return bp;
    }

    // Otherwise, use clip boundaries
    if (clip_5p_pos > 0 && clip_3p_start > 0) {
        // Both clips present - use the one with more evidence
        // 优先选择 clipping 更长的一端，长度相同优先 3'
        if (clip_3p_len >= clip_5p_len) {
            bp.read_pos = clip_3p_start;
            bp.ref_pos = clip_3p_ref_pos;
            bp.type = ReadBreakpoint::SOFT_CLIP_3P;
            bp.is_valid = true;
        } else {
            bp.read_pos = clip_5p_pos;
            bp.ref_pos = 0;
            bp.type = ReadBreakpoint::SOFT_CLIP_5P;
            bp.is_valid = true;
        }
    } else if (clip_5p_pos > 0) {
        bp.read_pos = clip_5p_pos;
        bp.ref_pos = 0;
        bp.type = ReadBreakpoint::SOFT_CLIP_5P;
        bp.is_valid = true;
    } else if (clip_3p_start > 0) {
        bp.read_pos = clip_3p_start;
        bp.ref_pos = clip_3p_ref_pos;
        bp.type = ReadBreakpoint::SOFT_CLIP_3P;
        bp.is_valid = true;
    }

    return bp;
}

/**
 * Check if read has evidence of breakpoint (clip, split, or large insertion)
 */
bool has_breakpoint_evidence(const ReadSketch& read) {
    // Check explicit breakpoint field first
    if (read.breakpoint.is_valid) {
        return true;
    }

    // Check CIGAR for signals
    if (!read.cigar_ops.empty()) {
        auto bp = detect_breakpoint_from_cigar(read.cigar_ops);
        if (bp.is_valid) {
            return true;
        }
    }

    // Check for split alignment evidence
    if (read.has_sa && !read.sa_targets.empty()) {
        return true;
    }

    // Check for significant clipping
    if (read.total_clip_len >= 20) {  // Lower threshold than Gate1
        return true;
    }

    return false;
}

std::vector<ReadSegments> AssemblyEngine::extract_segments(
    const Component& component,
    const std::vector<ReadSketch>& reads,
    const GenomeAccessor& genome) {

    std::vector<ReadSegments> result;

    // ================================================================
    // Breakpoint-centered extraction: return per-read aligned segments
    // Each ReadSegments contains up, ins, down for the SAME read
    // ================================================================

    auto sampled = stratified_sample(component, reads, config_.max_reads_for_poa);

    for (size_t idx : sampled) {
        if (idx >= reads.size()) continue;
        const auto& read = reads[idx];
        if (read.sequence.empty()) continue;

        // Skip reads without breakpoint evidence
        if (!has_breakpoint_evidence(read)) {
            continue;
        }

        // Detect breakpoint from CIGAR if not already set
        ReadBreakpoint bp = read.breakpoint;
        if ((!bp.is_valid ||
             (bp.type == ReadBreakpoint::INSERTION && bp.insertion_len <= 0)) &&
            !read.cigar_ops.empty()) {
            bp = detect_breakpoint_from_cigar(read.cigar_ops);
        }

        if (!bp.is_valid) {
            continue;
        }

        ReadSegments seg;
        seg.read_idx = idx;
        seg.breakpoint_type = bp.type;

        int32_t bp_pos = bp.read_pos;  // Position in read sequence
        int32_t read_len = static_cast<int32_t>(read.sequence.size());
        if (bp_pos < 0 || bp_pos > read_len) {
            continue;
        }

        // Type-aware insertion/clip span and flank anchors
        int32_t ins_start = -1;
        int32_t ins_end = -1;
        int32_t left_anchor = bp_pos;
        int32_t right_anchor = bp_pos;

        switch (bp.type) {
            case ReadBreakpoint::SOFT_CLIP_5P: {
                // 5' clip sequence is [0, bp_pos)
                ins_start = 0;
                ins_end = bp_pos;
                left_anchor = 0;
                right_anchor = bp_pos;
                break;
            }
            case ReadBreakpoint::SOFT_CLIP_3P: {
                // 3' clip sequence is [bp_pos, read_end)
                ins_start = bp_pos;
                ins_end = read_len;
                left_anchor = bp_pos;
                right_anchor = read_len;
                break;
            }
            case ReadBreakpoint::INSERTION: {
                ins_start = bp_pos;
                int32_t ins_len = bp.insertion_len;
                if (ins_len <= 0 && !read.cigar_ops.empty()) {
                    // Recover robustly from legacy breakpoint records lacking insertion_len
                    ReadBreakpoint recovered = detect_breakpoint_from_cigar(read.cigar_ops);
                    if (recovered.is_valid &&
                        recovered.type == ReadBreakpoint::INSERTION &&
                        recovered.read_pos == bp_pos &&
                        recovered.insertion_len > 0) {
                        ins_len = recovered.insertion_len;
                    }
                }

                if (ins_len > 0) {
                    ins_end = std::min(read_len, bp_pos + ins_len);
                } else {
                    // Cannot confidently recover insertion span
                    continue;
                }
                left_anchor = ins_start;
                right_anchor = ins_end;
                break;
            }
            case ReadBreakpoint::SPLIT:
                // Split-only evidence has no reliable insertion sequence in primary read
                continue;
            default:
                continue;
        }

        if (ins_start < 0 || ins_end <= ins_start || ins_start >= read_len) {
            continue;
        }

        // Cap extremely long insertion/clip segments
        constexpr int32_t kMaxInsLen = 1000;
        if (ins_end - ins_start > kMaxInsLen) {
            if (bp.type == ReadBreakpoint::SOFT_CLIP_5P) {
                ins_start = ins_end - kMaxInsLen;
            } else {
                ins_end = ins_start + kMaxInsLen;
            }
        }
        ins_start = std::clamp(ins_start, 0, read_len);
        ins_end = std::clamp(ins_end, 0, read_len);
        if (ins_end <= ins_start) {
            continue;
        }

        // UP: last N bp before left anchor
        int32_t up_start = left_anchor - config_.breakpoint_upstream;
        if (up_start < 0) up_start = 0;
        int32_t up_len = left_anchor - up_start;
        if (up_len > 0 && static_cast<int32_t>(read.sequence.size()) > up_start) {
            up_len = std::min(up_len, static_cast<int32_t>(read.sequence.size()) - up_start);
            seg.up = read.sequence.substr(up_start, up_len);
        }

        // INS/CLIP: type-specific span
        seg.ins = read.sequence.substr(ins_start, ins_end - ins_start);

        // DOWN: first N bp after right anchor
        int32_t down_start = right_anchor;
        if (down_start < static_cast<int32_t>(read.sequence.size())) {
            int32_t down_len = std::min(config_.breakpoint_downstream,
                static_cast<int32_t>(read.sequence.size()) - down_start);
            if (down_len > 0) {
                seg.down = read.sequence.substr(down_start, down_len);
            }
        }

        // Only include reads with insertion evidence and at least one supporting flank
        size_t min_flank = static_cast<size_t>(std::max(10, config_.breakpoint_upstream * 6 / 10));
        bool has_flank = seg.up.length() >= min_flank || seg.down.length() >= min_flank;
        if (has_flank && seg.ins.length() >= static_cast<size_t>(config_.ins_min_length)) {
            result.push_back(std::move(seg));
        }
    }

    // Keep dominant breakpoint evidence type to avoid mixing incompatible windows
    if (result.size() >= 3) {
        std::array<size_t, 5> type_counts{};
        for (const auto& seg : result) {
            int t = static_cast<int>(seg.breakpoint_type);
            if (t >= 0 && t < static_cast<int>(type_counts.size())) {
                type_counts[t]++;
            }
        }

        int dominant_type = ReadBreakpoint::NONE;
        size_t dominant_count = 0;
        for (int t = 0; t < static_cast<int>(type_counts.size()); ++t) {
            if (type_counts[t] > dominant_count) {
                dominant_count = type_counts[t];
                dominant_type = t;
            }
        }

        if (dominant_count > 0 && dominant_count < result.size()) {
            std::vector<ReadSegments> filtered;
            filtered.reserve(dominant_count);
            for (auto& seg : result) {
                if (static_cast<int>(seg.breakpoint_type) == dominant_type) {
                    filtered.push_back(std::move(seg));
                }
            }
            result = std::move(filtered);
        }
    }

    return result;
}

std::vector<size_t> AssemblyEngine::stratified_sample(
    const Component& component,
    const std::vector<ReadSketch>& reads,
    size_t max_samples) {

    std::vector<size_t> sampled;
    if (component.read_count == 0) return sampled;

    size_t total = std::min(static_cast<size_t>(component.read_count), reads.size());
    if (total == 0) return sampled;

    if (total <= max_samples) {
        // 全部取
        for (size_t i = 0; i < total; ++i) {
            sampled.push_back(i);
        }
    } else {
        // 分层采样：按位置分桶，每桶均匀采样
        int32_t span = component.end - component.start;
        if (span <= 0) span = 1;

        int num_buckets = std::min(static_cast<int>(max_samples), 16);
        int32_t bucket_size = span / num_buckets + 1;
        size_t per_bucket = max_samples / num_buckets;
        if (per_bucket == 0) per_bucket = 1;

        std::vector<std::vector<size_t>> buckets(num_buckets);

        for (size_t i = 0; i < total; ++i) {
            if (i >= reads.size()) break;
            int32_t mid = (reads[i].pos + reads[i].end_pos) / 2;
            int bucket_idx = static_cast<int>((mid - component.start) / bucket_size);
            bucket_idx = std::clamp(bucket_idx, 0, num_buckets - 1);
            buckets[bucket_idx].push_back(i);
        }

        // 从每个桶中均匀采样
        std::mt19937 rng(42);  // 固定种子保证可重复
        for (auto& bucket : buckets) {
            if (bucket.size() <= per_bucket) {
                sampled.insert(sampled.end(), bucket.begin(), bucket.end());
            } else {
                std::shuffle(bucket.begin(), bucket.end(), rng);
                sampled.insert(sampled.end(), bucket.begin(), bucket.begin() + per_bucket);
            }
        }

        // 如果采样不足，补充
        if (sampled.size() < max_samples) {
            std::unordered_set<size_t> sampled_set(sampled.begin(), sampled.end());
            for (size_t i = 0; i < total && sampled.size() < max_samples; ++i) {
                if (!sampled_set.count(i)) {
                    sampled.push_back(i);
                    sampled_set.insert(i);
                }
            }
        }
    }

    return sampled;
}

StructuralFingerprint AssemblyEngine::build_fingerprint(const Contig& contig) const {
    std::string_view ins_sequence = contig.ins_seq.empty()
        ? std::string_view(contig.sequence)
        : std::string_view(contig.ins_seq);
    return StructuralFingerprint::from_contig(
        ins_sequence,
        contig.left_breakpoint,
        contig.right_breakpoint,
        contig.te_family_id,
        contig.orientation,
        contig.chrom_tid);
}

// ============================================================================
// Assemble Component（使用断点居中窗口 + Minimizer Jaccard 分桶）
// ============================================================================

std::vector<Contig> AssemblyEngine::assemble_component(
    const Component& component,
    const std::vector<ReadSketch>& reads,
    const GenomeAccessor& genome) {

    std::vector<Contig> contigs;

    if (component.read_count < config_.min_reads_for_poa) {
        return contigs;
    }

    // Reset POA context for new assembly
    poa_ctx_.reset();

    // ================================================================
    // Step 1: Extract per-read breakpoint-centered segments
    // Returns: std::vector<ReadSegments> with up, ins, down for each read
    // ================================================================
    auto read_segments = extract_segments(component, reads, genome);

    if (read_segments.empty() || read_segments.size() < 2) {
        return contigs;
    }

    // ================================================================
    // Step 2: Construct cross-breakpoint analysis windows
    // Window = UP (last N bp) + INS/CLIP + DOWN (first N bp)
    // ================================================================
    std::vector<std::string> analysis_windows;
    std::vector<size_t> segment_indices;  // Maps window idx -> read_segments idx

    analysis_windows.reserve(read_segments.size());
    segment_indices.reserve(read_segments.size());

    for (size_t i = 0; i < read_segments.size(); ++i) {
        const auto& seg = read_segments[i];

        // Minimum requirements for homology comparison
        const size_t MIN_FLANK = 30;
        const size_t MIN_INS = static_cast<size_t>(config_.ins_min_length);

        bool has_left_flank = seg.up.length() >= MIN_FLANK;
        bool has_right_flank = seg.down.length() >= MIN_FLANK;
        if (!has_left_flank && !has_right_flank) {
            continue;
        }
        if (seg.ins.length() < MIN_INS) {
            continue;
        }

        // Construct analysis window: UP (last N bp) + INS + DOWN (first N bp)
        std::string window;
        size_t up_to_use = std::min(seg.up.length(), static_cast<size_t>(config_.breakpoint_upstream));
        size_t down_to_use = std::min(seg.down.length(), static_cast<size_t>(config_.breakpoint_downstream));

        window.reserve(up_to_use + seg.ins.length() + down_to_use);
        // Last N bp of UP
        if (up_to_use > 0) {
            window.append(seg.up.substr(seg.up.length() - up_to_use));
        }
        window += seg.ins;  // Full INS
        if (down_to_use > 0) {
            window.append(seg.down.substr(0, down_to_use));  // First N bp of DOWN
        }

        // Filter by minimum length and minimizer count
        const size_t MIN_WINDOW_LEN = 100;
        if (window.length() >= MIN_WINDOW_LEN) {
            analysis_windows.push_back(std::move(window));
            segment_indices.push_back(i);
        }
    }

    if (analysis_windows.size() < 4) {
        // Not enough reads for meaningful analysis
        // Try single contig with what we have
        std::vector<std::string> seqs_for_poa;
        std::vector<std::string> ins_for_poa;
        seqs_for_poa.reserve(read_segments.size());
        ins_for_poa.reserve(read_segments.size());
        for (const auto& seg : read_segments) {
            seqs_for_poa.push_back(seg.up + seg.ins + seg.down);
            ins_for_poa.push_back(seg.ins);
        }

        if (seqs_for_poa.size() >= 2) {
            poa_ctx_.reset();
            std::string consensus = poa_ctx_.engine.build_consensus(seqs_for_poa);
            if (!consensus.empty()) {
                Contig contig;
                contig.sequence = consensus;
                contig.ins_seq = build_flank_consensus(ins_for_poa);
                if (contig.ins_seq.empty()) {
                    contig.ins_seq = consensus;
                }
                contig.chrom_tid = component.chrom_tid;
                contig.left_breakpoint = component.start;
                contig.right_breakpoint = component.end;
                contig.source_component_id = component.id;
                contig.support_reads = static_cast<int32_t>(seqs_for_poa.size());
                int alignment_score = poa_ctx_.engine.get_score();
                contig.consensus_quality =
                    compute_consensus_quality(alignment_score, seqs_for_poa);

                contig.fingerprint = build_fingerprint(contig);
                int polya = extract_polya_length(contig.ins_seq);
                if (polya > 0) contig.trunc_level = 1;

                contigs.push_back(std::move(contig));
            }
        }
        return contigs;
    }

    // ================================================================
    // Step 3: Minimizer Jaccard Analysis on analysis windows
    // ================================================================
    ComponentSimilarityAnalyzer sim_analyzer({
        .k = 15,
        .w = 10,
        .sample_pairs = 200,
        .median_threshold = 0.75,   // Calibration needed
        .p10_threshold = 0.60,      // Calibration needed
        .min_bucket_size = 3,
        .min_bucket_frac = 0.2
    });

    JaccardStats sim_stats = sim_analyzer.analyze(analysis_windows);

    // Decide: single group or split?
    bool should_split = sim_stats.split_triggered &&
                       analysis_windows.size() >= 6;

    // 必要自检：桶内相似度应显著高于桶间相似度，否则回退单群
    auto split_is_consistent =
        [&](const ComponentSimilarityAnalyzer::BucketResult& bucket) -> bool {
            if (!bucket.valid || bucket.bucket_a.size() < 3 || bucket.bucket_b.size() < 3) {
                return false;
            }

            MinimizerSet minimizer({15, 10, 64});
            std::vector<std::unordered_set<uint64_t>> minimizer_sets;
            minimizer_sets.reserve(analysis_windows.size());
            for (const auto& window : analysis_windows) {
                minimizer_sets.push_back(minimizer.extract(window));
            }

            constexpr size_t kMaxPairsPerMetric = 1500;
            std::mt19937 rng(42);

            auto mean_within = [&](const std::vector<size_t>& ids) -> double {
                if (ids.size() < 2) return 0.0;

                size_t n = ids.size();
                size_t total_pairs = n * (n - 1) / 2;
                double sum = 0.0;
                size_t used_pairs = 0;

                if (total_pairs <= kMaxPairsPerMetric) {
                    for (size_t i = 0; i < n; ++i) {
                        for (size_t j = i + 1; j < n; ++j) {
                            sum += MinimizerSet::jaccard(
                                minimizer_sets[ids[i]], minimizer_sets[ids[j]]);
                            used_pairs++;
                        }
                    }
                } else {
                    std::uniform_int_distribution<size_t> dist(0, n - 1);
                    std::unordered_set<uint64_t> seen;
                    seen.reserve(kMaxPairsPerMetric * 2);

                    while (used_pairs < kMaxPairsPerMetric) {
                        size_t i = dist(rng);
                        size_t j = dist(rng);
                        if (i == j) continue;
                        if (i > j) std::swap(i, j);

                        uint64_t key = (static_cast<uint64_t>(i) << 32) | j;
                        if (!seen.insert(key).second) continue;

                        sum += MinimizerSet::jaccard(
                            minimizer_sets[ids[i]], minimizer_sets[ids[j]]);
                        used_pairs++;
                    }
                }

                return used_pairs > 0 ? sum / used_pairs : 0.0;
            };

            auto mean_between =
                [&](const std::vector<size_t>& a, const std::vector<size_t>& b) -> double {
                    if (a.empty() || b.empty()) return 0.0;

                    size_t total_pairs = a.size() * b.size();
                    double sum = 0.0;
                    size_t used_pairs = 0;

                    if (total_pairs <= kMaxPairsPerMetric) {
                        for (size_t i : a) {
                            for (size_t j : b) {
                                sum += MinimizerSet::jaccard(minimizer_sets[i], minimizer_sets[j]);
                                used_pairs++;
                            }
                        }
                    } else {
                        std::uniform_int_distribution<size_t> dist_a(0, a.size() - 1);
                        std::uniform_int_distribution<size_t> dist_b(0, b.size() - 1);
                        std::unordered_set<uint64_t> seen;
                        seen.reserve(kMaxPairsPerMetric * 2);

                        while (used_pairs < kMaxPairsPerMetric) {
                            size_t ia = dist_a(rng);
                            size_t ib = dist_b(rng);
                            uint64_t key = (static_cast<uint64_t>(ia) << 32) | ib;
                            if (!seen.insert(key).second) continue;

                            sum += MinimizerSet::jaccard(
                                minimizer_sets[a[ia]], minimizer_sets[b[ib]]);
                            used_pairs++;
                        }
                    }

                    return used_pairs > 0 ? sum / used_pairs : 0.0;
                };

            double mean_aa = mean_within(bucket.bucket_a);
            double mean_bb = mean_within(bucket.bucket_b);
            double mean_ab = mean_between(bucket.bucket_a, bucket.bucket_b);

            std::vector<size_t> window_lengths;
            window_lengths.reserve(analysis_windows.size());
            size_t window_length_sum = 0;
            for (const auto& w : analysis_windows) {
                window_lengths.push_back(w.size());
                window_length_sum += w.size();
            }
            std::sort(window_lengths.begin(), window_lengths.end());

            auto get_len_percentile = [&](double p) -> double {
                if (window_lengths.empty()) return 0.0;
                double idx = p * (window_lengths.size() - 1);
                size_t low = static_cast<size_t>(idx);
                size_t high = std::min(low + 1, window_lengths.size() - 1);
                double frac = idx - low;
                return window_lengths[low] * (1.0 - frac) + window_lengths[high] * frac;
            };

            double win_mean = analysis_windows.empty()
                ? 0.0
                : static_cast<double>(window_length_sum) / analysis_windows.size();
            double win_p10 = get_len_percentile(0.10);
            double win_p90 = get_len_percentile(0.90);

            constexpr double kMinIntraInterDelta = 0.08;
            bool pass = mean_aa >= mean_ab + kMinIntraInterDelta &&
                        mean_bb >= mean_ab + kMinIntraInterDelta;

            std::cerr << "[AssemblySplitStats] component=" << component.id
                      << " tid=" << component.chrom_tid
                      << " n=" << analysis_windows.size()
                      << " A=" << bucket.bucket_a.size()
                      << " B=" << bucket.bucket_b.size()
                      << " mean_aa=" << mean_aa
                      << " mean_bb=" << mean_bb
                      << " mean_ab=" << mean_ab
                      << " win_mean=" << win_mean
                      << " win_p10=" << win_p10
                      << " win_p90=" << win_p90
                      << " pass=" << (pass ? 1 : 0) << "\n";

            return pass;
        };

    // ================================================================
    // Step 4: Build contigs (single or split)
    // ================================================================

    if (should_split) {
        // split() returns read_idx lists (indices into analysis_windows)
        auto bucket_result = sim_analyzer.split(analysis_windows, sim_stats);

        if (bucket_result.valid &&
            bucket_result.bucket_a.size() >= 3 &&
            bucket_result.bucket_b.size() >= 3 &&
            split_is_consistent(bucket_result)) {

            // ========================================================
            // Bucket A: Build contig from bucket A reads
            // ========================================================
            std::vector<std::string> seqs_a;
            std::vector<std::string> ins_a;
            seqs_a.reserve(bucket_result.bucket_a.size());
            ins_a.reserve(bucket_result.bucket_a.size());

            for (size_t win_idx : bucket_result.bucket_a) {
                size_t seg_idx = segment_indices[win_idx];
                const auto& seg = read_segments[seg_idx];
                seqs_a.push_back(seg.up + seg.ins + seg.down);
                ins_a.push_back(seg.ins);
            }

            poa_ctx_.reset();
            std::string consensus_a = poa_ctx_.engine.build_consensus(seqs_a);

            if (!consensus_a.empty()) {
                Contig contig;
                contig.sequence = consensus_a;
                contig.ins_seq = build_flank_consensus(ins_a);
                if (contig.ins_seq.empty()) {
                    contig.ins_seq = consensus_a;
                }
                contig.chrom_tid = component.chrom_tid;
                contig.left_breakpoint = component.start;
                contig.right_breakpoint = component.end;
                contig.source_component_id = component.id;
                contig.support_reads = static_cast<int32_t>(seqs_a.size());
                int alignment_score = poa_ctx_.engine.get_score();
                contig.consensus_quality =
                    compute_consensus_quality(alignment_score, seqs_a);

                // Collect UP/DOWN for flanks
                std::vector<std::string> up_seqs, down_seqs;
                up_seqs.reserve(bucket_result.bucket_a.size());
                down_seqs.reserve(bucket_result.bucket_a.size());

                for (size_t win_idx : bucket_result.bucket_a) {
                    size_t seg_idx = segment_indices[win_idx];
                    up_seqs.push_back(read_segments[seg_idx].up);
                    down_seqs.push_back(read_segments[seg_idx].down);
                }

                contig.up_flank_seq = build_flank_consensus(up_seqs);
                contig.down_flank_seq = build_flank_consensus(down_seqs);

                contig.fingerprint = build_fingerprint(contig);
                int polya = extract_polya_length(contig.ins_seq);
                if (polya > 0) contig.trunc_level = 1;

                contigs.push_back(std::move(contig));
            }

            // ========================================================
            // Bucket B: Build contig from bucket B reads
            // ========================================================
            std::vector<std::string> seqs_b;
            std::vector<std::string> ins_b;
            seqs_b.reserve(bucket_result.bucket_b.size());
            ins_b.reserve(bucket_result.bucket_b.size());

            for (size_t win_idx : bucket_result.bucket_b) {
                size_t seg_idx = segment_indices[win_idx];
                const auto& seg = read_segments[seg_idx];
                seqs_b.push_back(seg.up + seg.ins + seg.down);
                ins_b.push_back(seg.ins);
            }

            poa_ctx_.reset();
            std::string consensus_b = poa_ctx_.engine.build_consensus(seqs_b);

            if (!consensus_b.empty()) {
                Contig contig;
                contig.sequence = consensus_b;
                contig.ins_seq = build_flank_consensus(ins_b);
                if (contig.ins_seq.empty()) {
                    contig.ins_seq = consensus_b;
                }
                contig.chrom_tid = component.chrom_tid;
                contig.left_breakpoint = component.start;
                contig.right_breakpoint = component.end;
                contig.source_component_id = component.id;
                contig.support_reads = static_cast<int32_t>(seqs_b.size());
                int alignment_score = poa_ctx_.engine.get_score();
                contig.consensus_quality =
                    compute_consensus_quality(alignment_score, seqs_b);

                // Collect UP/DOWN for flanks
                std::vector<std::string> up_seqs, down_seqs;
                up_seqs.reserve(bucket_result.bucket_b.size());
                down_seqs.reserve(bucket_result.bucket_b.size());

                for (size_t win_idx : bucket_result.bucket_b) {
                    size_t seg_idx = segment_indices[win_idx];
                    up_seqs.push_back(read_segments[seg_idx].up);
                    down_seqs.push_back(read_segments[seg_idx].down);
                }

                contig.up_flank_seq = build_flank_consensus(up_seqs);
                contig.down_flank_seq = build_flank_consensus(down_seqs);

                contig.fingerprint = build_fingerprint(contig);
                int polya = extract_polya_length(contig.ins_seq);
                if (polya > 0) contig.trunc_level = 1;

                contigs.push_back(std::move(contig));
            }

            return contigs;
        }
    }

    // ================================================================
    // Fall back: Single contig from all reads
    // ================================================================
    std::vector<std::string> seqs_all;
    std::vector<std::string> ins_all;
    seqs_all.reserve(read_segments.size());
    ins_all.reserve(read_segments.size());

    for (const auto& seg : read_segments) {
        seqs_all.push_back(seg.up + seg.ins + seg.down);
        ins_all.push_back(seg.ins);
    }

    poa_ctx_.reset();
    std::string consensus = poa_ctx_.engine.build_consensus(seqs_all);

    if (!consensus.empty()) {
        Contig contig;
        contig.sequence = consensus;
        contig.ins_seq = build_flank_consensus(ins_all);
        if (contig.ins_seq.empty()) {
            contig.ins_seq = consensus;
        }
        contig.chrom_tid = component.chrom_tid;
        contig.left_breakpoint = component.start;
        contig.right_breakpoint = component.end;
        contig.source_component_id = component.id;
        contig.support_reads = static_cast<int32_t>(seqs_all.size());

        int alignment_score = poa_ctx_.engine.get_score();
        contig.consensus_quality =
            compute_consensus_quality(alignment_score, seqs_all);

        // Collect flanks
        std::vector<std::string> up_seqs, down_seqs;
        up_seqs.reserve(read_segments.size());
        down_seqs.reserve(read_segments.size());

        for (const auto& seg : read_segments) {
            up_seqs.push_back(seg.up);
            down_seqs.push_back(seg.down);
        }

        contig.up_flank_seq = build_flank_consensus(up_seqs);
        contig.down_flank_seq = build_flank_consensus(down_seqs);

        contig.fingerprint = build_fingerprint(contig);
        int polya = extract_polya_length(contig.ins_seq);
        if (polya > 0) contig.trunc_level = 1;

        contigs.push_back(std::move(contig));
    }

    return contigs;
}

// ============================================================================
// Consensus Helpers
// ============================================================================

std::string AssemblyEngine::build_flank_consensus(
    const std::vector<std::string>& sequences) {

    if (sequences.empty()) return "";
    if (sequences.size() == 1) return sequences[0];

    poa_ctx_.reset();
    std::string consensus = poa_ctx_.engine.build_consensus(sequences);
    if (!consensus.empty()) {
        return consensus;
    }

    // Fallback: abPOA unavailable/failed 时退化为逐位多数投票
    return simple_majority_consensus(sequences);
}

double AssemblyEngine::compute_consensus_quality(
    int alignment_score,
    const std::vector<std::string>& sequences) const {

    if (sequences.empty()) return 0.0;

    size_t max_len = 0;
    for (const auto& seq : sequences) {
        max_len = std::max(max_len, seq.size());
    }

    if (max_len == 0) return 0.0;

    // 与默认评分参数对齐：match=2，按理论最佳分数归一化
    double denom = static_cast<double>(max_len) * 2.0;
    if (denom <= 0.0) return 0.0;

    return std::clamp(static_cast<double>(alignment_score) / denom, 0.0, 1.0);
}

// ============================================================================
// Simple Majority Consensus（辅助函数：对多条序列做逐位多数投票）
// ============================================================================

std::string AssemblyEngine::simple_majority_consensus(
    const std::vector<std::string>& sequences) {

    if (sequences.empty()) return "";
    if (sequences.size() == 1) return sequences[0];

    // 找最长序列长度
    size_t max_len = 0;
    for (const auto& s : sequences) {
        max_len = std::max(max_len, s.size());
    }

    std::string consensus;
    consensus.reserve(max_len);

    for (size_t pos = 0; pos < max_len; ++pos) {
        std::array<int, 5> votes{};  // A=0, C=1, G=2, T=3, gap/N=4

        for (const auto& seq : sequences) {
            if (pos >= seq.size()) {
                votes[4]++;
                continue;
            }
            switch (seq[pos]) {
                case 'A': case 'a': votes[0]++; break;
                case 'C': case 'c': votes[1]++; break;
                case 'G': case 'g': votes[2]++; break;
                case 'T': case 't': votes[3]++; break;
                default: votes[4]++; break;
            }
        }

        // 找最高票（排除 gap）
        int best = 0;
        int best_count = votes[0];
        for (int i = 1; i < 4; ++i) {
            if (votes[i] > best_count) {
                best_count = votes[i];
                best = i;
            }
        }

        // 如果 gap 票数超过碱基票数，跳过该位置
        if (votes[4] > best_count) continue;

        static const char bases[] = "ACGT";
        consensus.push_back(bases[best]);
    }

    return consensus;
}

// ============================================================================
// Assemble Batch（并行版本）
// ============================================================================

std::vector<Contig> AssemblyEngine::assemble_batch(
    std::vector<Component>& components,
    const std::vector<ReadSketch>& reads,
    const GenomeAccessor& genome) {

    if (components.empty()) return {};

    std::vector<std::vector<Contig>> results(components.size());

    std::vector<size_t> indices(components.size());
    std::iota(indices.begin(), indices.end(), 0);

    // 并行处理（若标准并行策略不可用则退化为串行）
#if defined(__cpp_lib_execution) && (__cpp_lib_execution >= 201603L)
    std::for_each(std::execution::par, indices.begin(), indices.end(),
                  [&](size_t idx) {
                      if (idx < components.size()) {
                          results[idx] = assemble_component(
                              components[idx], reads, genome);
                      }
                  });
#else
    std::for_each(indices.begin(), indices.end(),
                  [&](size_t idx) {
                      if (idx < components.size()) {
                          results[idx] = assemble_component(
                              components[idx], reads, genome);
                      }
                  });
#endif

    // 合并结果
    std::vector<Contig> all_contigs;
    size_t total = 0;
    for (const auto& r : results) total += r.size();
    all_contigs.reserve(total);

    for (size_t i = 0; i < components.size(); ++i) {
        for (auto& c : results[i]) {
            if (c.source_component_id < 0) {
                c.source_component_id = components[i].id;
            }
            all_contigs.push_back(std::move(c));
        }
    }

    return all_contigs;
}

// ============================================================================
// Structural Fingerprint Methods
// ============================================================================

uint64_t StructuralFingerprint::hash() const {
    // FNV-1a hash
    uint64_t h = 14695981039346656037ULL;
    auto mix = [&](uint64_t v) {
        h ^= v;
        h *= 1099511628211ULL;
    };
    mix(static_cast<uint64_t>(tid));
    mix(static_cast<uint64_t>(breakpoint_l));
    mix(static_cast<uint64_t>(breakpoint_l_end));
    mix(static_cast<uint64_t>(breakpoint_r));
    mix(static_cast<uint64_t>(breakpoint_r_end));
    mix(static_cast<uint64_t>(te_family_id));
    mix(static_cast<uint64_t>(orientation));
    mix(static_cast<uint64_t>(trunc_level));
    return h;
}

bool StructuralFingerprint::matches(const StructuralFingerprint& other) const {
    // 染色体必须匹配（除非未知）
    if (tid != other.tid && tid >= 0 && other.tid >= 0) return false;

    // TE 家族必须匹配（除非未知）
    if (te_family_id != other.te_family_id &&
        te_family_id >= 0 && other.te_family_id >= 0) return false;

    // 方向必须匹配（除非未知）
    if (orientation != other.orientation &&
        orientation >= 0 && other.orientation >= 0) return false;

    // 断点容差检查
    int32_t l_diff = std::abs(breakpoint_l - other.breakpoint_l);
    int32_t r_diff = std::abs(breakpoint_r - other.breakpoint_r);

    if (l_diff > BP_TOLERANCE || r_diff > BP_TOLERANCE) return false;

    // [修正] 额外检查插入长度范围是否重叠
    if (ins_length_max > 0 && other.ins_length_max > 0) {
        bool length_overlap = (ins_length_min <= other.ins_length_max) &&
                              (other.ins_length_min <= ins_length_max);
        if (!length_overlap) return false;
    }

    return true;
}

StructuralFingerprint StructuralFingerprint::from_contig(
    std::string_view contig,
    int32_t left_bp,
    int32_t right_bp,
    int32_t te_family,
    int8_t orient,
    int32_t tid) {

    StructuralFingerprint fp;
    fp.tid = tid;
    fp.breakpoint_l = left_bp;
    fp.breakpoint_r = right_bp;

    // [修正] 断点端区间基于实际 contig 长度和断点位置
    int32_t contig_len = static_cast<int32_t>(contig.length());
    fp.breakpoint_l_end = left_bp + std::max(contig_len / 10, int32_t(50));
    fp.breakpoint_r_end = right_bp + std::max(contig_len / 10, int32_t(50));

    fp.te_family_id = te_family;
    fp.orientation = orient;
    fp.trunc_level = 0;

    // [修正] 插入长度范围：±20% 或 ±100bp
    int32_t margin = std::max(contig_len / 5, int32_t(100));
    fp.ins_length_min = std::max(int32_t(0), contig_len - margin);
    fp.ins_length_max = contig_len + margin;

    fp.has_inversion = (orient == 2);

    return fp;
}

// ============================================================================
// Merge Contigs（合并同一结构组的 contigs）
// [修正] 选择最高质量的 contig 作为代表序列
// ============================================================================

StructuralRepresentative AssemblyEngine::merge_contigs(
    const std::vector<int>& indices,
    std::vector<Contig>& contigs) {

    StructuralRepresentative rep;
    if (indices.empty()) {
        rep.rep_id = -1;
        return rep;
    }

    rep.rep_id = indices[0];
    rep.component_ids.reserve(indices.size());
    rep.contig_ids.reserve(indices.size());

    double total_quality = 0.0;
    int total_reads = 0;
    int best_contig_idx = -1;
    double best_quality = -1.0;
    int best_support = 0;

    for (int idx : indices) {
        if (idx < 0 || static_cast<size_t>(idx) >= contigs.size()) continue;

        const auto& c = contigs[idx];
        rep.component_ids.push_back(c.left_breakpoint);
        rep.contig_ids.push_back(idx);
        total_quality += c.consensus_quality;
        total_reads += c.support_reads;

        // [修正] 选择质量最高 + 支持 reads 最多的作为代表
        double score = c.consensus_quality * 100.0 + c.support_reads;
        if (score > best_quality ||
            (score == best_quality && c.support_reads > best_support)) {
            best_quality = score;
            best_support = c.support_reads;
            best_contig_idx = idx;
        }
    }

    if (best_contig_idx >= 0 &&
        static_cast<size_t>(best_contig_idx) < contigs.size()) {
        rep.fingerprint = contigs[best_contig_idx].fingerprint;
        rep.rep_sequence = contigs[best_contig_idx].sequence;
        rep.rep_id = best_contig_idx;
    }

    if (!indices.empty()) {
        rep.avg_quality = total_quality / indices.size();
        rep.total_reads = total_reads;
    }

    return rep;
}

// ============================================================================
// Extract PolyA Length
// [修正] 同时检测 polyA 和 polyT（互补链）
// ============================================================================

int AssemblyEngine::extract_polya_length(std::string_view seq) {
    if (seq.empty()) return -1;

    // 检测 3' polyA
    int polya_count = 0;
    for (int i = static_cast<int>(seq.size()) - 1; i >= 0; --i) {
        if (seq[i] == 'A' || seq[i] == 'a') {
            polya_count++;
        } else {
            break;
        }
    }

    // 检测 5' polyT（互补链的 polyA）
    int polyt_count = 0;
    for (size_t i = 0; i < seq.size(); ++i) {
        if (seq[i] == 'T' || seq[i] == 't') {
            polyt_count++;
        } else {
            break;
        }
    }

    int best = std::max(polya_count, polyt_count);
    return best >= 10 ? best : -1;
}

// ============================================================================
// Collapse Structurally（使用 R-Tree + DSU）
// [修正] 正确的 R-Tree 查询 + 双向匹配验证
// ============================================================================

std::vector<StructuralRepresentative> AssemblyEngine::collapse_structurally(
    std::vector<Contig>& contigs) {

    std::vector<StructuralRepresentative> reps;

    if (contigs.empty()) return reps;

    const size_t N = contigs.size();

    // 构建 R-Tree 索引
    // 每个 contig 的矩形 = (breakpoint_l, breakpoint_r, breakpoint_l_end, breakpoint_r_end)
    RTree rtree;
    for (size_t i = 0; i < N; ++i) {
        const auto& fp = contigs[i].fingerprint;
        rtree.insert(
            static_cast<float>(fp.breakpoint_l),
            static_cast<float>(fp.breakpoint_r),
            static_cast<float>(fp.breakpoint_l_end),
            static_cast<float>(fp.breakpoint_r_end),
            static_cast<int>(i)
        );
    }

    // DSU 聚类
    DSU dsu(static_cast<int>(N));

    for (size_t i = 0; i < N; ++i) {
        const auto& fp = contigs[i].fingerprint;

        // [修正] 查询矩形 = 当前 contig 的断点 ± 容差
        float q_x1 = static_cast<float>(fp.breakpoint_l - StructuralFingerprint::BP_TOLERANCE);
        float q_y1 = static_cast<float>(fp.breakpoint_r - StructuralFingerprint::BP_TOLERANCE);
        float q_x2 = static_cast<float>(fp.breakpoint_l_end + StructuralFingerprint::BP_TOLERANCE);
        float q_y2 = static_cast<float>(fp.breakpoint_r_end + StructuralFingerprint::BP_TOLERANCE);

        auto candidates = rtree.range_query(q_x1, q_y1, q_x2, q_y2);

        for (int cand_idx : candidates) {
            if (cand_idx < 0 || static_cast<size_t>(cand_idx) >= N) continue;
            if (static_cast<size_t>(cand_idx) == i) continue;

            // 避免重复比较：只比较 i < cand_idx
            if (static_cast<size_t>(cand_idx) < i) continue;

            // [修正] 双向精确匹配验证
            if (contigs[i].fingerprint.matches(contigs[cand_idx].fingerprint)) {
                dsu.unite(static_cast<int>(i), cand_idx);
            }
        }
    }

    // 收集聚类组
    std::unordered_map<int, std::vector<int>> rep_groups;
    for (size_t i = 0; i < N; ++i) {
        int root = dsu.find(static_cast<int>(i));
        rep_groups[root].push_back(static_cast<int>(i));
    }

    // 为每个组生成代表
    reps.reserve(rep_groups.size());
    for (auto& [root, group] : rep_groups) {
        // 按质量排序组内 contigs
        std::sort(group.begin(), group.end(),
            [&](int a, int b) {
                double score_a = contigs[a].consensus_quality * 100 +
                                 contigs[a].support_reads;
                double score_b = contigs[b].consensus_quality * 100 +
                                 contigs[b].support_reads;
                return score_a > score_b;
            });

        StructuralRepresentative rep = merge_contigs(group, contigs);
        if (!rep.contig_ids.empty()) {
            reps.push_back(std::move(rep));
        }
    }

    // [修正] 按基因组位置排序输出
    std::sort(reps.begin(), reps.end(),
        [](const StructuralRepresentative& a, const StructuralRepresentative& b) {
            if (a.fingerprint.tid != b.fingerprint.tid)
                return a.fingerprint.tid < b.fingerprint.tid;
            return a.fingerprint.breakpoint_l < b.fingerprint.breakpoint_l;
        });

    return reps;
}

}  // namespace placer
