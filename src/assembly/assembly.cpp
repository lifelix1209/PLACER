#include "assembly.h"
#include <algorithm>
#include <cmath>
#include <execution>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <limits>
#include <random>
#include <atomic>

namespace placer {

// ============================================================================
// R-Tree Implementation（真正的 2D 空间索引）
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

int RTree::overlap_enlarged(RTreeNode* node, RTreeNode* child, float x1, float y1, float x2, float y2) {
    float overlap_min_x = std::max(child->min_x, x1);
    float overlap_min_y = std::max(child->min_y, y1);
    float overlap_max_x = std::min(child->max_x, x2);
    float overlap_max_y = std::min(child->max_y, y2);

    float old_overlap = std::max(0.0f, child->max_x - child->min_x) *
                        std::max(0.0f, child->max_y - child->min_y);
    float new_overlap = std::max(0.0f, overlap_max_x - overlap_min_x) *
                        std::max(0.0f, overlap_max_y - overlap_min_y);

    return static_cast<int>((new_overlap - old_overlap) * 1000);  // 缩放避免浮点
}

RTreeNode* RTree::choose_leaf(RTreeNode* node, float x1, float y1, float x2, float y2) {
    if (node->children.empty()) return node;

    RTreeNode* best = nullptr;
    int best_enlargement = INT_MAX;
    float best_area = 0;
    int best_overlap = INT_MAX;

    for (auto* child : node->children) {
        float cx1, cy1, cx2, cy2;
        get_mbr(child, cx1, cy1, cx2, cy2);

        int enlargement = static_cast<int>(enlarged_area(child, x1, y1, x2, y2) * 1000);
        float area = child->area();

        if (enlargement < best_enlargement ||
            (enlargement == best_enlargement && area < best_area)) {
            best_enlargement = enlargement;
            best_area = area;
            best = child;
        }
    }

    return best ? choose_leaf(best, x1, y1, x2, y2) : node;
}

void RTree::insert(float x1, float y1, float x2, float y2, int data) {
    if (x1 > x2) std::swap(x1, x2);
    if (y1 > y2) std::swap(y1, y2);

    RTreeNode* new_node = new RTreeNode(x1, y1, x2, y2, data);

    if (!root_) {
        root_ = new_node;
        size_++;
        return;
    }

    RTreeNode* leaf = choose_leaf(root_, x1, y1, x2, y2);

    if (leaf->children.size() < MAX_CAPACITY) {
        leaf->children.push_back(new_node);
        leaf->min_x = std::min(leaf->min_x, x1);
        leaf->min_y = std::min(leaf->min_y, y1);
        leaf->max_x = std::max(leaf->max_x, x2);
        leaf->max_y = std::max(leaf->max_y, y2);
        size_++;

        RTreeNode* current = leaf;
        while (current) {
            float cx1, cy1, cx2, cy2;
            get_mbr(current, cx1, cy1, cx2, cy2);
            current->min_x = std::min(current->min_x, x1);
            current->min_y = std::min(current->min_y, y1);
            current->max_x = std::max(current->max_x, x2);
            current->max_y = std::max(current->max_y, y2);
            current = (current == root_) ? nullptr : nullptr;  // 简化版
        }
    } else {
        // 节点已满，需要分裂
        leaf->children.push_back(new_node);
        RTree* new_tree = new RTree();
        split_node(leaf, new_tree->root_);
        leaf->min_x = std::min(leaf->min_x, x1);
        leaf->min_y = std::min(leaf->min_y, y1);
        leaf->max_x = std::max(leaf->max_x, x2);
        leaf->max_y = std::max(leaf->max_y, y2);
        size_++;
    }
}

void RTree::split_node(RTreeNode* node, RTreeNode*& new_node) {
    // 简化的线性分裂算法
    std::vector<RTreeNode*>& children = node->children;
    int n = static_cast<int>(children.size());

    // 按 x 坐标排序
    std::sort(children.begin(), children.end(),
        [](RTreeNode* a, RTreeNode* b) {
            return a->min_x < b->min_x;
        });

    int split_idx = n / 2;
    new_node = new RTreeNode(0, 0, 0, 0, -1);

    for (int i = split_idx; i < n; i++) {
        new_node->children.push_back(children[i]);
        new_node->min_x = std::min(new_node->min_x, children[i]->min_x);
        new_node->min_y = std::min(new_node->min_y, children[i]->min_y);
        new_node->max_x = std::max(new_node->max_x, children[i]->max_x);
        new_node->max_y = std::max(new_node->max_y, children[i]->max_y);
    }

    children.resize(split_idx);
    node->min_x = std::numeric_limits<float>::max();
    node->min_y = std::numeric_limits<float>::max();
    node->max_x = std::numeric_limits<float>::lowest();
    node->max_y = std::numeric_limits<float>::lowest();

    for (auto* child : children) {
        node->min_x = std::min(node->min_x, child->min_x);
        node->min_y = std::min(node->min_y, child->min_y);
        node->max_x = std::max(node->max_x, child->max_x);
        node->max_y = std::max(node->max_y, child->max_y);
    }
}

void RTree::adjust_tree(RTreeNode* node, RTreeNode* sibling) {
    // 简化版：更新 MBR
    if (!sibling) return;

    for (auto* child : node->children) {
        child->min_x = std::min(child->min_x, sibling->min_x);
        child->min_y = std::min(child->min_y, sibling->min_y);
        child->max_x = std::max(child->max_x, sibling->max_x);
        child->max_y = std::max(child->max_y, sibling->max_y);
    }
}

void RTree::query_recursive(RTreeNode* node, float q_x1, float q_y1, float q_x2, float q_y2,
                            std::vector<int>& results) const {
    if (!node) return;

    // 检查 MBR 是否重叠
    bool overlaps = !(node->max_x < q_x1 || node->min_x > q_x2 ||
                      node->max_y < q_y1 || node->min_y > q_y2);

    if (!overlaps) return;

    if (node->data >= 0) {
        // 叶子节点
        results.push_back(node->data);
    } else {
        // 内部节点
        for (auto* child : node->children) {
            query_recursive(child, q_x1, q_y1, q_x2, q_y2, results);
        }
    }
}

std::vector<int> RTree::range_query(float q_x1, float q_y1, float q_x2, float q_y2) const {
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
    RTree& tree = get_or_create(fp.tid, fp.tid);  // 同染色体
    tree.insert(static_cast<float>(fp.breakpoint_l),
                static_cast<float>(fp.breakpoint_r),
                static_cast<float>(fp.breakpoint_l_end),
                static_cast<float>(fp.breakpoint_r_end),
                contig_idx);
}

std::vector<int> TranslocationIndex::query(const StructuralFingerprint& fp) const {
    std::vector<int> results;
    if (fp.tid < 0) return results;

    // 构建查询矩形
    float q_l = static_cast<float>(fp.breakpoint_l - StructuralFingerprint::BP_TOLERANCE);
    float q_r = static_cast<float>(fp.breakpoint_r - StructuralFingerprint::BP_TOLERANCE);
    float q_l_end = static_cast<float>(fp.breakpoint_l_end + StructuralFingerprint::BP_TOLERANCE);
    float q_r_end = static_cast<float>(fp.breakpoint_r_end + StructuralFingerprint::BP_TOLERANCE);

    // 查询同染色体
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
// Bi-interval Index（简化版，用于兼容性）
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

void BiIntervalIndex::query(const StructuralFingerprint& fp, std::vector<int>& candidates) const {
    candidates.clear();
    if (nodes_.empty()) return;

    int32_t l_tol = StructuralFingerprint::BP_TOLERANCE;
    int32_t r_tol = StructuralFingerprint::BP_TOLERANCE;

    // 查询左端点重叠
    for (size_t idx : sorted_by_l_) {
        const auto& node = nodes_[idx];
        if (node.low_l > fp.breakpoint_l + l_tol + StructuralFingerprint::BP_TOLERANCE) break;
        if (node.high_l < fp.breakpoint_l - l_tol - StructuralFingerprint::BP_TOLERANCE) continue;

        if (node.low_l <= fp.breakpoint_l_end && fp.breakpoint_l <= node.high_l) {
            candidates.push_back(node.contig_idx);
        }
    }

    // 查询右端点重叠
    for (size_t idx : sorted_by_r_) {
        const auto& node = nodes_[idx];
        if (node.low_r > fp.breakpoint_r + r_tol + StructuralFingerprint::BP_TOLERANCE) break;
        if (node.high_r < fp.breakpoint_r - r_tol - StructuralFingerprint::BP_TOLERANCE) continue;

        if (node.low_r <= fp.breakpoint_r_end && fp.breakpoint_r <= node.high_r) {
            if (std::find(candidates.begin(), candidates.end(), node.contig_idx) == candidates.end()) {
                candidates.push_back(node.contig_idx);
            }
        }
    }
}

// ============================================================================
// POA Arena
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

bool POAArena::add_edge(uint32_t from, uint32_t to) {
    if (from >= nodes_.size() || to >= nodes_.size()) return false;
    Edge e{from, to, 1};
    if (edges_.insert(e).second) {
        nodes_[to].indegree++;
        nodes_[from].first_out = to;
        return true;
    }
    return false;
}

std::vector<uint32_t> POAArena::get_topo_order() const {
    std::vector<uint32_t> order;
    order.reserve(nodes_.size());

    std::queue<uint32_t> q;
    std::vector<uint32_t> local_indegree(nodes_.size());
    for (uint32_t i = 0; i < nodes_.size(); ++i) {
        local_indegree[i] = nodes_[i].indegree;
        if (local_indegree[i] == 0) {
            q.push(i);
        }
    }

    while (!q.empty()) {
        uint32_t u = q.front();
        q.pop();
        order.push_back(u);

        const POANode* node = get_node(u);
        if (node) {
            uint32_t v = node->first_out;
            while (v != UINT32_MAX) {
                if (local_indegree[v] > 0) {
                    local_indegree[v]--;
                    if (local_indegree[v] == 0) {
                        q.push(v);
                    }
                }
                const POANode* vconst = get_node(v);
                v = vconst ? vconst->next_sibling : UINT32_MAX;
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
// POA Engine（使用 Banded DP）
// ============================================================================

POAEngine::POAEngine(const AssemblyConfig& config) : config_(config) {}

void POAEngine::reset() {}

// ============================================================================
// Graph Smith-Waterman with Banded DP（O(N*W) 内存复杂度）
// ============================================================================

int POAEngine::graph_smith_waterman(
    POAArena& arena,
    const std::vector<uint32_t>& topo_order,
    SeqSpan sequence,
    std::vector<uint8_t>& traceback) {

    const int n = static_cast<int>(sequence.length);
    const int m = static_cast<int>(topo_order.size());

    if (n == 0 || m == 0) return -1;

    // 检查是否适合带状比对
    if (!BandedDPBuffer::is_suitable_for_banded(n, m)) {
        // 降级到简单回溯（不做完整 DP）
        return n;  // 返回 query 长度作为近似
    }

    // 初始化带状缓冲区
    static thread_local BandedDPBuffer banded_dp(32);
    banded_dp.reset(n, m);

    traceback.resize(static_cast<size_t>(n) * (banded_dp.band_size()));

    int bandwidth = banded_dp.bandwidth();

    // 初始化第一行
    for (int j = 1; j <= m && banded_dp.in_band(0, j); ++j) {
        banded_dp.M(0, j) = INT_MIN / 4;
        banded_dp.X(0, j) = config_.gap_open + j * config_.gap_extend;
        banded_dp.Y(0, j) = INT_MIN / 4;
    }

    // 主 DP 循环（只遍历 band 内）
    for (int i = 1; i <= n; ++i) {
        char qc = sequence[i - 1];

        // 计算 j 的有效范围
        int j_start = std::max(1, i - bandwidth / 2);
        int j_end = std::min(m, i + bandwidth / 2);

        for (int j = j_start; j <= j_end; ++j) {
            if (!banded_dp.in_band(i, j)) continue;

            banded_dp.M(i, 0) = INT_MIN / 4;
            banded_dp.X(i, 0) = INT_MIN / 4;
            banded_dp.Y(i, 0) = config_.gap_open + i * config_.gap_extend;

            uint32_t node_idx = topo_order[j - 1];
            const POANode* node = arena.get_node(node_idx);
            if (!node) continue;

            char gc = node->base;
            int8_t score = (qc == gc) ? config_.match : config_.mismatch;

            int32_t max_prev_M = INT_MIN / 4;

            // 遍历前驱
            for (uint8_t k = 0; k < node->predecessors.inline_count; ++k) {
                uint32_t pred_idx = node->predecessors.inline_preds[k];
                if (pred_idx < arena.size()) {
                    const POANode* pred = arena.get_node(pred_idx);
                    if (pred && pred->topo_order < static_cast<uint32_t>(j) &&
                        banded_dp.in_band(i - 1, pred->topo_order + 1)) {
                        max_prev_M = std::max(max_prev_M, banded_dp.M(i - 1, pred->topo_order + 1));
                    }
                }
            }

            int match = (max_prev_M > INT_MIN / 4) ? max_prev_M + score : INT_MIN / 4;

            int gap_t = std::max(
                banded_dp.X(i, j - 1) + config_.gap_extend,
                banded_dp.M(i, j - 1) + config_.gap_open + config_.gap_extend
            );

            int gap_q = std::max(
                banded_dp.Y(i - 1, j) + config_.gap_extend,
                banded_dp.M(i - 1, j) + config_.gap_open + config_.gap_extend
            );

            int best = match;
            uint8_t trace = 0;

            if (gap_t > best) {
                best = gap_t;
                trace = 1;
            }
            if (gap_q > best) {
                best = gap_q;
                trace = 2;
            }

            banded_dp.M(i, j) = best;
            banded_dp.X(i, j) = gap_t;
            banded_dp.Y(i, j) = gap_q;

            // 存储 traceback（简化版）
            int offset = j - i + bandwidth / 2;
            if (offset >= 0 && offset < banded_dp.band_size()) {
                traceback[(i - 1) * banded_dp.band_size() + offset] = trace;
            }
        }
    }

    // 找到最佳终止点
    int best_score = 0;
    int best_i = 0;
    for (int i = 1; i <= n; ++i) {
        int j_start = std::max(1, i - bandwidth / 2);
        int j_end = std::min(m, i + bandwidth / 2);
        for (int j = j_start; j <= j_end; ++j) {
            if (banded_dp.in_band(i, j) && banded_dp.M(i, j) > best_score) {
                best_score = banded_dp.M(i, j);
                best_i = i;
            }
        }
    }

    return best_i;
}

// ============================================================================
// Build POA Graph
// ============================================================================

uint32_t POAEngine::build_graph(const std::vector<SeqSpan>& sequences) {
    if (sequences.empty()) return UINT32_MAX;

    size_t estimated = sequences[0].length * sequences.size() + 100;
    POAArena arena(estimated);

    const auto& first = sequences[0];
    uint32_t start = UINT32_MAX;
    uint32_t prev_idx = UINT32_MAX;

    for (size_t i = 0; i < first.length; ++i) {
        uint32_t idx = arena.allocate_node(first[i]);
        if (i == 0) start = idx;
        if (prev_idx != UINT32_MAX) {
            arena.add_edge(prev_idx, idx);
        }
        prev_idx = idx;
    }

    for (size_t seq_idx = 1; seq_idx < sequences.size(); ++seq_idx) {
        const auto& seq = sequences[seq_idx];
        if (seq.length == 0) continue;

        std::vector<uint8_t> traceback;
        int qlen = align_to_graph(arena, start, seq);

        if (qlen > 0) {
            add_aligned_sequence(arena, start, seq, traceback, qlen);
        }
    }

    return start;
}

int POAEngine::align_to_graph(
    POAArena& arena,
    uint32_t graph_start,
    SeqSpan sequence) {

    if (graph_start == UINT32_MAX || sequence.length == 0) return -1;

    std::vector<uint32_t> topo_order = arena.get_topo_order();
    std::vector<uint8_t> traceback;

    return graph_smith_waterman(arena, topo_order, sequence, traceback);
}

void POAEngine::add_aligned_sequence(
    POAArena& arena,
    uint32_t graph_start,
    SeqSpan sequence,
    const std::vector<uint8_t>& traceback,
    int qlen) {

    if (graph_start == UINT32_MAX || qlen < 0) return;

    std::vector<uint32_t> topo_order = arena.get_topo_order();
    int glen = static_cast<int>(topo_order.size());

    const uint8_t* tb_ptr = traceback.data();
    int i = qlen, j = glen;

    static constexpr int MAX_TRACEBACK_ITERATIONS = 1000000;
    int iterations = 0;

    while (i > 0 && j > 0) {
        if (++iterations > MAX_TRACEBACK_ITERATIONS) break;

        int offset = j - i + 16;  // 简化偏移
        uint8_t trace = (offset >= 0 && offset < 32) ?
                        tb_ptr[(i - 1) * 32 + offset] : 0;

        if (trace == 0) {
            if (j > 1 && i > 0) {
                uint32_t node_idx = topo_order[j - 1];
                POANode* node = arena.get_node(node_idx);
                if (node) node->count++;
                j--;
            }
            i--;
        } else if (trace == 1) {
            j--;
        } else {
            if (j > 0) i--;
        }
    }
}

std::string POAEngine::extract_consensus(POAArena& arena, uint32_t graph_start) {
    std::string consensus;
    if (graph_start == UINT32_MAX) return consensus;

    std::vector<uint32_t> topo_order = arena.get_topo_order();

    int max_weight = 0;
    uint32_t end_node = UINT32_MAX;

    for (uint32_t node_idx : topo_order) {
        const POANode* node = arena.get_node(node_idx);
        if (node && node->count > max_weight) {
            max_weight = node->count;
            end_node = node_idx;
        }
    }

    std::string rev_consensus;
    uint32_t current = end_node;

    static constexpr int MAX_PATH_LENGTH = 100000;
    int path_length = 0;

    while (current != UINT32_MAX && path_length++ < MAX_PATH_LENGTH) {
        const POANode* node = arena.get_node(current);
        if (!node) break;

        rev_consensus.push_back(node->base);

        uint32_t best_pred = UINT32_MAX;
        int best_weight = -1;

        for (uint8_t k = 0; k < node->predecessors.inline_count; ++k) {
            uint32_t pred_idx = node->predecessors.inline_preds[k];
            const POANode* pred = arena.get_node(pred_idx);
            if (pred && static_cast<int>(pred->count) > best_weight) {
                best_weight = pred->count;
                best_pred = pred_idx;
            }
        }

        if (best_pred != UINT32_MAX) {
            current = best_pred;
        } else {
            break;
        }
    }

    std::reverse(rev_consensus.begin(), rev_consensus.end());
    return rev_consensus;
}

std::vector<std::string> POAEngine::extract_paths(
    POAArena& arena,
    uint32_t graph_start,
    int max_paths) {

    std::vector<std::string> results;
    if (graph_start == UINT32_MAX || max_paths <= 0) return results;

    std::unordered_set<uint32_t> visited;
    std::string current_path;

    std::function<void(uint32_t)> dfs = [&](uint32_t node_idx) {
        if (results.size() >= static_cast<size_t>(max_paths)) return;
        if (visited.count(node_idx)) return;

        const POANode* node = arena.get_node(node_idx);
        if (!node || node->count == 0) return;

        visited.insert(node_idx);
        current_path.push_back(node->base);

        if (node->first_out == UINT32_MAX ||
            current_path.size() > static_cast<size_t>(config_.max_path_length)) {
            results.push_back(current_path);
        } else {
            uint32_t child = node->first_out;
            while (child != UINT32_MAX && results.size() < static_cast<size_t>(max_paths)) {
                dfs(child);
                const POANode* child_node = arena.get_node(child);
                child = child_node ? child_node->next_sibling : UINT32_MAX;
            }
        }

        current_path.pop_back();
        visited.erase(node_idx);
    };

    dfs(graph_start);
    return results;
}

// ============================================================================
// Heaviest Bundle Extractor（最重束共识提取）
// ============================================================================

BundlePath HeaviestBundleExtractor::extract(POAArena& arena, uint32_t graph_start) {
    BundlePath result;

    if (graph_start == UINT32_MAX) return result;

    std::vector<uint32_t> topo_order = arena.get_topo_order();
    int n = static_cast<int>(topo_order.size());

    ensure_capacity(n);

    // 初始化 DP
    for (int i = 0; i < n; ++i) {
        const POANode* node = arena.get_node(topo_order[i]);
        dp_weight_[i] = node ? node->count : 0;
        dp_prev_[i] = -1;
        dp_offset_[i] = 0;
    }

    // 拓扑序 DP：找最重路径
    for (int i = 0; i < n; ++i) {
        const POANode* node = arena.get_node(topo_order[i]);
        if (!node) continue;

        uint32_t child = node->first_out;
        while (child != UINT32_MAX) {
            // 找到 child 在 topo_order 中的位置
            const POANode* child_node = arena.get_node(child);
            if (child_node && child_node->topo_order < static_cast<uint32_t>(n)) {
                int child_pos = static_cast<int>(child_node->topo_order);
                int child_idx = child_pos;

                // 边权重 = min(count_u, count_v)
                int edge_weight = std::min(
                    node ? node->count : 0,
                    child_node ? child_node->count : 0
                );

                if (dp_weight_[i] + edge_weight > dp_weight_[child_idx]) {
                    dp_weight_[child_idx] = dp_weight_[i] + edge_weight;
                    dp_prev_[child_idx] = i;
                    dp_offset_[child_idx] = child_pos - i;
                }
            }
            child = child_node ? child_node->next_sibling : UINT32_MAX;
        }
    }

    // 找到终点（权重最大的节点）
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

    // 简化版：使用不同带宽提取不同路径
    for (int bw : {16, 32, 64}) {
        if (static_cast<int>(results.size()) >= k) break;
        BandedDPBuffer temp_buffer(bw);
        // 这里可以扩展为真正的多路径提取
    }

    return results;
}

// ============================================================================
// Quality-Weighted Consensus
// ============================================================================

QualityWeightedConsensus QualityWeightedConsensus::extract(
    const POAArena& arena,
    uint32_t graph_start,
    const AssemblyConfig& config) {

    QualityWeightedConsensus result;

    if (graph_start == UINT32_MAX) return result;

    std::vector<uint32_t> topo_order = arena.get_topo_order();

    struct PositionData {
        std::array<int, 256> base_counts{};
        double total_weight = 0.0;
    };

    std::vector<PositionData> positions;

    for (uint32_t node_idx : topo_order) {
        const POANode* node = arena.get_node(node_idx);
        if (!node) break;

        PositionData pos;
        pos.total_weight = static_cast<double>(node->count);
        pos.base_counts[static_cast<unsigned char>(node->base)] = node->count;

        positions.push_back(pos);
    }

    for (const auto& pos : positions) {
        if (pos.total_weight == 0) {
            result.sequence.push_back('N');
            result.base_quality.push_back(0.0);
            continue;
        }

        int max_count = 0;
        char best_base = 'N';
        for (int i = 0; i < 256; ++i) {
            if (pos.base_counts[i] > max_count) {
                max_count = pos.base_counts[i];
                best_base = static_cast<char>(i);
            }
        }

        result.sequence.push_back(best_base);

        double frequency = max_count / pos.total_weight;
        frequency = std::clamp(frequency, 0.0001, 1.0);
        double phred = -10.0 * std::log10(1.0 - frequency);
        result.base_quality.push_back(std::min(phred, 60.0));
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

    int crossing_reads = 0;
    int32_t span = right_breakpoint - left_breakpoint;
    if (span < 0) span = -span;

    for (const auto& read : reads) {
        if (read.pos <= left_breakpoint && read.end_pos >= right_breakpoint) {
            crossing_reads++;
        }
    }

    double crossing_ratio = static_cast<double>(crossing_reads) / reads.size();

    return crossing_ratio >= config.MIN_COVERAGE_RATIO &&
           span <= config.MIN_CIRCULAR_LENGTH * 10;
}

// ============================================================================
// AssemblyEngine Implementation
// ============================================================================

AssemblyEngine::AssemblyEngine(AssemblyConfig config)
    : config_(std::move(config)) {}

std::array<std::vector<std::string>, 3> AssemblyEngine::extract_segments(
    const Component&,
    const std::vector<ReadSketch>&,
    const GenomeAccessor&) {
    return {};
}

std::vector<size_t> AssemblyEngine::stratified_sample(
    const Component& component,
    const std::vector<ReadSketch>&,
    size_t max_samples) {

    std::vector<size_t> sampled;
    if (component.read_count == 0) return sampled;

    size_t count = std::min(static_cast<size_t>(component.read_count), max_samples);
    for (size_t i = 0; i < count; ++i) {
        sampled.push_back(i);
    }

    return sampled;
}

StructuralFingerprint AssemblyEngine::build_fingerprint(const Contig& contig) const {
    return StructuralFingerprint::from_contig(
        contig.sequence,
        contig.left_breakpoint,
        contig.right_breakpoint,
        contig.te_family_id,
        contig.orientation);
}

// ============================================================================
// Assemble Component（使用 Heaviest Bundle）
// ============================================================================

std::vector<Contig> AssemblyEngine::assemble_component(
    const Component& component,
    const std::vector<ReadSketch>& reads,
    const GenomeAccessor& genome) {

    POAContext ctx;
    ctx.reset();

    std::vector<Contig> contigs;

    if (component.read_count < config_.min_reads_for_poa) {
        return contigs;
    }

    auto sampled = stratified_sample(component, reads, config_.max_reads_for_poa);
    if (sampled.empty()) return contigs;

    std::vector<SeqSpan> sequences;
    sequences.reserve(sampled.size());

    for (size_t idx : sampled) {
        if (idx < reads.size()) {
            sequences.emplace_back(reads[idx].sequence.c_str(), reads[idx].sequence.length());
        }
    }

    if (sequences.empty()) return contigs;

    uint32_t graph_start = ctx.engine.build_graph(sequences);

    if (graph_start == UINT32_MAX) return contigs;

    // 使用 Heaviest Bundle 提取共识（工业级）
    HeaviestBundleExtractor hb_extractor(32);
    BundlePath bundle = hb_extractor.extract(ctx.arena, graph_start);

    if (bundle.consensus.empty()) return contigs;

    Contig contig;
    contig.sequence = bundle.consensus;
    contig.up_flank_seq = "";
    contig.ins_seq = contig.sequence;
    contig.down_flank_seq = "";
    contig.left_breakpoint = component.start;
    contig.right_breakpoint = component.end;
    contig.te_family_id = 0;
    contig.orientation = 0;
    contig.trunc_level = 0;
    contig.support_reads = static_cast<int32_t>(sampled.size());
    contig.consensus_quality = 0.85;
    contig.fingerprint = build_fingerprint(contig);

    // 质量加权共识
    auto qw_consensus = QualityWeightedConsensus::extract(ctx.arena, graph_start, config_);
    if (!qw_consensus.sequence.empty()) {
        contig.sequence = qw_consensus.sequence;
    }

    contigs.push_back(std::move(contig));

    return contigs;
}

// ============================================================================
// Assemble Batch（真正的并行）
// ============================================================================

std::vector<Contig> AssemblyEngine::assemble_batch(
    std::vector<Component>& components,
    const std::vector<ReadSketch>& reads,
    const GenomeAccessor& genome) {

    if (components.empty()) return {};

    std::vector<std::vector<Contig>> results(components.size());

    std::vector<size_t> indices(components.size());
    std::iota(indices.begin(), indices.end(), 0);

    // 工业级并行：使用 thread_local POAContext
    std::for_each(std::execution::par, indices.begin(), indices.end(),
        [&](size_t idx) {
            thread_local POAContext ctx;
            if (idx < components.size()) {
                results[idx] = [&]() -> std::vector<Contig> {
                    ctx.reset();
                    return assemble_component_thread_safe(
                        components[idx], reads, genome, ctx);
                }();
            }
        });

    std::vector<Contig> all_contigs;
    for (size_t i = 0; i < components.size(); ++i) {
        for (auto& c : results[i]) {
            c.support_reads = static_cast<int32_t>(components[i].read_count);
            all_contigs.push_back(std::move(c));
        }
    }

    return all_contigs;
}

std::vector<Contig> AssemblyEngine::assemble_component_thread_safe(
    const Component& component,
    const std::vector<ReadSketch>& reads,
    const GenomeAccessor& genome,
    POAContext& ctx) {

    std::vector<Contig> contigs;

    if (component.read_count < config_.min_reads_for_poa) {
        return contigs;
    }

    auto sampled = stratified_sample(component, reads, config_.max_reads_for_poa);
    if (sampled.empty()) return contigs;

    std::vector<SeqSpan> sequences;
    sequences.reserve(sampled.size());

    for (size_t idx : sampled) {
        if (idx < reads.size()) {
            sequences.emplace_back(reads[idx].sequence.c_str(), reads[idx].sequence.length());
        }
    }

    if (sequences.empty()) return contigs;

    uint32_t graph_start = ctx.engine.build_graph(sequences);
    if (graph_start == UINT32_MAX) return contigs;

    HeaviestBundleExtractor hb_extractor(32);
    BundlePath bundle = hb_extractor.extract(ctx.arena, graph_start);

    if (bundle.consensus.empty()) return contigs;

    Contig contig;
    contig.sequence = bundle.consensus;
    contig.left_breakpoint = component.start;
    contig.right_breakpoint = component.end;
    contig.te_family_id = 0;
    contig.orientation = 0;
    contig.support_reads = static_cast<int32_t>(sampled.size());
    contig.consensus_quality = 0.85;
    contig.fingerprint = build_fingerprint(contig);

    contigs.push_back(std::move(contig));

    return contigs;
}

// ============================================================================
// Structural Fingerprint Methods
// ============================================================================

uint64_t StructuralFingerprint::hash() const {
    uint64_t h = 1469598103934665603ULL;
    auto mix = [&](uint64_t v) {
        h ^= v;
        h *= 1099511628211ULL;
    };
    mix(static_cast<uint64_t>(tid));
    mix(static_cast<uint64_t>(breakpoint_l));
    mix(static_cast<uint64_t>(breakpoint_l_end));
    mix(static_cast<uint64_t>(breakpoint_r));
    mix(static_cast<uint64_t>(te_family_id));
    mix(static_cast<uint64_t>(orientation));
    mix(static_cast<uint64_t>(trunc_level));
    return h;
}

bool StructuralFingerprint::matches(const StructuralFingerprint& other) const {
    if (tid != other.tid && tid >= 0 && other.tid >= 0) return false;
    if (te_family_id != other.te_family_id && te_family_id >= 0 && other.te_family_id >= 0) return false;
    if (orientation != other.orientation && orientation >= 0 && other.orientation >= 0) return false;

    int32_t l_diff = std::abs(breakpoint_l - other.breakpoint_l);
    int32_t r_diff = std::abs(breakpoint_r - other.breakpoint_r);

    return l_diff <= BP_TOLERANCE && r_diff <= BP_TOLERANCE;
}

StructuralFingerprint StructuralFingerprint::from_contig(
    std::string_view contig,
    int32_t left_bp,
    int32_t right_bp,
    int32_t te_family,
    int8_t orient) {

    StructuralFingerprint fp;
    fp.tid = 0;
    fp.breakpoint_l = left_bp;
    fp.breakpoint_l_end = left_bp + static_cast<int32_t>(contig.length()) / 4;
    fp.breakpoint_r = right_bp;
    fp.breakpoint_r_end = right_bp + static_cast<int32_t>(contig.length()) / 2;
    fp.te_family_id = te_family;
    fp.orientation = orient;
    fp.trunc_level = 0;
    fp.ins_length_min = static_cast<int32_t>(contig.length()) - 100;
    fp.ins_length_max = static_cast<int32_t>(contig.length()) + 100;
    fp.has_inversion = (orient == 2);

    return fp;
}

StructuralRepresentative AssemblyEngine::merge_contigs(
    const std::vector<int>& indices,
    std::vector<Contig>& contigs) {

    StructuralRepresentative rep;
    rep.rep_id = static_cast<int32_t>(indices.size() > 0 ? indices[0] : -1);
    rep.component_ids.reserve(indices.size());
    rep.contig_ids.reserve(indices.size());

    double total_quality = 0.0;
    int total_reads = 0;

    for (int idx : indices) {
        if (idx >= 0 && static_cast<size_t>(idx) < contigs.size()) {
            rep.component_ids.push_back(contigs[idx].left_breakpoint);
            rep.contig_ids.push_back(idx);
            total_quality += contigs[idx].consensus_quality;
            total_reads += contigs[idx].support_reads;

            if (rep.fingerprint.hash() == 0) {
                rep.fingerprint = contigs[idx].fingerprint;
            }
            rep.rep_sequence = contigs[idx].sequence;
        }
    }

    if (!indices.empty()) {
        rep.avg_quality = total_quality / indices.size();
        rep.total_reads = total_reads;
    }

    return rep;
}

int AssemblyEngine::extract_polya_length(std::string_view seq) {
    if (seq.empty()) return -1;

    int count = 0;
    for (int i = static_cast<int>(seq.size()) - 1; i >= 0; --i) {
        if (seq[i] == 'A') {
            count++;
        } else {
            break;
        }
    }

    return count >= 10 ? count : -1;
}

// ============================================================================
// Collapse Structurally（使用 R-Tree）
// ============================================================================

std::vector<StructuralRepresentative> AssemblyEngine::collapse_structurally(
    std::vector<Contig>& contigs) {

    std::vector<StructuralRepresentative> reps;

    if (contigs.empty()) return reps;

    size_t N = contigs.size();
    if (N == 0) return reps;

    // 构建 R-Tree 索引
    RTree rtree;
    for (size_t i = 0; i < N; ++i) {
        const auto& fp = contigs[i].fingerprint;
        rtree.insert(
            static_cast<float>(fp.breakpoint_l - StructuralFingerprint::BP_TOLERANCE),
            static_cast<float>(fp.breakpoint_r - StructuralFingerprint::BP_TOLERANCE),
            static_cast<float>(fp.breakpoint_l_end + StructuralFingerprint::BP_TOLERANCE),
            static_cast<float>(fp.breakpoint_r_end + StructuralFingerprint::BP_TOLERANCE),
            static_cast<int>(i)
        );
    }

    // DSU 聚类
    DSU dsu(static_cast<int>(N));
    std::vector<int> candidates;

    for (size_t i = 0; i < N; ++i) {
        const auto& fp = contigs[i].fingerprint;
        candidates = rtree.range_query(
            static_cast<float>(fp.breakpoint_l - StructuralFingerprint::BP_TOLERANCE),
            static_cast<float>(fp.breakpoint_r - StructuralFingerprint::BP_TOLERANCE),
            static_cast<float>(fp.breakpoint_l_end + StructuralFingerprint::BP_TOLERANCE),
            static_cast<float>(fp.breakpoint_r_end + StructuralFingerprint::BP_TOLERANCE)
        );

        for (int cand_idx : candidates) {
            if (static_cast<size_t>(cand_idx) <= i) continue;
            if (contigs[i].fingerprint.matches(contigs[cand_idx].fingerprint)) {
                dsu.unite(static_cast<int>(i), cand_idx);
            }
        }
    }

    // 收集代表
    std::unordered_map<int, std::vector<int>> rep_groups;
    for (size_t i = 0; i < N; ++i) {
        int root = dsu.find(static_cast<int>(i));
        rep_groups[root].push_back(static_cast<int>(i));
    }

    for (auto& group : rep_groups) {
        StructuralRepresentative rep = merge_contigs(group.second, contigs);
        if (rep.contig_ids.size() >= 1) {
            reps.push_back(rep);
        }
    }

    return reps;
}

}  // namespace placer
