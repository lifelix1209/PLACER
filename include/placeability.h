#ifndef PLACER_PLACEABILITY_H
#define PLACER_PLACEABILITY_H

#include "local_realign.h"
#include "assembly.h"
#include <vector>
#include <string>
#include <cstdint>
#include <string_view>
#include <optional>
#include <array>
#include <limits>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <unordered_map>
#include <memory>
#include <sstream>
#include <iomanip>

namespace placer {

// ============================================================================
// Placeability Configuration
// ============================================================================

struct PlaceabilityConfig {
    // Delta Score 阈值
    double delta_score_threshold = 30.0;     // Tier 1 需要的高分差
    double delta_score_tier2 = 10.0;         // Tier 2 的最低分差

    // 侧翼一致性阈值
    int side_consistency_gap = 50;            // 侧翼落点差异容忍 (bp)

    // 候选位点约束
    int min_locus_support = 2;                // 最少支持 reads
    int max_locus_for_tier1 = 5;              // Tier 1 的最大候选数
    int max_candidate_locus = 20;             // 整体最大候选数（超过降级）

    // 支持度一致性阈值
    double min_support_consistency = 0.5;      // 最低支持度一致性

    // TSD 检测参数
    int min_tsd_length = 3;                   // 最小 TSD 长度
    int max_tsd_length = 50;                 // 最大 TSD 长度
    double tsd_mismatch_ratio = 0.1;          // TSD 匹配允许的错配率
    double tsd_bg_threshold = 0.05;           // TSD 显著性背景频率阈值

    // 额外参数
    int max_tsd_mismatches = 2;               // TSD 模糊匹配最大错配
    int tsd_search_window = 200;               // TSD 搜索窗口
    double min_tsd_complexity = 0.5;          // TSD 最小复杂度
    double max_tsd_mismatch_ratio = 0.2;       // TSD 最大错配比例
};

// ============================================================================
// Tier Enumeration
// ============================================================================

enum class Tier : int8_t {
    TIER1 = 1,
    TIER2 = 2,
    TIER3 = 3,
    TIER4 = 4,
    TIER5 = 5,
    UNTYPED = -1
};

// ============================================================================
// Extended Placeability Report
// ============================================================================

struct ExtendedPlaceabilityReport {
    // 核心指标
    int32_t best_locus = -1;              // 最佳落点
    int32_t second_best_locus = -1;       // 第二佳落点
    double delta_score = 0.0;              // 最佳-第二佳分差

    // 侧翼一致性
    bool side_consistent = false;          // 侧翼一致性
    bool side_consistency_verified = false; // 侧向一致性是否经过 read_pos 验证
    double support_consistency = 0.0;       // 支持度一致性 [0, 1]

    // 候选统计
    int candidate_count = 0;               // 候选位点数
    int support_reads = 0;                 // 支持 reads 数

    // 综合评分
    double overall_score = 0.0;            // 综合评分
    Tier tier = Tier::UNTYPED;             // 分类结果

    // 详细信息
    std::vector<int32_t> all_loci;         // 所有候选位点
    std::vector<double> locus_scores;      // 各候选位点得分
    std::vector<int> locus_support;        // 各候选位点支持数

    // 链锁一致性
    bool chain_consistent = false;         // Read 链一致性
    double chain_ratio = 0.0;              // 链锁比例

    // 额外字段
    int unique_support_reads = 0;
    bool is_ambiguous = false;
    int forward_count = 0;
    int reverse_count = 0;
    bool strand_balanced = false;
    float strand_ratio = 0.0f;
    std::vector<int32_t> locus_unique_reads;
};

// ============================================================================
// TSD Result
// ============================================================================

struct TSDResult {
    bool found = false;
    int tsd_length = 0;
    std::string tsd_seq;
    int mismatch_count = 0;
    double mismatch_ratio = 0.0;

    // 显著性评估
    double background_freq = 0.0;          // 在背景中的频率
    bool is_significant = false;           // 是否显著（非随机重复）

    // 断点位置
    int32_t left_bp = 0;
    int32_t right_bp = 0;
    int32_t bp_offset = 0;
    std::string detection_method;
};

// ============================================================================
// TSD Detector
// ============================================================================

class TSDDetector {
public:
    explicit TSDDetector(const PlaceabilityConfig& config = PlaceabilityConfig());

    // 检测 TSD
    TSDResult detect(
        std::string_view left_flank,
        std::string_view right_flank,
        int32_t left_bp,
        int32_t right_bp);

    // 评估 TSD 显著性
    bool is_significant(
        const TSDResult& tsd,
        const PlaceabilityConfig& config);

    // 设置背景序列（用于频率计算）
    void set_background(std::string_view bg_sequence);

private:
    PlaceabilityConfig config_;

    // 二核苷酸背景频率
    double dinuc_bg_freq_[16];
    std::string background_;

    // 找最长公共前缀
    int find_lcp(std::string_view a, std::string_view b) const;

    // 找最长公共后缀
    int find_lcs(std::string_view a, std::string_view b) const;

    // 带容错的最长公共前缀
    int find_lcp_fuzzy(std::string_view a, std::string_view b,
                       int max_mismatches, int& out_mismatches) const;

    // 带容错的最长公共后缀
    int find_lcs_fuzzy(std::string_view a, std::string_view b,
                       int max_mismatches, int& out_mismatches) const;

    // 二核苷酸索引
    int dinuc_index(char a, char b) const;

    // 快速大写转换
    char toupper_fast(char c) const;

    // 计算序列在背景中的频率
    double calculate_background_freq(std::string_view seq);

    // 计算序列复杂度
    double calculate_sequence_complexity(std::string_view seq) const;

    // 评分 TSD 候选
    double score_tsd_candidate(const TSDResult& tsd) const;
};

// ============================================================================
// Placeability Scorer
// ============================================================================

class PlaceabilityScorer {
public:
    explicit PlaceabilityScorer(const PlaceabilityConfig& config = PlaceabilityConfig());

    // 完整 Placeability 评估
    ExtendedPlaceabilityReport calculate_full_placeability(
        const std::vector<LocusEvidence>& evidence);

    // 计算 Delta Score
    static double calculate_delta(double best_score, double second_best_score);

    // 检查侧翼一致性
    static bool check_side_consistency(
        int32_t left_best, int32_t left_second,
        int32_t right_best, int32_t right_second,
        int gap_threshold);

    // 计算支持度一致性
    static double calculate_support_consistency(const std::vector<double>& scores);

    // Tier 判定
    static Tier determine_tier(
        const ExtendedPlaceabilityReport& report,
        const PlaceabilityConfig& config);

    // 从 LocalRealigner 的结果计算扩展 Placeability
    ExtendedPlaceabilityReport extend_placeability_report(
        const PlaceabilityReport& existing_report,
        const std::vector<LocusEvidence>& evidence);

private:
    PlaceabilityConfig config_;

    // 收集每个位点的证据
    std::vector<LocusEvidence> aggregate_evidence(
        const std::vector<LocusEvidence>& evidence);

    // 计算各候选位点得分
    std::vector<double> calculate_locus_scores(
        const std::vector<LocusEvidence>& evidence,
        int32_t& best_idx,
        int32_t& second_best_idx);

    // 计算侧翼一致性
    bool calculate_side_consistency(const std::vector<LocusEvidence>& evidence);

    // 提取所有候选位点
    std::vector<int32_t> extract_candidate_loci(
        const std::vector<LocusEvidence>& evidence);

    // 计算综合评分
    double calculate_overall_score(const ExtendedPlaceabilityReport& report) const;
};

// ============================================================================
// Placeability 输出生成器
// ============================================================================

class PlaceabilityOutput {
public:
    explicit PlaceabilityOutput(const PlaceabilityConfig& config = PlaceabilityConfig());

    // 生成 INFO 字段字符串
    std::string generate_info_fields(const ExtendedPlaceabilityReport& report);

    // 生成 Tier 描述
    std::string get_tier_description(Tier tier);

    // 生成 TSD 相关 INFO
    std::string generate_tsd_info(const TSDResult& tsd);

    // 生成完整信息
    std::string generate_full_info(const ExtendedPlaceabilityReport& report,
                                   const TSDResult& tsd);

    // 生成 JSON
    std::string generate_json(const ExtendedPlaceabilityReport& report,
                              const TSDResult& tsd);

private:
    PlaceabilityConfig config_;
};

}  // namespace placer

#endif  // PLACER_PLACEABILITY_H
