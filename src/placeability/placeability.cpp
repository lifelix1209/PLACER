#include "placeability.h"
#include "local_realign.h"
#include <algorithm>
#include <numeric>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <array>
#include <cassert>

namespace placer {

// ============================================================================
// TSD Detector Implementation
// ============================================================================

TSDDetector::TSDDetector(const PlaceabilityConfig& config)
    : config_(config) {
    for (int i = 0; i < 16; ++i) {
        dinuc_bg_freq_[i] = 1.0 / 16.0;
    }
}

int TSDDetector::find_lcp(std::string_view a, std::string_view b) const {
    int len = std::min(static_cast<int>(a.size()),
                       static_cast<int>(b.size()));
    int i = 0;
    while (i < len && toupper_fast(a[i]) == toupper_fast(b[i])) {
        ++i;
    }
    return i;
}

int TSDDetector::find_lcs(std::string_view a, std::string_view b) const {
    int len = std::min(static_cast<int>(a.size()),
                       static_cast<int>(b.size()));
    int i = 0;
    while (i < len &&
           toupper_fast(a[a.size() - 1 - i]) ==
           toupper_fast(b[b.size() - 1 - i])) {
        ++i;
    }
    return i;
}

// [5.5 修正] fuzzy 搜索长度上限不再漂移，capped 到 max_tsd_length
int TSDDetector::find_lcp_fuzzy(std::string_view a, std::string_view b,
                                 int max_mismatches,
                                 int& out_mismatches) const {
    int len = std::min(static_cast<int>(a.size()),
                       static_cast<int>(b.size()));
    // [5.5 修正] 搜索范围仍然 capped 到 max_tsd_length
    // 错配只影响是否 accept，不扩展搜索长度
    len = std::min(len, config_.max_tsd_length);

    int mismatches = 0;
    int best_len = 0;
    int best_mm = 0;

    for (int i = 0; i < len; ++i) {
        if (toupper_fast(a[i]) != toupper_fast(b[i])) {
            mismatches++;
            if (mismatches > max_mismatches) break;
        }
        int current_len = i + 1;
        double mm_ratio = static_cast<double>(mismatches) / current_len;
        if (mm_ratio <= config_.max_tsd_mismatch_ratio &&
            current_len >= config_.min_tsd_length) {
            best_len = current_len;
            best_mm = mismatches;
        }
    }

    // [5.5 修正] 最终结果也 clamp 到 max_tsd_length
    best_len =    std::min(best_len, config_.max_tsd_length);

    out_mismatches = best_mm;
    return best_len;
}

int TSDDetector::find_lcs_fuzzy(std::string_view a, std::string_view b,
                                 int max_mismatches,
                                 int& out_mismatches) const {
    int len = std::min(static_cast<int>(a.size()),
                       static_cast<int>(b.size()));
    // [5.5 修正] 同上
    len = std::min(len, config_.max_tsd_length);

    int mismatches = 0;
    int best_len = 0;
    int best_mm = 0;

    for (int i = 0; i < len; ++i) {
        char ca = toupper_fast(a[a.size() - 1 - i]);
        char cb = toupper_fast(b[b.size() - 1 - i]);
        if (ca != cb) {
            mismatches++;
            if (mismatches > max_mismatches) break;
        }
        int current_len = i + 1;
        double mm_ratio = static_cast<double>(mismatches) / current_len;
        if (mm_ratio <= config_.max_tsd_mismatch_ratio &&
            current_len >= config_.min_tsd_length) {
            best_len = current_len;
            best_mm = mismatches;
        }
    }

    // [5.5 修正] clamp
    best_len = std::min(best_len, config_.max_tsd_length);

    out_mismatches = best_mm;
    return best_len;
}

double TSDDetector::calculate_background_freq(std::string_view seq) {
    if (seq.size() < 2) {
        return 1.0;
    }

    double log_prob = 0.0;

    for (size_t i = 0; i + 1 < seq.size(); ++i) {
        int idx = dinuc_index(toupper_fast(seq[i]),
                              toupper_fast(seq[i + 1]));
        if (idx >= 0 && idx < 16) {
            double freq = dinuc_bg_freq_[idx];
            log_prob += std::log(std::max(freq, 1e-10));
        } else {
            log_prob += std::log(0.25);
        }
    }

    return std::exp(log_prob);
}

int TSDDetector::dinuc_index(char a, char b) const {
    auto base_idx = [](char c) -> int {
        switch (c) {
            case 'A': return 0;
            case 'C': return 1;
            case 'G': return 2;
            case 'T': return 3;
            default: return -1;
        }
    };
    int ai = base_idx(a);
    int bi = base_idx(b);
    if (ai < 0 || bi < 0) return -1;
    return ai * 4 + bi;
}

char TSDDetector::toupper_fast(char c) const {
    return (c >= 'a' && c <= 'z') ? (c - 32) : c;
}

void TSDDetector::set_background(std::string_view bg_sequence) {
    background_ = std::string(bg_sequence);

    if (bg_sequence.size() < 100) {
        for (int i = 0; i < 16; ++i) {
            dinuc_bg_freq_[i] = 1.0 / 16.0;
        }
        return;
    }

    std::array<int, 16> counts{};
    int total = 0;

    for (size_t i = 0; i + 1 < bg_sequence.size(); ++i) {
        int idx = dinuc_index(toupper_fast(bg_sequence[i]),
                              toupper_fast(bg_sequence[i + 1]));
        if (idx >= 0 && idx < 16) {
            counts[idx]++;
            total++;
        }
    }

    double pseudo = 1.0;
    double denom = total + 16.0 * pseudo;

    for (int i = 0; i < 16; ++i) {
        dinuc_bg_freq_[i] = (counts[i] + pseudo) / denom;
    }
}

// ============================================================================
// TSD Detection Main Logic
// [5.4 修正] sliding window 对称更新 left_bp 和 right_bp
// ============================================================================

TSDResult TSDDetector::detect(
    std::string_view left_flank,
    std::string_view right_flank,
    int32_t left_bp,
    int32_t right_bp) {

    TSDResult result;
    result.found = false;
    result.tsd_length = 0;
    result.mismatch_count = 0;
    result.mismatch_ratio = 0.0;
    result.background_freq = 1.0;
    result.is_significant = false;
    result.left_bp = left_bp;
    result.right_bp = right_bp;
    result.bp_offset = 0;

    if (left_flank.empty() || right_flank.empty()) {
        return result;
    }

    TSDResult best_result = result;
    double best_significance_score = 0.0;

    // === 方法 1: 精确前缀匹配 ===
    {
        int lcp = find_lcp(left_flank, right_flank);
        lcp = std::min(lcp, config_.max_tsd_length);

        if (lcp >= config_.min_tsd_length) {
            TSDResult r;
            r.found = true;
            r.tsd_length = lcp;
            r.tsd_seq = std::string(left_flank.substr(0, lcp));
            r.mismatch_count = 0;
            r.mismatch_ratio = 0.0;
            r.left_bp = left_bp;
            r.right_bp = right_bp;
            r.bp_offset = 0;
            r.background_freq = calculate_background_freq(r.tsd_seq);
            r.is_significant = is_significant(r, config_);
            r.detection_method = "exact_prefix";

            double sig_score = score_tsd_candidate(r);
            if (sig_score > best_significance_score) {
                best_significance_score = sig_score;
                best_result = r;
            }
        }
    }

    // === 方法 2: 精确后缀匹配 ===
    {
        int lcs = find_lcs(left_flank, right_flank);
        lcs = std::min(lcs, config_.max_tsd_length);

        if (lcs >= config_.min_tsd_length) {
            TSDResult r;
            r.found = true;
            r.tsd_length = lcs;
            r.tsd_seq = std::string(
                left_flank.substr(left_flank.size() - lcs));
            r.mismatch_count = 0;
            r.mismatch_ratio = 0.0;
            r.left_bp = left_bp;
            r.right_bp = right_bp;
            r.bp_offset = 0;
            r.background_freq = calculate_background_freq(r.tsd_seq);
            r.is_significant = is_significant(r, config_);
            r.detection_method = "exact_suffix";

            double sig_score = score_tsd_candidate(r);
            if (sig_score > best_significance_score) {
                best_significance_score = sig_score;
                best_result = r;
            }
        }
    }

    // === 方法 3: 模糊前缀匹配 ===
    if (!best_result.found ||
        best_result.tsd_length < config_.min_tsd_length + 2) {
        int mm = 0;
        int fuzzy_lcp = find_lcp_fuzzy(
            left_flank, right_flank,
            config_.max_tsd_mismatches, mm);

        if (fuzzy_lcp >= config_.min_tsd_length &&
            fuzzy_lcp > best_result.tsd_length) {
            TSDResult r;
            r.found = true;
            r.tsd_length = fuzzy_lcp;
            r.tsd_seq = std::string(left_flank.substr(0, fuzzy_lcp));
            r.mismatch_count = mm;
            r.mismatch_ratio = static_cast<double>(mm) / fuzzy_lcp;
            r.left_bp = left_bp;
            r.right_bp = right_bp;
            r.bp_offset = 0;
            r.background_freq = calculate_background_freq(r.tsd_seq);
            r.is_significant = is_significant(r, config_);
            r.detection_method = "fuzzy_prefix";

            double sig_score = score_tsd_candidate(r);
            if (sig_score > best_significance_score) {
                best_significance_score = sig_score;
                best_result = r;
            }
        }
    }

    // === 方法 4: 模糊后缀匹配 ===
    if (!best_result.found ||
        best_result.tsd_length < config_.min_tsd_length + 2) {
        int mm = 0;
        int fuzzy_lcs = find_lcs_fuzzy(
            left_flank, right_flank,
            config_.max_tsd_mismatches, mm);

        if (fuzzy_lcs >= config_.min_tsd_length &&
            fuzzy_lcs > best_result.tsd_length) {
            TSDResult r;
            r.found = true;
            r.tsd_length = fuzzy_lcs;
            r.tsd_seq = std::string(
                left_flank.substr(left_flank.size() - fuzzy_lcs));
            r.mismatch_count = mm;
            r.mismatch_ratio = static_cast<double>(mm) / fuzzy_lcs;
            r.left_bp = left_bp;
            r.right_bp = right_bp;
            r.bp_offset = 0;
            r.background_freq = calculate_background_freq(r.tsd_seq);
            r.is_significant = is_significant(r, config_);
            r.detection_method = "fuzzy_suffix";

            double sig_score = score_tsd_candidate(r);
            if (sig_score > best_significance_score) {
                best_significance_score = sig_score;
                best_result = r;
            }
        }
    }

    // === 方法 5: 滑动窗口搜索 ===
    // [5.4 修正] offset 定义明确为"断点偏移量"，对称更新两侧断点
    if (!best_result.found &&
        left_flank.size() >= 10 && right_flank.size() >= 10) {

        int search_range = std::min(
            config_.tsd_search_window,
            static_cast<int>(
                std::min(left_flank.size(), right_flank.size())) / 2);

        for (int offset = -search_range; offset <= search_range; ++offset) {
            if (offset == 0) continue;  // 已在方法 1-4 中处理

            std::string_view lf_shifted, rf_shifted;

            // offset > 0: 左侧 flank 向右偏移（左断点右移）
            // offset < 0: 右侧 flank 向右偏移（右断点左移）
            if (offset > 0) {
                if (static_cast<size_t>(offset) >= left_flank.size()) continue;
                lf_shifted = left_flank.substr(offset);
                rf_shifted = right_flank;
            } else {
                if (static_cast<size_t>(-offset) >= right_flank.size()) continue;
                lf_shifted = left_flank;
                rf_shifted = right_flank.substr(-offset);
            }

            int mm = 0;
            int shifted_lcp = find_lcp_fuzzy(
                lf_shifted, rf_shifted, 1, mm);

            if (shifted_lcp >= config_.min_tsd_length + 1) {
                TSDResult r;
                r.found = true;
                r.tsd_length = shifted_lcp;
                r.tsd_seq = std::string(lf_shifted.substr(0, shifted_lcp));
                r.mismatch_count = mm;
                r.mismatch_ratio = static_cast<double>(mm) / shifted_lcp;
                r.bp_offset = offset;

                // [5.4 修正] 对称更新断点位置
                // offset > 0: 左断点右移 offset bp
                // offset < 0: 右断点左移 |offset| bp
                if (offset > 0) {
                    r.left_bp = left_bp + offset;
                    r.right_bp = right_bp;
                } else {
                    r.left_bp = left_bp;
                    r.right_bp = right_bp + offset;  // offset 是负数，右断点左移
                }

                r.background_freq = calculate_background_freq(r.tsd_seq);
                r.is_significant = is_significant(r, config_);
                r.detection_method = "sliding_window";

                double sig_score = score_tsd_candidate(r);
                if (sig_score > best_significance_score) {
                    best_significance_score = sig_score;
                    best_result = r;
                }
            }
        }
    }

    return best_result;
}

double TSDDetector::score_tsd_candidate(const TSDResult& tsd) const {
    if (!tsd.found) return 0.0;

    double len_score = std::min(
        static_cast<double>(tsd.tsd_length) / config_.max_tsd_length, 1.0);

    double mm_penalty = 1.0 - tsd.mismatch_ratio * 2.0;
    mm_penalty = std::max(mm_penalty, 0.0);

    double bg_score = 0.0;
    if (tsd.background_freq > 0 && tsd.background_freq < 1.0) {
        bg_score = -std::log10(tsd.background_freq);
        bg_score = std::min(bg_score / 10.0, 1.0);
    }

    double complexity = calculate_sequence_complexity(tsd.tsd_seq);

    // [新增] 偏移惩罚：偏移越大越不可信
    double offset_penalty = 1.0;
    if (tsd.bp_offset != 0) {
        offset_penalty = 1.0 / (1.0 + std::abs(tsd.bp_offset) * 0.1);
    }

    return (len_score * 0.25 + mm_penalty * 0.2 +
            bg_score * 0.3 + complexity * 0.15) * offset_penalty;
}

double TSDDetector::calculate_sequence_complexity(
    std::string_view seq) const {

    if (seq.size() < 2) return 0.0;

    std::array<int, 16> counts{};
    int total = 0;

    for (size_t i = 0; i + 1 < seq.size(); ++i) {
        int idx = dinuc_index(toupper_fast(seq[i]),
                              toupper_fast(seq[i + 1]));
        if (idx >= 0) {
            counts[idx]++;
            total++;
        }
    }

    if (total == 0) return 0.0;

    double entropy = 0.0;
    for (int i = 0; i < 16; ++i) {
        if (counts[i] > 0) {
            double p = static_cast<double>(counts[i]) / total;
            entropy -= p * std::log2(p);
        }
    }

    return entropy / 4.0;
}

bool TSDDetector::is_significant(
    const TSDResult& tsd,
    const PlaceabilityConfig& config) {

    if (!tsd.found || tsd.tsd_length < config.min_tsd_length) {
        return false;
    }

    if (tsd.mismatch_ratio > config.max_tsd_mismatch_ratio) {
        return false;
    }

    const std::string& seq = tsd.tsd_seq;

    // 单核苷酸重复过滤
    std::array<int, 4> base_counts{};
    for (char c : seq) {
        switch (toupper_fast(c)) {
            case 'A': base_counts[0]++; break;
            case 'C': base_counts[1]++; break;
            case 'G': base_counts[2]++; break;
            case 'T': base_counts[3]++; break;
        }
    }

    int max_single = *std::max_element(
        base_counts.begin(), base_counts.end());
    double max_single_ratio =
        static_cast<double>(max_single) / seq.size();

    if (max_single_ratio > 0.80) {
        return false;
    }

    // 二核苷酸重复过滤
    if (seq.size() >= 6) {
        bool is_dinuc_repeat = true;
        char d1 = toupper_fast(seq[0]);
        char d2 = toupper_fast(seq[1]);

        for (size_t i = 0; i < seq.size(); ++i) {
            char expected = (i % 2 == 0) ? d1 : d2;
            if (toupper_fast(seq[i]) != expected) {
                is_dinuc_repeat = false;
                break;
            }
        }

        if (is_dinuc_repeat && d1 != d2) {
            if (tsd.tsd_length < 12) {
                return false;
            }
        }
    }

    // 序列复杂度
    double complexity = calculate_sequence_complexity(seq);
    if (complexity < config.min_tsd_complexity) {
        return false;
    }

    // 背景频率 + 长度动态阈值
    double length_adjusted_threshold = config.tsd_bg_threshold;
    if (tsd.tsd_length <= 6) {
        length_adjusted_threshold *= 0.1;
    } else if (tsd.tsd_length <= 10) {
        length_adjusted_threshold *= 0.5;
    }

    if (tsd.background_freq >= length_adjusted_threshold) {
        return false;
    }

    return true;
}

// ============================================================================
// PlaceabilityScorer Implementation
// ============================================================================

PlaceabilityScorer::PlaceabilityScorer(const PlaceabilityConfig& config)
    : config_(config) {}

double PlaceabilityScorer::calculate_delta(
    double best_score, double second_best_score) {

    double second = std::max(second_best_score, 0.0);
    return std::max(best_score - second, 0.0);
}

bool PlaceabilityScorer::check_side_consistency(
    int32_t left_best, int32_t left_second,
    int32_t right_best, int32_t right_second,
    int gap_threshold) {

    if (left_best < 0 || right_best < 0) {
        return true;
    }

    int32_t distance = std::abs(left_best - right_best);
    if (distance > gap_threshold) {
        return false;
    }

    if (left_second >= 0 && right_second >= 0) {
        int32_t second_distance = std::abs(left_second - right_second);
        if (second_distance <= gap_threshold &&
            std::abs(left_second - left_best) > gap_threshold) {
            return false;
        }
    }

    return true;
}

double PlaceabilityScorer::calculate_support_consistency(
    const std::vector<double>& scores) {

    if (scores.empty()) return 0.0;
    if (scores.size() == 1) return 1.0;

    double sum = std::accumulate(scores.begin(), scores.end(), 0.0);
    double mean = sum / scores.size();

    if (mean <= 0) return 0.0;

    double sq_sum = 0.0;
    for (double s : scores) {
        sq_sum += (s - mean) * (s - mean);
    }

    double variance = sq_sum / scores.size();
    double std_dev = std::sqrt(variance);
    double cv = std_dev / mean;

    double consistency = 1.0 / (1.0 + cv * cv);
    return consistency;
}

std::vector<int32_t> PlaceabilityScorer::extract_candidate_loci(
    const std::vector<LocusEvidence>& evidence) {

    std::vector<int32_t> loci;
    loci.reserve(evidence.size());

    for (const auto& e : evidence) {
        if (e.locus_pos >= 0) {
            loci.push_back(e.locus_pos);
        }
    }

    std::sort(loci.begin(), loci.end());
    loci.erase(std::unique(loci.begin(), loci.end()), loci.end());

    return loci;
}

std::vector<double> PlaceabilityScorer::calculate_locus_scores(
    const std::vector<LocusEvidence>& evidence,
    int32_t& best_pos,
    int32_t& second_best_pos) {

    best_pos = -1;
    second_best_pos = -1;

    if (evidence.empty()) return {};

    struct LocusStats {
        double score_sum = 0.0;
        double max_score = 0.0;
        int count = 0;
        std::unordered_set<size_t> unique_reads;
        int forward = 0;
        int reverse = 0;
    };

    std::unordered_map<int32_t, LocusStats> locus_map;

    for (const auto& e : evidence) {
        if (e.locus_pos < 0) continue;

        auto& stats = locus_map[e.locus_pos];
        stats.score_sum += e.normalized_score;
        stats.max_score = std::max(stats.max_score, e.normalized_score);
        stats.count++;
        stats.unique_reads.insert(e.read_idx);

        if (e.is_reverse) stats.reverse++;
        else stats.forward++;
    }

    struct ScoredLocus {
        int32_t pos;
        double combined_score;
        int unique_read_count;
    };

    std::vector<ScoredLocus> scored;
    scored.reserve(locus_map.size());

    for (const auto& [pos, stats] : locus_map) {
        double avg = stats.score_sum / std::max(stats.count, 1);
        int unique_reads = static_cast<int>(stats.unique_reads.size());

        double strand_bonus = 1.0;
        if (stats.forward > 0 && stats.reverse > 0) {
            int min_s = std::min(stats.forward, stats.reverse);
            int max_s = std::max(stats.forward, stats.reverse);
            strand_bonus = 1.0 + 0.3 *
                (static_cast<double>(min_s) / max_s);
        }

        double combined = avg *
            std::log2(unique_reads + 1) *
            strand_bonus;

        scored.push_back({pos, combined, unique_reads});
    }

    std::sort(scored.begin(), scored.end(),
        [](const ScoredLocus& a, const ScoredLocus& b) {
            return a.combined_score > b.combined_score;
        });

    if (!scored.empty()) {
        best_pos = scored[0].pos;
    }
    if (scored.size() > 1) {
        second_best_pos = scored[1].pos;
    }

    std::vector<double> result;
    result.reserve(scored.size());
    for (const auto& s : scored) {
        result.push_back(s.combined_score);
    }

    return result;
}

// ============================================================================
// [5.2 修正] side_consistency: 使用 LocusEvidence 中的 side 字段
// 而非依赖未定义的 evidence_bits 0x10/0x20
// ============================================================================

bool PlaceabilityScorer::calculate_side_consistency(
    const std::vector<LocusEvidence>& evidence) {

    if (evidence.empty()) return true;

    // [5.2 修正] 不再依赖 evidence_bits 的 0x10/0x20
    // 而是使用 LocusEvidence 中的 read_pos 与 locus_pos 的关系来判断侧向
    //
    // 策略：
    //   read_pos < locus_pos → 该 read 提供左侧约束
    //   read_pos >= locus_pos → 该 read 提供右侧约束
    //   如果 LocusEvidence 没有 read_pos 字段，
    //   则使用 is_reverse 作为近似（但标记为低置信）

    std::unordered_map<int32_t, std::vector<std::pair<double, bool>>>
        locus_sides;

    bool has_position_info = false;

    for (const auto& e : evidence) {
        if (e.locus_pos < 0) continue;

        bool is_left = false;
        bool is_right = false;

        // 优先使用 read_pos（如果可用）
        if (e.read_pos >= 0) {
            has_position_info = true;
            is_left = (e.read_pos < e.locus_pos);
            is_right = (e.read_pos >= e.locus_pos);
        }
        // 其次使用 read_end_pos（如果可用）
        else if (e.read_end_pos > 0) {
            has_position_info = true;
            int32_t read_mid = (e.read_pos + e.read_end_pos) / 2;
            is_left = (read_mid < e.locus_pos);
            is_right = (read_mid >= e.locus_pos);
        }
        // 最后回退到 is_reverse 近似
        else {
            is_left = e.is_reverse;
            is_right = !e.is_reverse;
        }

        if (is_left) {
            locus_sides[e.locus_pos].push_back(
                {e.normalized_score, true});
        }
        if (is_right) {
            locus_sides[e.locus_pos].push_back(
                {e.normalized_score, false});
        }
    }

    // [5.2 修正] 如果完全没有位置信息，不让 side_consistency 影响 Tier
    // 返回 true 但标记为"未验证"
    if (!has_position_info) {
        // 无法可靠判断，返回 true 但调用方应知道这是默认值
        return true;
    }

    int consistent_loci = 0;
    int total_loci = 0;

    for (const auto& [pos, sides] : locus_sides) {
        bool has_left = false;
        bool has_right = false;
        double left_score_sum = 0.0;
        double right_score_sum = 0.0;
        int left_count = 0;
        int right_count = 0;

        for (const auto& [score, is_left] : sides) {
            if (is_left) {
                has_left = true;
                left_score_sum += score;
                left_count++;
            } else {
                has_right = true;
                right_score_sum += score;
                right_count++;
            }
        }

        total_loci++;

        if (has_left && has_right) {
            double left_avg = left_count > 0 ?
                left_score_sum / left_count : 0;
            double right_avg = right_count > 0 ?
                right_score_sum / right_count : 0;

            double ratio = std::min(left_avg, right_avg) /
                           std::max(left_avg + 1e-10, right_avg + 1e-10);

            if (ratio >= 0.3) {
                consistent_loci++;
            }
        }
    }

    if (total_loci == 0) return true;

    return consistent_loci >= (total_loci + 1) / 2;
}

// ============================================================================
// [5.1 修正] aggregate_evidence: old_bits 用 uint32_t 而非 double
// ============================================================================

std::vector<LocusEvidence> PlaceabilityScorer::aggregate_evidence(
    const std::vector<LocusEvidence>& evidence) {

    if (evidence.empty()) return {};

    struct EvidenceKey {
        size_t read_idx;
        int32_t locus_pos;

        bool operator==(const EvidenceKey& other) const {
            return read_idx == other.read_idx &&
                   locus_pos == other.locus_pos;
        }
    };

    struct EvidenceKeyHash {
        size_t operator()(const EvidenceKey& k) const {
            return std::hash<size_t>()(k.read_idx) ^
                   (std::hash<int32_t>()(k.locus_pos) << 16);
        }
    };

    std::unordered_map<EvidenceKey, LocusEvidence, EvidenceKeyHash> merged;

    for (const auto& e : evidence) {
        EvidenceKey key{e.read_idx, e.locus_pos};
        auto it = merged.find(key);

        if (it == merged.end()) {
            merged[key] = e;
        } else {
            auto& existing = it->second;

            // [5.1 修正] 使用 uint32_t 而非 double 存储 bitmask
            uint32_t old_bits = existing.evidence_bits;

            if (e.normalized_score > existing.normalized_score) {
                existing = e;
                existing.evidence_bits |= old_bits;
            } else {
                existing.evidence_bits |= e.evidence_bits;
            }
            existing.total_score = std::max(
                existing.total_score, e.total_score);
        }
    }

    std::vector<LocusEvidence> result;
    result.reserve(merged.size());
    for (auto& [key, ev] : merged) {
        result.push_back(std::move(ev));
    }

    return result;
}

// ============================================================================
// Tier Determination
// [5.2 修正] 当 side_consistency 基于 is_reverse 近似时，
// 不给 Tier1 加成，避免虚高
// ============================================================================

Tier PlaceabilityScorer::determine_tier(
    const ExtendedPlaceabilityReport& report,
    const PlaceabilityConfig& config) {

    // Tier 1: 高置信唯一定位
    if (report.delta_score >= config.delta_score_threshold &&
        report.side_consistent &&
        report.side_consistency_verified &&  // [5.2 修正] 必须是经过验证的
        report.candidate_count <= config.max_locus_for_tier1 &&
        report.support_consistency >= config.min_support_consistency &&
        report.unique_support_reads >= config.min_locus_support &&
        !report.is_ambiguous) {
        return Tier::TIER1;
    }

    // Tier 1 降级版：side_consistency 未验证但其他指标极好
    if (report.delta_score >= config.delta_score_threshold * 1.5 &&
        report.candidate_count <= config.max_locus_for_tier1 &&
        report.support_consistency >= config.min_support_consistency * 1.2 &&
        report.unique_support_reads >= config.min_locus_support * 2 &&
        report.strand_balanced &&
        !report.is_ambiguous) {
        return Tier::TIER1;
    }

    // Tier 2
    if (report.delta_score >= config.delta_score_tier2 &&
        (report.side_consistent || report.support_consistency >= 0.6) &&
        report.candidate_count <= config.max_candidate_locus &&
        report.unique_support_reads >= config.min_locus_support) {

        if (report.strand_balanced || report.unique_support_reads >= 5) {
            return Tier::TIER2;
        }
    }

    // Tier 3
    if (report.delta_score > 0 &&
        report.unique_support_reads >= 2) {
        return Tier::TIER3;
    }

    // Tier 4
    if (report.unique_support_reads >= 1) {
        return Tier::TIER4;
    }

    // Tier 5
    return Tier::TIER5;
}

// ============================================================================
// [5.3 修正] strand_ratio 计算
// [5.6 修正] best_locus 初值处理
// ============================================================================

ExtendedPlaceabilityReport PlaceabilityScorer::extend_placeability_report(
    const PlaceabilityReport& existing_report,
    const std::vector<LocusEvidence>& raw_evidence) {

    ExtendedPlaceabilityReport report;

    // 聚合证据
    auto evidence = aggregate_evidence(raw_evidence);

    // [5.6 修正] 记录是否有外部提供的 best_locus
    // 如果有，后续新计算会覆盖它（明确语义：重新评估）
    bool had_external_best = (existing_report.best_locus_pos >= 0);
    report.best_locus = existing_report.best_locus_pos;
    report.delta_score = existing_report.delta_score;

    // unique reads 统计
    std::unordered_set<size_t> unique_reads;
    for (const auto& e : evidence) {
        unique_reads.insert(e.read_idx);
    }
    report.support_reads = static_cast<int>(evidence.size());
    report.unique_support_reads = static_cast<int>(unique_reads.size());

    // 候选位点
    report.all_loci = extract_candidate_loci(evidence);
    report.candidate_count = static_cast<int>(report.all_loci.size());

    if (report.candidate_count > config_.max_candidate_locus) {
        report.tier = Tier::TIER3;
        report.overall_score = 0.0;
        report.is_ambiguous = true;
        return report;
    }

    // 位点得分
    int32_t best_pos = -1, second_best_pos = -1;
    report.locus_scores = calculate_locus_scores(
        evidence, best_pos, second_best_pos);

    // [5.6 修正] 新计算的 best_pos 覆盖外部值（重新评估语义）
    // 如果新计算有效，使用新值；否则保留外部值
    if (best_pos >= 0) {
        report.best_locus = best_pos;
    }
    report.second_best_locus = second_best_pos;

    // Delta Score
    if (report.locus_scores.size() >= 2) {
        report.delta_score = calculate_delta(
            report.locus_scores[0],
            report.locus_scores[1]);
    } else if (report.locus_scores.size() == 1) {
        report.delta_score = report.locus_scores[0];
    }

    // 支持度一致性
    if (!evidence.empty()) {
        std::vector<double> scores;
        scores.reserve(evidence.size());
        for (const auto& e : evidence) {
            if (e.normalized_score > 0) {
                scores.push_back(e.normalized_score);
            }
        }
        report.support_consistency =
            calculate_support_consistency(scores);
    }

    // 侧翼一致性
    // [5.2 修正] 检查是否有 read_pos 信息来判断 side_consistency 是否可靠
    bool has_position_info = false;
    for (const auto& e : evidence) {
        if (e.read_pos >= 0 || e.read_end_pos > 0) {
            has_position_info = true;
            break;
        }
    }

    report.side_consistent = calculate_side_consistency(evidence);
    report.side_consistency_verified = has_position_info;

    // [5.3 修正] strand_ratio = min(fwd,rev) / max(fwd,rev)
    // 而非 min(fwd,rev) / (fwd+rev)
    int fwd = 0, rev = 0;
    for (const auto& e : evidence) {
        if (e.normalized_score < 0.1) continue;
        if (e.is_reverse) rev++;
        else fwd++;
    }
    report.forward_count = fwd;
    report.reverse_count = rev;
    report.strand_balanced = (fwd > 0 && rev > 0);

    // [5.3 修正] 正确的平衡比：min/max，而非 min/(fwd+rev)
    if (fwd > 0 || rev > 0) {
        report.strand_ratio = static_cast<float>(
            std::min(fwd, rev)) /
            static_cast<float>(std::max(std::max(fwd, rev), 1));
    } else {
        report.strand_ratio = 0.0f;
    }

    // 歧义检测
    if (report.locus_scores.size() >= 2) {
        double ratio = report.locus_scores[1] /
            std::max(report.locus_scores[0], 1e-6);
        report.is_ambiguous = (ratio > 0.8);
    }

    // 各候选位点的支持数
    std::unordered_map<int32_t, int> locus_count;
    std::unordered_map<int32_t, std::unordered_set<size_t>> locus_reads;
    for (const auto& e : evidence) {
        if (e.locus_pos >= 0) {
            locus_count[e.locus_pos]++;
            locus_reads[e.locus_pos].insert(e.read_idx);
        }
    }

    report.locus_support.clear();
    report.locus_unique_reads.clear();
    for (int32_t locus : report.all_loci) {
        report.locus_support.push_back(locus_count[locus]);
        report.locus_unique_reads.push_back(
            static_cast<int>(locus_reads[locus].size()));
    }

    // Tier
    report.tier = determine_tier(report, config_);

    // 综合评分
    report.overall_score = calculate_overall_score(report);

    return report;
}

ExtendedPlaceabilityReport PlaceabilityScorer::calculate_full_placeability(
    const std::vector<LocusEvidence>& raw_evidence) {

    PlaceabilityReport base_report;
    base_report.best_locus_pos = -1;
    base_report.delta_score = 0.0;

    return extend_placeability_report(base_report, raw_evidence);
}

// [5.2/5.3 修正] overall_score: side_consistent 加成仅在 verified 时生效
double PlaceabilityScorer::calculate_overall_score(
    const ExtendedPlaceabilityReport& report) const {

    double base = report.delta_score * report.support_consistency;

    // [5.2 修正] 侧翼一致加成仅在经过验证时生效
    if (report.side_consistent && report.side_consistency_verified) {
        base *= 1.3;
    }
    // 未验证的 side_consistent 不加成也不惩罚

    // [5.3 修正] 链平衡加成使用修正后的 strand_ratio
    if (report.strand_balanced) {
        base *= (1.0 + 0.2 * report.strand_ratio);
    }

    // 支持数量加成
    if (report.unique_support_reads > 0) {
        base *= std::log2(report.unique_support_reads + 1) / 3.0;
    }

    // 歧义惩罚
    if (report.is_ambiguous) {
        base *= 0.5;
    }

    // 候选数量惩罚
    if (report.candidate_count > 1) {
        base *= 1.0 / std::log2(report.candidate_count + 1);
    }

    return std::max(base, 0.0);
}

// ============================================================================
// PlaceabilityOutput Implementation
// ============================================================================

PlaceabilityOutput::PlaceabilityOutput(const PlaceabilityConfig& config)
    : config_(config) {}

std::string PlaceabilityOutput::generate_info_fields(
    const ExtendedPlaceabilityReport& report) {

    std::ostringstream oss;

    oss << "PLACEABILITY=" << std::fixed << std::setprecision(3)
        << report.overall_score;

    oss << ";TIER=" << static_cast<int>(report.tier);

    oss << ";DELTA=" << std::fixed << std::setprecision(2)
        << report.delta_score;

    oss << ";CONSISTENCY=" << std::fixed << std::setprecision(3)
        << report.support_consistency;

    oss << ";CANDIDATES=" << report.candidate_count;

    oss << ";SUPPORT=" << report.unique_support_reads;

    if (report.side_consistent) {
        if (report.side_consistency_verified) {
            oss << ";SIDES=CONSISTENT_VERIFIED";
        } else {
            oss << ";SIDES=CONSISTENT_UNVERIFIED";
        }
    } else {
        oss << ";SIDES=INCONSISTENT";
    }

    if (report.strand_balanced) {
        oss << ";STRAND=BALANCED";
        // [5.3 修正] 输出修正后的 strand_ratio
        oss << ";STRAND_RATIO=" << std::fixed << std::setprecision(2)
            << report.strand_ratio;
    } else {
        oss << ";STRAND=UNBALANCED";
    }

    if (report.is_ambiguous) {
        oss << ";AMBIGUOUS=1";
    }

    if (report.best_locus >= 0) {
        oss << ";BEST_LOCUS=" << report.best_locus;
    }
    if (report.second_best_locus >= 0) {
        oss << ";SECOND_LOCUS=" << report.second_best_locus;
    }

    if (!report.locus_support.empty()) {
        oss << ";LOCUS_SUPPORT=";
        for (size_t i = 0; i < report.locus_support.size(); ++i) {
            if (i > 0) oss << ",";
            oss << report.all_loci[i] << ":"
                << report.locus_support[i];
            // [新增] 也输出 unique reads
            if (i < report.locus_unique_reads.size()) {
                oss << "/" << report.locus_unique_reads[i];
            }
        }
    }

    return oss.str();
}

std::string PlaceabilityOutput::get_tier_description(Tier tier) {
    switch (tier) {
        case Tier::TIER1:
            return "Unique placement with high confidence, "
                   "verified consistent flanks and strand balance";
        case Tier::TIER2:
            return "Acceptable placement with moderate confidence, "
                   "consistent structure";
        case Tier::TIER3:
            return "Low confidence placement, possible ambiguity "
                   "or inconsistency";
        case Tier::TIER4:
            return "Very low confidence, minimal supporting evidence";
        case Tier::TIER5:
            return "Unable to place, no reliable evidence";
        default:
            return "Unknown tier";
    }
}

std::string PlaceabilityOutput::generate_tsd_info(const TSDResult& tsd) {
    std::ostringstream oss;

    if (!tsd.found) {
        oss << "TSD=NOT_FOUND";
        return oss.str();
    }

    oss << "TSD_FOUND=1";
    oss << ";TSD_LEN=" << tsd.tsd_length;

    if (!tsd.tsd_seq.empty()) {
        oss << ";TSD_SEQ=" << tsd.tsd_seq;
    }

    if (tsd.mismatch_count > 0) {
        oss << ";TSD_MISMATCH=" << tsd.mismatch_count;
        oss << ";TSD_MM_RATIO=" << std::fixed << std::setprecision(3)
            << tsd.mismatch_ratio;
    }

    oss << ";TSD_SIG=" << (tsd.is_significant ? "1" : "0");

    if (!tsd.detection_method.empty()) {
        oss << ";TSD_METHOD=" << tsd.detection_method;
    }

    if (tsd.bp_offset != 0) {
        oss << ";TSD_BP_OFFSET=" << tsd.bp_offset;
    }

    oss << ";TSD_BG_FREQ=" << std::scientific << std::setprecision(2)
        << tsd.background_freq;

    return oss.str();
}

std::string PlaceabilityOutput::generate_full_info(
    const ExtendedPlaceabilityReport& report,
    const TSDResult& tsd) {

    std::string info = generate_info_fields(report);
    info += ";";
    info += generate_tsd_info(tsd);
    return info;
}

std::string PlaceabilityOutput::generate_json(
    const ExtendedPlaceabilityReport& report,
    const TSDResult& tsd) {

    std::ostringstream oss;
    oss << "{";

    oss << "\"placeability\":{";
    oss << "\"overall_score\":" << std::fixed << std::setprecision(4)
        << report.overall_score;
    oss << ",\"tier\":" << static_cast<int>(report.tier);
    oss << ",\"tier_desc\":\"" << get_tier_description(report.tier)
        << "\"";
    oss << ",\"delta_score\":" << std::setprecision(3)
        << report.delta_score;
    oss << ",\"consistency\":" << report.support_consistency;
    oss << ",\"candidates\":" << report.candidate_count;
    oss << ",\"unique_reads\":" << report.unique_support_reads;
    oss << ",\"side_consistent\":"
        << (report.side_consistent ? "true" : "false");
    oss << ",\"side_verified\":"
        << (report.side_consistency_verified ? "true" : "false");
    oss << ",\"strand_balanced\":"
        << (report.strand_balanced ? "true" : "false");
    oss << ",\"strand_ratio\":" << std::setprecision(3)
        << report.strand_ratio;
    oss << ",\"ambiguous\":"
        << (report.is_ambiguous ? "true" : "false");

    if (report.best_locus >= 0) {
        oss << ",\"best_locus\":" << report.best_locus;
    }
    if (report.second_best_locus >= 0) {
        oss << ",\"second_locus\":" << report.second_best_locus;
    }

    if (!report.all_loci.empty()) {
        oss << ",\"loci\":[";
        for (size_t i = 0; i < report.all_loci.size(); ++i) {
            if (i > 0) oss << ",";
            oss << "{\"pos\":" << report.all_loci[i];
            if (i < report.locus_scores.size()) {
                oss << ",\"score\":" << std::setprecision(3)
                    << report.locus_scores[i];
            }
            if (i < report.locus_support.size()) {
                oss << ",\"evidence_count\":"
                    << report.locus_support[i];
            }
            if (i < report.locus_unique_reads.size()) {
                oss << ",\"unique_reads\":"
                    << report.locus_unique_reads[i];
            }
            oss << "}";
        }
        oss << "]";
    }

    oss << "}";

    oss << ",\"tsd\":{";
    oss << "\"found\":" << (tsd.found ? "true" : "false");
    if (tsd.found) {
        oss << ",\"length\":" << tsd.tsd_length;
        oss << ",\"seq\":\"" << tsd.tsd_seq << "\"";
        oss << ",\"mismatches\":" << tsd.mismatch_count;
        oss << ",\"mismatch_ratio\":" << std::fixed
            << std::setprecision(3) << tsd.mismatch_ratio;
        oss << ",\"significant\":"
            << (tsd.is_significant ? "true" : "false");
        if (!tsd.detection_method.empty()) {
            oss << ",\"method\":\"" << tsd.detection_method << "\"";
        }
        if (tsd.bp_offset != 0) {
            oss << ",\"bp_offset\":" << tsd.bp_offset;
        }
        oss << ",\"left_bp\":" << tsd.left_bp;
        oss << ",\"right_bp\":" << tsd.right_bp;
        oss << ",\"bg_freq\":" << std::scientific
            << std::setprecision(3) << tsd.background_freq;
    }
    oss << "}";

    oss << "}";
    return oss.str();
}

}  // namespace placer


