#ifndef PLACER_WINDOW_STATS_H
#define PLACER_WINDOW_STATS_H

#include <string>
#include <unordered_map>
#include <cmath>

namespace placer {

/**
 * 简单的滑动窗口分位数追踪器
 * 不需要存储所有观察值，内存 O(1)
 */
class QuantileTracker {
public:
    explicit QuantileTracker(double quantile = 0.99)
        : quantile_(quantile), count_(0), sum_(0.0), min_(0.0), max_(0.0) {}

    void add(double value) {
        if (count_ == 0) {
            min_ = max_ = value;
        } else {
            min_ = std::min(min_, value);
            max_ = std::max(max_, value);
        }
        sum_ += value;
        ++count_;
    }

    double estimate() const {
        if (count_ == 0) return 0.0;
        // 简单估计：假设值近似均匀分布
        return min_ + quantile_ * (max_ - min_);
    }

    size_t count() const { return count_; }
    double sum() const { return sum_; }
    double min() const { return min_; }
    double max() const { return max_; }
    double mean() const { return count_ > 0 ? sum_ / count_ : 0.0; }

    void reset() {
        count_ = 0;
        sum_ = 0.0;
        min_ = max_ = 0.0;
    }

private:
    double quantile_;
    size_t count_;
    double sum_;
    double min_;
    double max_;
};

/**
 * 窗口统计收集器
 * 按染色体收集分位数估计
 */
class WindowStatsCollector {
public:
    void add_clip_bp(const std::string& chrom, int64_t value) {
        clip_[chrom].add(static_cast<double>(value));
    }
    void add_sa_reads(const std::string& chrom, int64_t value) {
        sa_[chrom].add(static_cast<double>(value));
    }
    void add_ins_events(const std::string& chrom, int64_t value) {
        ins_[chrom].add(static_cast<double>(value));
    }

    double get_q99_clip(const std::string& chrom) const {
        auto it = clip_.find(chrom);
        return (it != clip_.end() && it->second.count() > 0) ? it->second.estimate() : 1000.0;
    }
    double get_q99_sa(const std::string& chrom) const {
        auto it = sa_.find(chrom);
        return (it != sa_.end() && it->second.count() > 0) ? it->second.estimate() : 10.0;
    }
    double get_q99_ins(const std::string& chrom) const {
        auto it = ins_.find(chrom);
        return (it != ins_.end() && it->second.count() > 0) ? it->second.estimate() : 5.0;
    }

    void reset_chrom(const std::string& chrom) {
        clip_.erase(chrom);
        sa_.erase(chrom);
        ins_.erase(chrom);
    }

    void reset_all() {
        clip_.clear();
        sa_.clear();
        ins_.clear();
    }

private:
    std::unordered_map<std::string, QuantileTracker> clip_;
    std::unordered_map<std::string, QuantileTracker> sa_;
    std::unordered_map<std::string, QuantileTracker> ins_;
};

}  // namespace placer

#endif  // PLACER_WINDOW_STATS_H
