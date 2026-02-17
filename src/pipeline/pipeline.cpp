#include "pipeline.h"

#include <algorithm>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <iostream>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

#include <htslib/sam.h>

namespace placer {
namespace {

constexpr int32_t kSoftClipSignalMin = 20;
constexpr int32_t kLongInsertionSignalMin = 50;
constexpr int32_t kLargeIndelEvidenceMin = 40;
constexpr double kSaHintWeight = 0.35;
constexpr int32_t kHistBinSize = 20;
constexpr int32_t kHistPadding = 120;
constexpr double kPeakMinWeight = 1.2;
constexpr int32_t kPeakMergeDistance = 120;
constexpr int32_t kWindowMergeGap = 80;
constexpr int32_t kWindowMaxExpandBins = 30;
constexpr double kWindowDropRatio = 0.35;
constexpr int32_t kMinWindowSpan = 120;
constexpr double kReadAssignMinScore = 0.8;
constexpr double kReadAmbiguousRatio = 0.75;
constexpr double kReadAssignDecayBp = 150.0;

enum class EvidenceKind : uint8_t {
    kSoftClip = 0,
    kIndel = 1,
    kSAHint = 2
};

struct EvidencePoint {
    size_t read_index = 0;
    int32_t pos = -1;
    double weight = 0.0;
    EvidenceKind kind = EvidenceKind::kSoftClip;
    int32_t signal_len = 0;
    uint8_t class_mask = 0;
};

struct ReadSignalSummary {
    std::string read_id;
    bool is_reverse = false;
    bool has_sa_or_supp = false;
    int32_t max_soft_clip = 0;
    int32_t max_ins = 0;
    uint8_t class_mask = 0;
};

struct EvidenceBundle {
    std::vector<EvidencePoint> points;
    std::vector<ReadSignalSummary> read_summaries;
};

struct DensityPeak {
    int32_t pos = -1;
    double weight = 0.0;
};

struct CandidateWindow {
    int32_t start = -1;
    int32_t end = -1;
    int32_t center = -1;
    double peak_weight = 0.0;
};

int32_t classify_tier(double delta_score) {
    if (delta_score >= 30.0) {
        return 1;
    }
    if (delta_score >= 10.0) {
        return 2;
    }
    return 3;
}

const char* theta_tag(InsertionTheta theta) {
    switch (theta) {
        case InsertionTheta::kFwd: return "FWD";
        case InsertionTheta::kRev: return "REV";
        default: return "NA";
    }
}

template <typename T>
class SafeQueue {
public:
    void push(T&& value) {
        {
            std::lock_guard<std::mutex> lock(mu_);
            queue_.push(std::move(value));
        }
        cv_.notify_one();
    }

    bool pop(T& out) {
        std::unique_lock<std::mutex> lock(mu_);
        cv_.wait(lock, [this]() { return finished_ || !queue_.empty(); });
        if (queue_.empty()) {
            return false;
        }
        out = std::move(queue_.front());
        queue_.pop();
        return true;
    }

    void close() {
        {
            std::lock_guard<std::mutex> lock(mu_);
            finished_ = true;
        }
        cv_.notify_all();
    }

private:
    std::queue<T> queue_;
    std::mutex mu_;
    std::condition_variable cv_;
    bool finished_ = false;
};

bool is_match_like(int op) {
    return op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF;
}

int find_first_non_hard_clip(const uint32_t* cigar, int32_t n_cigar) {
    for (int32_t i = 0; i < n_cigar; ++i) {
        if (bam_cigar_op(cigar[i]) != BAM_CHARD_CLIP) {
            return i;
        }
    }
    return -1;
}

int find_last_non_hard_clip(const uint32_t* cigar, int32_t n_cigar) {
    for (int32_t i = n_cigar - 1; i >= 0; --i) {
        if (bam_cigar_op(cigar[i]) != BAM_CHARD_CLIP) {
            return i;
        }
    }
    return -1;
}

int32_t weighted_median_position(std::vector<std::pair<int32_t, double>> pos_weights) {
    pos_weights.erase(
        std::remove_if(pos_weights.begin(), pos_weights.end(), [](const auto& row) {
            return row.second <= 0.0;
        }),
        pos_weights.end());
    if (pos_weights.empty()) {
        return -1;
    }

    std::sort(pos_weights.begin(), pos_weights.end(), [](const auto& a, const auto& b) {
        return a.first < b.first;
    });

    double total_weight = 0.0;
    for (const auto& row : pos_weights) {
        total_weight += row.second;
    }
    if (total_weight <= 0.0) {
        return pos_weights[pos_weights.size() / 2].first;
    }

    const double half_weight = 0.5 * total_weight;
    double acc_weight = 0.0;
    for (const auto& row : pos_weights) {
        acc_weight += row.second;
        if (acc_weight >= half_weight) {
            return row.first;
        }
    }
    return pos_weights.back().first;
}

std::vector<double> smooth_histogram(const std::vector<double>& hist) {
    std::vector<double> smooth(hist.size(), 0.0);
    for (size_t i = 0; i < hist.size(); ++i) {
        double value = hist[i];
        if (i >= 1) {
            value += 0.60 * hist[i - 1];
        }
        if (i + 1 < hist.size()) {
            value += 0.60 * hist[i + 1];
        }
        if (i >= 2) {
            value += 0.25 * hist[i - 2];
        }
        if (i + 2 < hist.size()) {
            value += 0.25 * hist[i + 2];
        }
        smooth[i] = value;
    }
    return smooth;
}

EvidenceBundle extract_evidence_points(
    const std::vector<BamRecordPtr>& bin_records,
    int32_t expected_tid) {
    EvidenceBundle bundle;
    bundle.read_summaries.resize(bin_records.size());

    for (size_t idx = 0; idx < bin_records.size(); ++idx) {
        if (!bin_records[idx]) {
            continue;
        }

        ReadView view(bin_records[idx].get());
        if (view.tid() != expected_tid) {
            continue;
        }

        ReadSignalSummary summary;
        summary.read_id = std::string(view.qname());
        summary.is_reverse = (view.flag() & BAM_FREVERSE) != 0;
        summary.has_sa_or_supp = view.has_sa_tag() || ((view.flag() & BAM_FSUPPLEMENTARY) != 0);

        const uint32_t* cigar = view.cigar();
        const int32_t n_cigar = view.n_cigar();
        if (!cigar || n_cigar <= 0) {
            if (summary.has_sa_or_supp) {
                summary.class_mask |= kCandidateSplitSaSupplementary;
            }
            bundle.read_summaries[idx] = std::move(summary);
            continue;
        }

        int32_t leading_soft = 0;
        int32_t trailing_soft = 0;
        const int first = find_first_non_hard_clip(cigar, n_cigar);
        const int last = find_last_non_hard_clip(cigar, n_cigar);
        if (first >= 0 && bam_cigar_op(cigar[first]) == BAM_CSOFT_CLIP) {
            leading_soft = static_cast<int32_t>(bam_cigar_oplen(cigar[first]));
        }
        if (last >= 0 && bam_cigar_op(cigar[last]) == BAM_CSOFT_CLIP) {
            trailing_soft = static_cast<int32_t>(bam_cigar_oplen(cigar[last]));
        }

        std::vector<EvidencePoint> local_points;
        local_points.reserve(8);

        int32_t ref_pos = view.pos();
        for (int32_t ci = 0; ci < n_cigar; ++ci) {
            const int op = bam_cigar_op(cigar[ci]);
            const int32_t len = static_cast<int32_t>(bam_cigar_oplen(cigar[ci]));

            if (op == BAM_CSOFT_CLIP) {
                summary.max_soft_clip = std::max(summary.max_soft_clip, len);
            } else if (op == BAM_CINS) {
                summary.max_ins = std::max(summary.max_ins, len);
            }

            if (op == BAM_CINS && len >= kLargeIndelEvidenceMin) {
                EvidencePoint point;
                point.read_index = idx;
                point.pos = ref_pos;
                point.weight = 0.80 + (static_cast<double>(std::min(len, 400)) / 200.0);
                point.kind = EvidenceKind::kIndel;
                point.signal_len = len;
                local_points.push_back(point);
            } else if (op == BAM_CDEL && len >= kLargeIndelEvidenceMin) {
                const double weight = 0.60 + (static_cast<double>(std::min(len, 400)) / 260.0);

                EvidencePoint left_point;
                left_point.read_index = idx;
                left_point.pos = ref_pos;
                left_point.weight = weight;
                left_point.kind = EvidenceKind::kIndel;
                left_point.signal_len = len;
                local_points.push_back(left_point);

                EvidencePoint right_point = left_point;
                right_point.pos = ref_pos + len;
                local_points.push_back(right_point);
            }

            if ((bam_cigar_type(op) & 2) != 0) {
                ref_pos += len;
            }
        }
        const int32_t ref_end = ref_pos;

        if (leading_soft >= kSoftClipSignalMin) {
            EvidencePoint point;
            point.read_index = idx;
            point.pos = view.pos();
            point.weight = 1.00 + (static_cast<double>(std::min(leading_soft, 500)) / 180.0);
            point.kind = EvidenceKind::kSoftClip;
            point.signal_len = leading_soft;
            local_points.push_back(point);
        }
        if (trailing_soft >= kSoftClipSignalMin) {
            EvidencePoint point;
            point.read_index = idx;
            point.pos = ref_end;
            point.weight = 1.00 + (static_cast<double>(std::min(trailing_soft, 500)) / 180.0);
            point.kind = EvidenceKind::kSoftClip;
            point.signal_len = trailing_soft;
            local_points.push_back(point);
        }

        if (summary.has_sa_or_supp) {
            EvidencePoint left_hint;
            left_hint.read_index = idx;
            left_hint.pos = view.pos();
            left_hint.weight = kSaHintWeight;
            left_hint.kind = EvidenceKind::kSAHint;
            local_points.push_back(left_hint);

            EvidencePoint right_hint = left_hint;
            right_hint.pos = ref_end;
            local_points.push_back(right_hint);
        }

        if (summary.has_sa_or_supp) {
            summary.class_mask |= kCandidateSplitSaSupplementary;
        }
        if (summary.max_soft_clip >= kSoftClipSignalMin) {
            summary.class_mask |= kCandidateSoftClip;
        }
        if (summary.max_ins >= kLongInsertionSignalMin) {
            summary.class_mask |= kCandidateLongInsertion;
        }

        for (auto& point : local_points) {
            point.class_mask = summary.class_mask;
            bundle.points.push_back(point);
        }
        bundle.read_summaries[idx] = std::move(summary);
    }

    return bundle;
}

std::vector<CandidateWindow> build_density_windows(
    const std::vector<EvidencePoint>& evidence,
    int32_t bin_start,
    int32_t bin_end) {
    if (evidence.empty()) {
        return {};
    }

    int32_t min_pos = evidence.front().pos;
    int32_t max_pos = evidence.front().pos;
    for (const auto& point : evidence) {
        min_pos = std::min(min_pos, point.pos);
        max_pos = std::max(max_pos, point.pos);
    }

    const int32_t hist_start = std::min(min_pos, bin_start) - kHistPadding;
    const int32_t hist_end = std::max(max_pos, bin_end) + kHistPadding;
    const int32_t num_bins = std::max(1, ((hist_end - hist_start) / kHistBinSize) + 1);

    std::vector<double> hist(static_cast<size_t>(num_bins), 0.0);
    for (const auto& point : evidence) {
        int32_t idx = (point.pos - hist_start) / kHistBinSize;
        idx = std::max(0, std::min(num_bins - 1, idx));
        hist[static_cast<size_t>(idx)] += point.weight;
    }
    const std::vector<double> smooth = smooth_histogram(hist);

    std::vector<DensityPeak> peaks;
    for (int32_t i = 0; i < num_bins; ++i) {
        const double center = smooth[static_cast<size_t>(i)];
        if (center < kPeakMinWeight) {
            continue;
        }
        const double left = (i > 0) ? smooth[static_cast<size_t>(i - 1)] : center;
        const double right = (i + 1 < num_bins) ? smooth[static_cast<size_t>(i + 1)] : center;
        if (center >= left && center >= right) {
            DensityPeak peak;
            peak.pos = hist_start + (i * kHistBinSize) + (kHistBinSize / 2);
            peak.weight = center;
            peaks.push_back(peak);
        }
    }

    if (peaks.empty()) {
        std::vector<std::pair<int32_t, double>> pos_weights;
        pos_weights.reserve(evidence.size());
        for (const auto& point : evidence) {
            pos_weights.push_back({point.pos, point.weight});
        }
        const int32_t fallback_center = weighted_median_position(std::move(pos_weights));
        if (fallback_center >= 0) {
            DensityPeak fallback_peak;
            fallback_peak.pos = fallback_center;
            fallback_peak.weight = 1.0;
            peaks.push_back(fallback_peak);
        }
    }
    if (peaks.empty()) {
        return {};
    }

    std::sort(peaks.begin(), peaks.end(), [](const DensityPeak& a, const DensityPeak& b) {
        return a.pos < b.pos;
    });

    std::vector<DensityPeak> merged_peaks;
    for (const auto& peak : peaks) {
        if (merged_peaks.empty() ||
            (peak.pos - merged_peaks.back().pos) > kPeakMergeDistance) {
            merged_peaks.push_back(peak);
            continue;
        }

        DensityPeak& back = merged_peaks.back();
        const double total_weight = back.weight + peak.weight;
        if (total_weight > 0.0) {
            const double weighted_center =
                ((back.weight * static_cast<double>(back.pos)) +
                 (peak.weight * static_cast<double>(peak.pos))) / total_weight;
            back.pos = static_cast<int32_t>(std::llround(weighted_center));
        }
        back.weight = total_weight;
    }

    std::vector<CandidateWindow> windows;
    windows.reserve(merged_peaks.size());

    for (const auto& peak : merged_peaks) {
        int32_t center_idx = (peak.pos - hist_start) / kHistBinSize;
        center_idx = std::max(0, std::min(num_bins - 1, center_idx));
        const double center_weight = smooth[static_cast<size_t>(center_idx)];
        const double expand_threshold = std::max(0.50, center_weight * kWindowDropRatio);

        int32_t left = center_idx;
        int32_t right = center_idx;

        while (left > 0 &&
               (center_idx - left) < kWindowMaxExpandBins &&
               smooth[static_cast<size_t>(left - 1)] >= expand_threshold) {
            --left;
        }
        while ((right + 1) < num_bins &&
               (right - center_idx) < kWindowMaxExpandBins &&
               smooth[static_cast<size_t>(right + 1)] >= expand_threshold) {
            ++right;
        }

        int32_t start = hist_start + (left * kHistBinSize);
        int32_t end = hist_start + ((right + 1) * kHistBinSize);
        if ((end - start) < kMinWindowSpan) {
            start = peak.pos - (kMinWindowSpan / 2);
            end = start + kMinWindowSpan;
        }

        start = std::max(bin_start, start);
        end = std::min(bin_end, end);
        if (end <= start) {
            continue;
        }

        CandidateWindow window;
        window.start = start;
        window.end = end;
        window.center = std::max(window.start, std::min(window.end - 1, peak.pos));
        window.peak_weight = std::max(peak.weight, center_weight);
        windows.push_back(window);
    }

    if (windows.empty()) {
        return windows;
    }

    std::sort(windows.begin(), windows.end(), [](const CandidateWindow& a, const CandidateWindow& b) {
        if (a.start != b.start) {
            return a.start < b.start;
        }
        return a.end < b.end;
    });

    std::vector<CandidateWindow> merged_windows;
    for (const auto& window : windows) {
        if (merged_windows.empty() ||
            window.start > (merged_windows.back().end + kWindowMergeGap)) {
            merged_windows.push_back(window);
            continue;
        }

        CandidateWindow& back = merged_windows.back();
        const double total_weight = back.peak_weight + window.peak_weight;
        if (total_weight > 0.0) {
            const double center =
                ((back.peak_weight * static_cast<double>(back.center)) +
                 (window.peak_weight * static_cast<double>(window.center))) / total_weight;
            back.center = static_cast<int32_t>(std::llround(center));
        }
        back.start = std::min(back.start, window.start);
        back.end = std::max(back.end, window.end);
        back.peak_weight = std::max(back.peak_weight, window.peak_weight);
    }

    return merged_windows;
}

}  // namespace

std::vector<ComponentCall> LinearBinComponentModule::build(
    const std::vector<BamRecordPtr>& bin_records,
    const std::string& chrom,
    int32_t tid,
    int32_t bin_start,
    int32_t bin_end) const {
    if (bin_records.empty()) {
        return {};
    }

    const EvidenceBundle evidence_bundle = extract_evidence_points(bin_records, tid);
    const auto windows = build_density_windows(evidence_bundle.points, bin_start, bin_end);
    if (windows.empty()) {
        return {};
    }

    std::vector<std::vector<const EvidencePoint*>> evidence_by_read(bin_records.size());
    for (const auto& point : evidence_bundle.points) {
        if (point.read_index < evidence_by_read.size()) {
            evidence_by_read[point.read_index].push_back(&point);
        }
    }

    std::vector<int32_t> assigned_window(bin_records.size(), -1);
    std::vector<int32_t> assigned_count(windows.size(), 0);
    std::vector<int32_t> ambiguous_count(windows.size(), 0);

    for (size_t read_idx = 0; read_idx < evidence_by_read.size(); ++read_idx) {
        const auto& read_points = evidence_by_read[read_idx];
        if (read_points.empty()) {
            continue;
        }

        int32_t best_window = -1;
        double best_score = 0.0;
        double second_score = 0.0;

        for (size_t wi = 0; wi < windows.size(); ++wi) {
            const auto& window = windows[wi];
            double score = 0.0;
            for (const EvidencePoint* point : read_points) {
                if (!point) {
                    continue;
                }
                const int32_t dist = std::abs(point->pos - window.center);
                const double proximity =
                    (point->pos >= window.start && point->pos <= window.end)
                    ? 1.0
                    : std::exp(-static_cast<double>(dist) / kReadAssignDecayBp);
                score += point->weight * proximity;
            }

            if (score > best_score) {
                second_score = best_score;
                best_score = score;
                best_window = static_cast<int32_t>(wi);
            } else if (score > second_score) {
                second_score = score;
            }
        }

        if (best_window < 0 || best_score < kReadAssignMinScore) {
            continue;
        }
        if (second_score >= best_score * kReadAmbiguousRatio) {
            ambiguous_count[static_cast<size_t>(best_window)] += 1;
            continue;
        }

        assigned_window[read_idx] = best_window;
        assigned_count[static_cast<size_t>(best_window)] += 1;
    }

    std::vector<ComponentCall> components;
    components.reserve(windows.size());

    for (size_t wi = 0; wi < windows.size(); ++wi) {
        if (assigned_count[wi] <= 0) {
            continue;
        }

        ComponentCall call;
        call.chrom = chrom;
        call.tid = tid;
        call.bin_start = windows[wi].start;
        call.bin_end = windows[wi].end;
        call.anchor_pos = windows[wi].center;
        call.peak_weight = windows[wi].peak_weight;

        const int32_t denom = assigned_count[wi] + ambiguous_count[wi];
        if (denom > 0) {
            call.ambiguous_frac = static_cast<double>(ambiguous_count[wi]) / static_cast<double>(denom);
        }

        std::vector<std::pair<int32_t, double>> anchor_points;
        for (size_t read_idx = 0; read_idx < assigned_window.size(); ++read_idx) {
            if (assigned_window[read_idx] != static_cast<int32_t>(wi)) {
                continue;
            }

            call.read_indices.push_back(read_idx);
            if (read_idx >= evidence_bundle.read_summaries.size()) {
                continue;
            }

            const auto& summary = evidence_bundle.read_summaries[read_idx];
            if (summary.max_soft_clip >= kSoftClipSignalMin) {
                call.soft_clip_read_indices.push_back(read_idx);
            }
            if (summary.has_sa_or_supp) {
                call.split_sa_read_indices.push_back(read_idx);
            }
            if (summary.max_ins >= kLongInsertionSignalMin) {
                call.insertion_read_indices.push_back(read_idx);
            }
        }

        for (const auto& point : evidence_bundle.points) {
            if (point.read_index >= assigned_window.size() ||
                assigned_window[point.read_index] != static_cast<int32_t>(wi)) {
                continue;
            }

            anchor_points.push_back({point.pos, point.weight});
            switch (point.kind) {
                case EvidenceKind::kSoftClip:
                    call.evidence_soft_clip_count += 1;
                    break;
                case EvidenceKind::kIndel:
                    call.evidence_indel_count += 1;
                    break;
                case EvidenceKind::kSAHint:
                    call.evidence_sa_hint_count += 1;
                    break;
            }

            BreakpointCandidate candidate;
            candidate.chrom = chrom;
            candidate.pos = point.pos;
            candidate.read_index = point.read_index;
            candidate.class_mask = point.class_mask;
            if (point.read_index < evidence_bundle.read_summaries.size()) {
                const auto& summary = evidence_bundle.read_summaries[point.read_index];
                candidate.read_id = summary.read_id;
                candidate.is_reverse = summary.is_reverse;
                candidate.anchor_len = summary.max_soft_clip;
            }
            if (point.kind == EvidenceKind::kSoftClip) {
                candidate.clip_len = point.signal_len;
            } else if (point.kind == EvidenceKind::kIndel) {
                candidate.ins_len = point.signal_len;
            }
            call.breakpoint_candidates.push_back(std::move(candidate));
        }

        if (!anchor_points.empty()) {
            const int32_t anchor = weighted_median_position(std::move(anchor_points));
            if (anchor >= 0) {
                call.anchor_pos = anchor;
            }
        }

        if (!call.read_indices.empty()) {
            components.push_back(std::move(call));
        }
    }

    return components;
}

Pipeline::Pipeline(PipelineConfig config, std::unique_ptr<BamStreamReader> bam_reader)
    : config_(std::move(config)),
      bam_reader_(std::move(bam_reader)),
      gate1_module_(),
      component_module_(),
      ins_fragment_module_(config_),
      te_classifier_module_(config_),
      anchor_locked_module_(config_) {}

PipelineResult Pipeline::run() const {
    if (!bam_reader_ || !bam_reader_->is_valid()) {
        throw std::runtime_error("BAM reader is not valid");
    }

    if (config_.enable_parallel) {
        return run_parallel();
    }
    return run_streaming();
}

PipelineResult Pipeline::run_streaming() const {
    PipelineResult result;
    StreamingState state;

    const auto progress_cb = [this](int64_t processed, int32_t tid) {
        std::cerr << "[Pipeline] processed=" << processed << " current_tid=" << tid << '\n';
        return true;
    };

    result.total_reads = bam_reader_->stream(
        [this, &state, &result](BamRecordPtr&& record) {
            consume_record(std::move(record), state, result);
        },
        progress_cb,
        config_.progress_interval);

    flush_current_bin(state, result);
    return result;
}

PipelineResult Pipeline::run_parallel() const {
    PipelineResult result;
    SafeQueue<std::vector<BamRecordPtr>> batch_queue;

    PipelineResult worker_result;
    std::thread worker([this, &batch_queue, &worker_result]() {
        StreamingState worker_state;
        std::vector<BamRecordPtr> batch;
        while (batch_queue.pop(batch)) {
            for (auto& record : batch) {
                consume_record(std::move(record), worker_state, worker_result);
            }
            batch.clear();
        }
        flush_current_bin(worker_state, worker_result);
    });

    std::vector<BamRecordPtr> current_batch;
    current_batch.reserve(config_.batch_size);

    const auto progress_cb = [this](int64_t processed, int32_t tid) {
        std::cerr << "[Pipeline] processed=" << processed << " current_tid=" << tid << '\n';
        return true;
    };

    const int64_t total_reads = bam_reader_->stream(
        [this, &current_batch, &batch_queue](BamRecordPtr&& record) {
            current_batch.push_back(std::move(record));
            if (current_batch.size() >= config_.batch_size) {
                batch_queue.push(std::move(current_batch));
                current_batch.clear();
                current_batch.reserve(config_.batch_size);
            }
        },
        progress_cb,
        config_.progress_interval);

    if (!current_batch.empty()) {
        batch_queue.push(std::move(current_batch));
    }

    batch_queue.close();
    worker.join();

    result = std::move(worker_result);
    result.total_reads = total_reads;
    return result;
}

void Pipeline::consume_record(
    BamRecordPtr&& record,
    StreamingState& state,
    PipelineResult& result) const {
    if (!record) {
        return;
    }

    ReadView view(record.get());
    if (!gate1_module_.pass_preliminary(view)) {
        return;
    }

    result.gate1_passed += 1;

    while (!state.active_window.empty()) {
        const auto& front = state.active_window.front();
        const bool different_tid = (view.tid() != front.tid);
        const bool beyond_window = (!different_tid && (view.pos() - front.pos > config_.window_size));
        if (!different_tid && !beyond_window) {
            break;
        }
        state.active_window.pop_front();
    }
    state.active_window.push_back({view.tid(), view.pos()});

    const int32_t bin_index = (config_.bin_size > 0) ? (view.pos() / config_.bin_size) : 0;

    if (state.current_tid < 0) {
        state.current_tid = view.tid();
        state.current_bin_index = bin_index;
    }

    if (view.tid() != state.current_tid || bin_index != state.current_bin_index) {
        flush_current_bin(state, result);
        state.current_tid = view.tid();
        state.current_bin_index = bin_index;
    }

    state.current_bin_records.push_back(std::move(record));
}

void Pipeline::flush_current_bin(
    StreamingState& state,
    PipelineResult& result) const {
    if (state.current_bin_records.empty()) {
        return;
    }

    process_bin_records(
        std::move(state.current_bin_records),
        state.current_tid,
        state.current_bin_index,
        result);

    state.current_bin_records.clear();
}

std::vector<LocusEvidence> Pipeline::collect_evidence(
    const ComponentCall& component,
    const std::vector<BamRecordPtr>& bin_records) const {
    std::vector<LocusEvidence> evidence;
    evidence.reserve(component.read_indices.size());

    const int32_t center = component.anchor_pos;

    for (size_t idx : component.read_indices) {
        if (idx >= bin_records.size() || !bin_records[idx]) {
            continue;
        }

        ReadView view(bin_records[idx].get());
        LocusEvidence e;
        e.tid = view.tid();
        e.pos = view.pos();
        e.read_index = idx;

        const int32_t distance = std::abs(view.pos() - center);
        e.normalized_score = std::max(0.0, 1.0 - static_cast<double>(distance) / 10000.0);
        evidence.push_back(e);
    }

    return evidence;
}

AssemblyCall Pipeline::assemble_component(
    const ComponentCall& component,
    const std::vector<BamRecordPtr>& bin_records) const {
    AssemblyCall call;
    call.tid = component.tid;
    call.pos = component.anchor_pos;

    if (!component.read_indices.empty()) {
        const size_t idx = component.read_indices.front();
        if (idx < bin_records.size() && bin_records[idx]) {
            ReadView view(bin_records[idx].get());
            std::string seq = view.decode_sequence();
            if (seq.size() > 300) {
                seq.resize(300);
            }
            call.consensus = std::move(seq);
        }
    }

    return call;
}

PlaceabilityReport Pipeline::score_placeability(
    const AssemblyCall& assembly,
    const std::vector<LocusEvidence>& evidence) const {
    PlaceabilityReport report;
    report.tid = assembly.tid;
    report.pos = assembly.pos;
    report.support_reads = static_cast<int32_t>(evidence.size());

    double score_sum = 0.0;
    for (const auto& row : evidence) {
        score_sum += row.normalized_score;
    }

    report.delta_score = report.support_reads > 0
        ? (score_sum / static_cast<double>(report.support_reads)) * 100.0
        : 0.0;
    report.tier = classify_tier(report.delta_score);
    return report;
}

GenotypeCall Pipeline::genotype_call(
    const AssemblyCall& assembly,
    const PlaceabilityReport& placeability) const {
    GenotypeCall call;
    call.tid = assembly.tid;
    call.pos = assembly.pos;

    const double af = std::clamp(static_cast<double>(placeability.support_reads) / 20.0, 0.0, 1.0);
    call.af = af;
    call.gq = static_cast<int32_t>(std::round(placeability.delta_score));

    if (af >= 0.8) {
        call.genotype = "1/1";
    } else if (af >= 0.2) {
        call.genotype = "0/1";
    } else {
        call.genotype = "0/0";
    }

    return call;
}

void Pipeline::process_bin_records(
    std::vector<BamRecordPtr>&& bin_records,
    int32_t tid,
    int32_t bin_index,
    PipelineResult& result) const {
    if (bin_records.empty()) {
        return;
    }

    result.processed_bins += 1;

    const int32_t bin_start = bin_index * config_.bin_size;
    const int32_t bin_end = bin_start + config_.bin_size;
    const std::string chrom = bam_reader_->chromosome_name(tid);

    auto components = component_module_.build(bin_records, chrom, tid, bin_start, bin_end);
    result.built_components += static_cast<int64_t>(components.size());

    for (const auto& component : components) {
        std::vector<InsertionFragment> fragments = ins_fragment_module_.extract(component, bin_records);
        std::vector<FragmentTEHit> hits;
        ClusterTECall te_call;

        if (te_classifier_module_.is_enabled()) {
            hits = te_classifier_module_.classify(fragments);
            te_call = te_classifier_module_.vote_cluster(hits);
        }
        // Only emit TE-passing calls.
        if (!te_call.passed || te_call.te_name.empty()) {
            continue;
        }

        AnchorLockedReport anchor_report;
        if (anchor_locked_module_.is_enabled()) {
            anchor_report = anchor_locked_module_.resolve(component, fragments, hits, te_call);
        }

        auto evidence = collect_evidence(component, bin_records);
        result.evidence_rows += static_cast<int64_t>(evidence.size());

        AssemblyCall assembly = assemble_component(component, bin_records);
        result.assembled_calls += 1;

        PlaceabilityReport placeability = score_placeability(assembly, evidence);
        result.placeability_calls += 1;

        GenotypeCall genotype = genotype_call(assembly, placeability);
        result.genotype_calls += 1;

        FinalCall call;
        call.chrom = component.chrom;
        call.tid = assembly.tid;
        call.pos = assembly.pos;
        call.window_start = component.bin_start;
        call.window_end = component.bin_end;
        call.te_name = te_call.te_name;
        call.te_vote_fraction = te_call.vote_fraction;
        call.te_median_identity = te_call.median_identity;
        call.te_fragment_count = te_call.fragment_count;
        call.te_theta = theta_tag(anchor_report.theta0);
        call.te_mad_fwd = anchor_report.mad_fwd;
        call.te_mad_rev = anchor_report.mad_rev;
        call.te_breakpoint_core = anchor_report.te_breakpoint_core;
        call.te_breakpoint_window_start = anchor_report.te_breakpoint_window_start;
        call.te_breakpoint_window_end = anchor_report.te_breakpoint_window_end;
        call.te_core_candidates = anchor_report.core_candidate_count;
        call.te_core_set = anchor_report.core_set_count;
        call.te_split_sa_core_frac = anchor_report.split_sa_core_frac;
        call.te_ref_junc_pos_min = anchor_report.ref_junc_pos_min;
        call.te_ref_junc_pos_max = anchor_report.ref_junc_pos_max;
        if (!anchor_report.enabled) {
            call.te_qc = "DISABLED";
        } else if (anchor_report.fail_theta_uncertain) {
            call.te_qc = "FAIL_THETA_UNCERTAIN";
        } else if (!anchor_report.has_result) {
            call.te_qc = "NO_CORE_RESULT";
        } else {
            call.te_qc = "PASS";
        }

        call.tier = placeability.tier;
        call.support_reads = placeability.support_reads;
        call.genotype = genotype.genotype;
        call.af = genotype.af;
        call.gq = genotype.gq;
        result.final_calls.push_back(std::move(call));
    }
}

std::unique_ptr<Pipeline> build_default_pipeline(const PipelineConfig& config) {
    auto reader = make_bam_reader(config.bam_path, config.bam_threads);
    return std::make_unique<Pipeline>(config, std::move(reader));
}

}  // namespace placer
