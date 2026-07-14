template <typename T>
class SafeQueue {
public:
    explicit SafeQueue(size_t max_size = 0) : max_size_(max_size) {}

    void push(T&& value) {
        std::unique_lock<std::mutex> lock(mu_);
        if (max_size_ > 0) {
            not_full_cv_.wait(lock, [this]() { return finished_ || queue_.size() < max_size_; });
        }
        if (finished_) {
            return;
        }
        queue_.push(std::move(value));
        lock.unlock();
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
        lock.unlock();
        if (max_size_ > 0) {
            not_full_cv_.notify_one();
        }
        return true;
    }

    void close() {
        {
            std::lock_guard<std::mutex> lock(mu_);
            finished_ = true;
        }
        cv_.notify_all();
        if (max_size_ > 0) {
            not_full_cv_.notify_all();
        }
    }

private:
    std::queue<T> queue_;
    std::mutex mu_;
    std::condition_variable cv_;
    std::condition_variable not_full_cv_;
    bool finished_ = false;
    size_t max_size_ = 0;
};

struct BinTask {
    std::vector<const bam1_t*> records;
    std::array<SharedRecordBatch, static_cast<size_t>(kCrossBinContextBins) + 1> owners;
    size_t owner_count = 0;
    int32_t tid = -1;
    int32_t bin_index = -1;
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

int32_t compute_ref_end(const ReadView& read) {
    const uint32_t* cigar = read.cigar();
    const int32_t n_cigar = read.n_cigar();
    if (!cigar || n_cigar <= 0) {
        return read.pos();
    }
    int32_t ref_pos = read.pos();
    for (int32_t i = 0; i < n_cigar; ++i) {
        const int op = bam_cigar_op(cigar[i]);
        const int32_t len = static_cast<int32_t>(bam_cigar_oplen(cigar[i]));
        if ((bam_cigar_type(op) & 2) != 0) {
            ref_pos += len;
        }
    }
    return ref_pos;
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
    const std::vector<const bam1_t*>& bin_records,
    int32_t expected_tid) {
    EvidenceBundle bundle;
    bundle.read_summaries.resize(bin_records.size());

    for (size_t idx = 0; idx < bin_records.size(); ++idx) {
        if (!bin_records[idx]) {
            continue;
        }

        ReadView view(bin_records[idx]);
        if (view.tid() != expected_tid) {
            continue;
        }
        if ((view.flag() & BAM_FSUPPLEMENTARY) != 0) {
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
                point.class_mask = kCandidateLongInsertion;
                local_points.push_back(point);
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
            point.class_mask = kCandidateSoftClip;
            local_points.push_back(point);
        }
        if (trailing_soft >= kSoftClipSignalMin) {
            EvidencePoint point;
            point.read_index = idx;
            point.pos = ref_end;
            point.weight = 1.00 + (static_cast<double>(std::min(trailing_soft, 500)) / 180.0);
            point.kind = EvidenceKind::kSoftClip;
            point.signal_len = trailing_soft;
            point.class_mask = kCandidateSoftClip;
            local_points.push_back(point);
        }

        if (summary.has_sa_or_supp) {
            const SaHintSide hint_side = choose_sa_hint_side(leading_soft, trailing_soft);
            if (hint_side != SaHintSide::kNone) {
                EvidencePoint hint;
                hint.read_index = idx;
                hint.pos = (hint_side == SaHintSide::kLeft) ? view.pos() : ref_end;
                hint.weight = kSaHintWeight;
                hint.kind = EvidenceKind::kSAHint;
                hint.class_mask = kCandidateSplitSaSupplementary;
                local_points.push_back(hint);
            }
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

    // Preserve direct long-indel loci even when nearby clip-driven density
    // dominates the smoothed histogram. These points are breakpoint-specific
    // enough to seed their own candidate windows.
    std::vector<const EvidencePoint*> indel_points;
    indel_points.reserve(evidence.size());
    for (const auto& point : evidence) {
        if (point.kind == EvidenceKind::kIndel) {
            indel_points.push_back(&point);
        }
    }
    std::sort(indel_points.begin(), indel_points.end(), [](const EvidencePoint* a, const EvidencePoint* b) {
        if (a == nullptr || b == nullptr) {
            return a < b;
        }
        return a->pos < b->pos;
    });
    size_t indel_cluster_begin = 0;
    while (indel_cluster_begin < indel_points.size()) {
        size_t indel_cluster_end = indel_cluster_begin + 1;
        while (indel_cluster_end < indel_points.size() &&
               std::abs(indel_points[indel_cluster_end]->pos -
                        indel_points[indel_cluster_end - 1]->pos) <= kPeakMergeDistance) {
            ++indel_cluster_end;
        }

        double total_weight = 0.0;
        double weighted_pos = 0.0;
        for (size_t i = indel_cluster_begin; i < indel_cluster_end; ++i) {
            total_weight += indel_points[i]->weight;
            weighted_pos += indel_points[i]->weight * static_cast<double>(indel_points[i]->pos);
        }
        if (total_weight >= kPeakMinWeight) {
            DensityPeak peak;
            peak.pos = static_cast<int32_t>(std::llround(weighted_pos / total_weight));
            peak.weight = total_weight;
            peaks.push_back(peak);
        }
        indel_cluster_begin = indel_cluster_end;
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

        const int32_t relaxed_bin_start = bin_start - kWindowBinSlackBp;
        const int32_t relaxed_bin_end = bin_end + kWindowBinSlackBp;
        start = std::max(relaxed_bin_start, start);
        end = std::min(relaxed_bin_end, end);
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
