#include "local_interval_cache.h"

#include <algorithm>

namespace placer {

std::vector<CanonicalLocalInterval> build_canonical_local_intervals(
    const std::vector<LocalIntervalRequest>& requests,
    int32_t merge_gap_bp) {
    std::vector<LocalIntervalRequest> ordered = requests;
    std::sort(ordered.begin(), ordered.end(), [](const auto& a, const auto& b) {
        if (a.chrom != b.chrom) {
            return a.chrom < b.chrom;
        }
        if (a.start != b.start) {
            return a.start < b.start;
        }
        return a.end < b.end;
    });

    std::vector<CanonicalLocalInterval> out;
    for (const auto& req : ordered) {
        if (out.empty() ||
            out.back().chrom != req.chrom ||
            req.start > (out.back().end + merge_gap_bp)) {
            CanonicalLocalInterval interval;
            interval.chrom = req.chrom;
            interval.start = req.start;
            interval.end = req.end;
            interval.request_ids.push_back(req.request_id);
            out.push_back(std::move(interval));
            continue;
        }
        CanonicalLocalInterval& back = out.back();
        back.end = std::max(back.end, req.end);
        back.request_ids.push_back(req.request_id);
    }
    return out;
}

LocalIntervalProjection project_cached_interval_reads(
    const LocalIntervalRequest& request,
    const std::vector<LocalIntervalCacheEntry>& cache_entries) {
    LocalIntervalProjection out;
    for (const auto& entry : cache_entries) {
        if (entry.interval.chrom != request.chrom) {
            continue;
        }
        if (std::find(
                entry.interval.request_ids.begin(),
                entry.interval.request_ids.end(),
                request.request_id) == entry.interval.request_ids.end()) {
            continue;
        }
        for (size_t i = 0; i < entry.records.size() && i < entry.read_spans.size(); ++i) {
            const ReadReferenceSpan& span = entry.read_spans[i];
            if (!span.valid) {
                continue;
            }
            if (span.end <= request.start || span.start >= request.end) {
                continue;
            }
            out.records.push_back(entry.records[i]);
            out.read_spans.push_back(span);
        }
        break;
    }
    return out;
}

double local_interval_reuse_ratio(const LocalIntervalReuseStats& stats) {
    if (stats.canonical_interval_count == 0) {
        return 1.0;
    }
    return static_cast<double>(stats.request_count) /
        static_cast<double>(stats.canonical_interval_count);
}

}  // namespace placer
