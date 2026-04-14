#ifndef PLACER_LOCAL_INTERVAL_CACHE_H
#define PLACER_LOCAL_INTERVAL_CACHE_H

#include "pipeline.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace placer {

struct LocalIntervalRequest {
    std::string chrom;
    int32_t start = -1;
    int32_t end = -1;
    size_t request_id = 0;
};

struct CanonicalLocalInterval {
    std::string chrom;
    int32_t start = -1;
    int32_t end = -1;
    std::vector<size_t> request_ids;
};

struct LocalIntervalCacheEntry {
    CanonicalLocalInterval interval;
    std::vector<BamRecordPtr> owned_records;
    std::vector<const bam1_t*> records;
    std::vector<ReadReferenceSpan> read_spans;
};

struct LocalIntervalProjection {
    std::vector<const bam1_t*> records;
    std::vector<ReadReferenceSpan> read_spans;
};

struct LocalIntervalReuseStats {
    size_t request_count = 0;
    size_t canonical_interval_count = 0;
};

std::vector<CanonicalLocalInterval> build_canonical_local_intervals(
    const std::vector<LocalIntervalRequest>& requests,
    int32_t merge_gap_bp);

LocalIntervalProjection project_cached_interval_reads(
    const LocalIntervalRequest& request,
    const std::vector<LocalIntervalCacheEntry>& cache_entries);

double local_interval_reuse_ratio(const LocalIntervalReuseStats& stats);

}  // namespace placer

#endif
