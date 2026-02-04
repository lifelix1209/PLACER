#include "gate1_filter.h"
#include <algorithm>

namespace placer {

Gate1Filter::Gate1Filter(Gate1FilterConfig config)
    : config_(std::move(config)) {}

bool Gate1Filter::build_index() {
    if (!config_.enabled || config_.te_fasta_path.empty()) {
        // Filtering disabled or no TE FASTA
        return true;
    }

    HashTEIndexConfig idx_config;
    auto index = HashTEIndex::build_from_fasta(config_.te_fasta_path, idx_config);
    if (!index) {
        return false;
    }

    index_ = std::move(index);
    gate1_ = std::make_unique<Gate1>(index_, Gate1Config());
    return true;
}

void Gate1Filter::batch_update_stats(int64_t processed, int64_t passed, int64_t filtered) {
    total_reads_.fetch_add(processed, std::memory_order_relaxed);
    total_passed_.fetch_add(passed, std::memory_order_relaxed);
    total_filtered_.fetch_add(filtered, std::memory_order_relaxed);
}

// ============= WindowProcessor 实现 =============

WindowProcessor::WindowProcessor(Config config) : config_(std::move(config)) {}

void WindowProcessor::process_triggered_window(std::unique_ptr<Window> window) {
    if (!window || (window->priority_reads.empty() && window->normal_reads.empty())) {
        return;
    }

    // Result container: reserve space to avoid reallocations
    std::vector<ReadSketch> passing_reads;
    passing_reads.reserve(window->priority_reads.size() + window->normal_reads.size());

    // Statistics tracking
    int64_t processed = 0;
    int64_t passed = 0;

    if (config_.gate1_filter && config_.gate1_filter->is_enabled()) {
        auto* filter = config_.gate1_filter;

        // ============= Priority Reads =============
        // Priority reads (SA/MD events): process ALL, only filter by Gate 1
        // Note: These are moved directly, no intermediate copies
        for (auto& read : window->priority_reads) {
            ++processed;
            if (filter->passes(read)) {
                passing_reads.push_back(std::move(read));
                ++passed;
            }
        }

        // ============= Normal Reads =============
        // Normal reads: apply sampling if window is too large
        // Sampling is deterministic: process reads in order, stop when limit reached
        // This assumes WindowBuffer has already done stratified sampling during fill
        size_t max_normal = static_cast<size_t>(filter->config().max_reads_per_window);
        size_t normal_count = 0;

        for (auto& read : window->normal_reads) {
            if (normal_count >= max_normal) {
                // Sampling limit reached, skip remaining reads
                break;
            }

            ++processed;
            if (filter->passes(read)) {
                passing_reads.push_back(std::move(read));
                ++passed;
                ++normal_count;
            }
        }

        // Batch update statistics (reduce atomic contention)
        int64_t filtered = processed - passed;
        filter->batch_update_stats(processed, passed, filtered);
        reads_filtered_.fetch_add(filtered, std::memory_order_relaxed);
    } else {
        // ============= Pass-through Mode =============
        // No Gate 1 filtering: move all reads directly
        for (auto& read : window->priority_reads) {
            passing_reads.push_back(std::move(read));
        }
        for (auto& read : window->normal_reads) {
            passing_reads.push_back(std::move(read));
        }
        processed = window->priority_reads.size() + window->normal_reads.size();
        passed = processed;
    }

    // ============= Task Submission =============
    // Submit only if there are passing reads
    if (!passing_reads.empty() && config_.task_queue) {
        // Generate window ID
        std::string window_id = std::to_string(window->chrom_tid) + ":" +
                                std::to_string(window->start) + "-" +
                                std::to_string(window->end);

        // submit_serialized takes by value, which is move-optimized
        // The vector and its contents are moved, not copied
        bool submitted = config_.task_queue->submit_serialized(
            TaskType::COMPONENT_BUILD,
            window_id,
            std::move(passing_reads)
        );

        if (submitted) {
            reads_submitted_.fetch_add(passed, std::memory_order_relaxed);
        }
    }

    windows_processed_.fetch_add(1, std::memory_order_relaxed);
}

void WindowProcessor::process_triggered_windows(std::vector<std::unique_ptr<Window>>& windows) {
    for (auto& window : windows) {
        process_triggered_window(std::move(window));
    }
}

}  // namespace placer
