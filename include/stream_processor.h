#ifndef PLACER_STREAM_PROCESSOR_H
#define PLACER_STREAM_PROCESSOR_H

#include "bam_reader.h"
#include "window_buffer.h"
#include "component_builder.h"
#include "task_queue.h"
#include <atomic>
#include <chrono>
#include <functional>
#include <optional>
#include <vector>

namespace placer {

/**
 * StreamProcessor: True streaming processor for TE detection
 *
 * Key principles:
 * - No full genome scan - single pass through BAM
 * - Windows trigger and process as reads arrive
 * - Memory bounded by window_size and trigger frequency
 * - Reference genome determines chromosome count
 */
class StreamProcessor {
public:
    struct Config {
        // Window config
        int window_size = 10000;
        int window_step = 5000;
        int max_priority_reads = 50;
        int max_normal_reads = 200;
        int min_clip_bp = 50;
        int min_sa_reads = 3;

        // Component building
        int min_density = 0;
        int max_cluster_span = 200;
        int max_recursive_depth = 4;

        // Progress reporting
        bool verbose = true;
        int progress_interval = 100000;
    };

    struct Result {
        int64_t total_reads = 0;
        int64_t triggered_windows = 0;
        int64_t built_components = 0;
        double elapsed_seconds = 0.0;
        std::string chrom_name;  // Reference genome's chromosome count
    };

    explicit StreamProcessor(Config config);

    /**
     * Process BAM in streaming mode
     * Windows trigger immediately when threshold exceeded
     *
     * @param reader BAM reader (must be valid)
     * @param on_window_triggered Called when a window triggers (before flush)
     * @return Processing result
     */
    Result process(
        BamReader& reader,
        std::function<void(std::unique_ptr<Window>)> on_window_triggered = nullptr);

    /**
     * Get the chromosome count from reference genome
     */
    int32_t get_chromosome_count() const { return chromosome_count_; }

    /**
     * Get chromosome name from index
     */
    std::string get_chromosome_name(int32_t tid) const;

private:
    Config config_;
    WindowBuffer window_buffer_;
    ComponentBuilder component_builder_;
    std::atomic<int64_t> processed_reads_{0};
    std::atomic<int64_t> triggered_count_{0};
    std::vector<std::string> chromosome_names_;
    int32_t chromosome_count_ = 0;

    // Progress reporting
    std::chrono::steady_clock::time_point start_time_;
    bool verbose_ = true;

    bool progress_callback_(int64_t processed, int32_t current_chrom);
};

/**
 * Simple in-memory stream processor for quick iteration
 * Accumulates components during streaming, processes at the end
 *
 * Use this when you need components for further processing
 */
class AccumulatingStreamProcessor {
public:
    struct Config {
        int window_size = 10000;
        int window_step = 5000;
        int min_clip_bp = 50;
        int min_sa_reads = 3;
        int min_density = 0;
        int max_cluster_span = 200;
        int max_recursive_depth = 4;
        bool verbose = true;
        int progress_interval = 500000;
    };

    struct Result {
        int64_t total_reads = 0;
        int64_t triggered_windows = 0;
        int64_t built_components = 0;
        std::vector<Component> components;
        double elapsed_seconds = 0.0;
    };

    explicit AccumulatingStreamProcessor(Config config);

    Result process(BamReader& reader);

    const std::vector<std::string>& get_chromosome_names() const {
        return chromosome_names_;
    }

private:
    Config config_;
    WindowBuffer window_buffer_;
    ComponentBuilder component_builder_;
    std::vector<std::string> chromosome_names_;
    std::chrono::steady_clock::time_point start_time_;
    bool verbose_ = true;
};

/**
 * Estimate processing time based on BAM size
 */
inline double estimate_processing_speed(int64_t reads_per_second) {
    // Typical long-read throughput: 1M-5M reads/hour for ONT
    // This is just informational
    return reads_per_second;
}

}  // namespace placer

#endif  // PLACER_STREAM_PROCESSOR_H
