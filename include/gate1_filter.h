#ifndef PLACER_GATE1_FILTER_H
#define PLACER_GATE1_FILTER_H

#include "gate1.h"
#include "window_buffer.h"
#include "task_queue.h"
#include "bam_reader.h"
#include <memory>
#include <functional>
#include <atomic>
#include <vector>

namespace placer {

/**
 * Gate1FilterConfig: Configuration for Gate 1 filtering
 */
struct Gate1FilterConfig {
    // TE index configuration
    std::string te_fasta_path;  // Path to TE library FASTA

    // Filter mode
    bool enabled = true;         // Enable/disable Gate 1 filtering

    // Sampling control (for performance in high-depth regions)
    int max_reads_per_window = 100;  // Max reads to evaluate per triggered window
    double sample_ratio = 1.0;         // Sample ratio for high-depth windows [0, 1]

    // Reporting
    bool report_filtered = false;  // Log reads that were filtered out
};

/**
 * Gate1Filter: TE-proxy filtering for triggered windows
 *
 * Industrial-grade implementation:
 * - Zero-allocation hot path (passes())
 * - Move semantics throughout to avoid deep copies
 * - Configurable sampling for high-depth regions
 * - Atomic counters for thread-safe statistics
 */
class Gate1Filter {
public:
    explicit Gate1Filter(Gate1FilterConfig config);

    /**
     * Check if filtering is enabled
     */
    bool is_enabled() const { return config_.enabled && index_ != nullptr; }

    /**
     * Fast path: Check if a single read passes Gate 1
     * Zero allocation - uses string_view from ReadSketch
     *
     * Returns: true if read should be kept
     */
    bool passes(const ReadSketch& read) const {
        if (!is_enabled()) return true;
        if (!index_) return true;
        return gate1_->passes(read);
    }

    /**
     * Statistics update methods (atomic, thread-safe)
     */
    void increment_total() { total_reads_.fetch_add(1, std::memory_order_relaxed); }
    void increment_passed() { total_passed_.fetch_add(1, std::memory_order_relaxed); }
    void batch_update_stats(int64_t processed, int64_t passed, int64_t filtered);

    /**
     * Statistics accessors
     */
    int64_t get_total_reads_processed() const { return total_reads_.load(std::memory_order_relaxed); }
    int64_t get_total_reads_passed() const { return total_passed_.load(std::memory_order_relaxed); }
    int64_t get_total_reads_filtered() const { return total_filtered_.load(std::memory_order_relaxed); }
    int64_t get_total_windows_processed() const { return total_windows_.load(std::memory_order_relaxed); }

    /**
     * Get configuration
     */
    const Gate1FilterConfig& config() const { return config_; }

    /**
     * Build the TE index from FASTA
     * Must be called before filtering
     */
    bool build_index();

    /**
     * Check if index is ready
     */
    bool is_index_ready() const { return index_ != nullptr; }

private:
    Gate1FilterConfig config_;
    std::shared_ptr<HashTEIndex> index_;
    std::unique_ptr<Gate1> gate1_;

    // Thread-safe statistics (relaxed ordering for performance)
    std::atomic<int64_t> total_reads_{0};
    std::atomic<int64_t> total_passed_{0};
    std::atomic<int64_t> total_filtered_{0};
    std::atomic<int64_t> total_windows_{0};
};

/**
 * WindowProcessor: Orchestrates window processing with optional Gate 1 filtering
 *
 * Industrial-grade implementation:
 * - Zero-copy: Direct ownership transfer from Window to result
 * - Move semantics: ReadSketch objects are moved, not copied
 * - Batch submission: Efficient task queue integration
 *
 * Data flow:
 *   WindowBuffer.seal_and_flush() → std::unique_ptr<Window>
 *       → WindowProcessor::process_triggered_window()
 *           → Gate1Filter (optional filtering)
 *               → TaskQueue.submit_serialized() (by move)
 */
class WindowProcessor {
public:
    struct Config {
        TaskQueue* task_queue = nullptr;       // Required: task queue for output
        Gate1Filter* gate1_filter = nullptr;  // Optional: Gate 1 filter (nullptr = pass-through)
        int32_t batch_size = 10;               // Batch size for task submission
    };

    explicit WindowProcessor(Config config);

    /**
     * Process a triggered window
     *
     * Performance characteristics:
     * - Takes ownership of Window via std::unique_ptr
     * - Moves ReadSketch objects directly (no intermediate copies)
     * - Single pass through reads
     *
     * @param window Owned window to process (moved from)
     */
    void process_triggered_window(std::unique_ptr<Window> window);

    /**
     * Process multiple triggered windows
     * @param windows Vector of owned windows (moved from)
     */
    void process_triggered_windows(std::vector<std::unique_ptr<Window>>& windows);

    /**
     * Statistics accessors
     */
    int64_t get_total_windows_processed() const { return windows_processed_.load(std::memory_order_relaxed); }
    int64_t get_total_reads_submitted() const { return reads_submitted_.load(std::memory_order_relaxed); }
    int64_t get_total_reads_filtered() const { return reads_filtered_.load(std::memory_order_relaxed); }

private:
    Config config_;

    // Thread-safe statistics
    std::atomic<int64_t> windows_processed_{0};
    std::atomic<int64_t> reads_submitted_{0};
    std::atomic<int64_t> reads_filtered_{0};
};

}  // namespace placer

#endif  // PLACER_GATE1_FILTER_H
