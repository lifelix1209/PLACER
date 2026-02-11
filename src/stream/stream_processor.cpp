#include "stream_processor.h"
#include "local_realign.h"
#include <iostream>
#include <iomanip>
#include <thread>

namespace placer {

StreamProcessor::StreamProcessor(Config config)
    : config_(std::move(config)),
      window_buffer_({
          config_.window_size,
          config_.window_step,
          config_.max_priority_reads,
          config_.max_normal_reads,
          config_.min_clip_bp,
          config_.min_sa_reads
      }),
      component_builder_({
          20,   // min_clip_len
          50,   // min_ins_len
          50,   // min_del_len
          50,   // cluster_gap
          config_.max_cluster_span,
          20,   // anchor_merge_distance
          20,   // min_mapq
          config_.min_density / 100.0,  // min_density (scaled)
          true,  // parse_all_sa_splits
          config_.max_recursive_depth,
          0.3,   // min_split_ratio
          3,     // min_anchors_per_cluster
          2      // min_reads_per_component
      }),
      verbose_(config_.verbose) {}

StreamProcessor::Result StreamProcessor::process(
    BamReader& reader,
    std::function<void(std::unique_ptr<Window>)> on_window_triggered) {

    Result result;
    start_time_ = std::chrono::steady_clock::now();
    processed_reads_ = 0;
    triggered_count_ = 0;

    // Get chromosome names from BAM header
    chromosome_count_ = reader.get_num_chromosomes();
    chromosome_names_.reserve(chromosome_count_);
    for (int32_t i = 0; i < chromosome_count_; ++i) {
        chromosome_names_.push_back(reader.get_chrom_name(i));
    }

    if (verbose_) {
        std::cout << "[StreamProcessor] Starting streaming with "
                  << chromosome_count_ << " chromosomes" << std::endl;
    }

    // Create the read callback that integrates with window buffer
    auto read_callback = [&](const ReadSketch& read) {
        window_buffer_.add_read(read);
        processed_reads_++;
    };

    // Create progress callback
    StreamProgressCallback progress_cb = nullptr;
    if (config_.progress_interval > 0) {
        progress_cb = [this](int64_t processed, int32_t current_chrom) -> bool {
            return progress_callback_(processed, current_chrom);
        };
    }

    // Stream through BAM
    int64_t total = reader.stream_with_progress(
        read_callback, progress_cb, config_.progress_interval);

    result.total_reads = total;
    result.elapsed_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start_time_).count();

    if (verbose_) {
        std::cout << "\n[StreamProcessor] Streaming complete: "
                  << result.total_reads << " reads in "
                  << result.elapsed_seconds << "s" << std::endl;
    }

    return result;
}

std::string StreamProcessor::get_chromosome_name(int32_t tid) const {
    if (tid >= 0 && tid < static_cast<int32_t>(chromosome_names_.size())) {
        return chromosome_names_[tid];
    }
    return "";
}

bool StreamProcessor::progress_callback_(int64_t processed, int32_t current_chrom) {
    double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start_time_).count();
    double rate = processed / elapsed;

    std::cout << "\r  Processed " << processed << " reads ("
              << std::fixed << std::setprecision(0) << rate << " reads/sec)";

    if (current_chrom >= 0 && current_chrom < chromosome_count_) {
        std::cout << " - chrom " << chromosome_names_[current_chrom];
    }
    std::cout << std::flush;

    return true;  // Continue streaming
}

// ============================================================================
// AccumulatingStreamProcessor
// ============================================================================

AccumulatingStreamProcessor::AccumulatingStreamProcessor(Config config)
    : config_(std::move(config)),
      window_buffer_({
          config_.window_size,
          config_.window_step,
          50,    // max_priority_reads
          200,   // max_normal_reads
          config_.min_clip_bp,
          config_.min_sa_reads
      }),
      component_builder_({
          20,   // min_clip_len
          50,   // min_ins_len
          50,   // min_del_len
          50,   // cluster_gap
          config_.max_cluster_span,
          20,   // anchor_merge_distance
          20,   // min_mapq
          config_.min_density / 100.0,
          true,
          config_.max_recursive_depth,
          0.3,
          3,
          2
      }),
      verbose_(config_.verbose) {}

AccumulatingStreamProcessor::Result AccumulatingStreamProcessor::process(
    BamReader& reader) {

    Result result;
    start_time_ = std::chrono::steady_clock::now();
    int64_t processed = 0;

    // Get chromosome names
    chromosome_names_.clear();
    int32_t num_chroms = reader.get_num_chromosomes();
    chromosome_names_.reserve(num_chroms);
    for (int32_t i = 0; i < num_chroms; ++i) {
        chromosome_names_.push_back(reader.get_chrom_name(i));
    }

    if (verbose_) {
        std::cout << "[AccumulatingStreamProcessor] Streaming with "
                  << num_chroms << " chromosomes" << std::endl;
    }

    // Callback for each read
    auto read_callback = [&](const ReadSketch& read) {
        window_buffer_.add_read(read);
        processed++;
    };

    // Progress callback
    StreamProgressCallback progress_cb = nullptr;
    if (config_.progress_interval > 0) {
        progress_cb = [&](int64_t p, int32_t tid) -> bool {
            double elapsed = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - start_time_).count();
            double rate = p / elapsed;

            std::cout << "\r  Loaded " << p << " reads ("
                      << std::fixed << std::setprecision(0)
                      << rate << " reads/sec)" << std::flush;

            if (tid >= 0 && tid < num_chroms) {
                std::cout << " - " << chromosome_names_[tid];
            }

            return true;
        };
    }

    // Stream BAM
    reader.stream_with_progress(read_callback, progress_cb, config_.progress_interval);
    result.total_reads = processed;

    if (verbose_) {
        std::cout << "\n  Loaded " << result.total_reads << " reads" << std::endl;
    }

    // Flush remaining windows and build components
    auto sealed_windows = window_buffer_.flush_current_chromosome();
    result.triggered_windows = sealed_windows.size();

    if (verbose_) {
        std::cout << "  Triggered " << result.triggered_windows << " windows" << std::endl;
    }

    // Build components from sealed windows
    for (auto& win : sealed_windows) {
        // Merge priority and normal reads
        std::vector<ReadSketch> window_reads;
        window_reads.reserve(win->priority_reads.size() + win->normal_reads.size());
        window_reads.insert(window_reads.end(),
            win->priority_reads.begin(), win->priority_reads.end());
        window_reads.insert(window_reads.end(),
            win->normal_reads.begin(), win->normal_reads.end());

        if (window_reads.size() < 2) continue;

        // Build components
        auto components = component_builder_.build(window_reads, win->chrom_tid);
        for (auto& comp : components) {
            comp.id = result.components.size();
            result.components.push_back(std::move(comp));
        }
    }

    result.built_components = result.components.size();
    result.elapsed_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start_time_).count();

    if (verbose_) {
        std::cout << "  Built " << result.built_components << " components in "
                  << result.elapsed_seconds << "s" << std::endl;
    }

    return result;
}

}  // namespace placer
