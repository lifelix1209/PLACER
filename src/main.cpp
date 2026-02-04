/**
 * PLACER - Phase 2: Stream Layer with Gate 1 Integration
 *
 * Demonstrates:
 * - Single-pass BAM streaming
 * - Window buffering with triggering
 * - Gate 1 TE-proxy filtering (optional)
 * - Task queue integration
 */

#include <iostream>
#include <chrono>
#include <memory>
#include "bam_reader.h"
#include "window_buffer.h"
#include "window_stats.h"
#include "trigger.h"
#include "task_queue.h"
#include "gate1_filter.h"

using namespace placer;

void print_usage(const char* prog) {
    std::cout << "Usage: " << prog << " <bam_file> [te_fasta]" << std::endl;
    std::cout << std::endl;
    std::cout << "PLACER Phase 2 - Stream Layer + Gate 1" << std::endl;
    std::cout << "  Streams BAM, buffers windows, triggers on anomalies" << std::endl;
    std::cout << "  Optionally filters reads through Gate 1 TE-proxy" << std::endl;
    std::cout << std::endl;
    std::cout << "Arguments:" << std::endl;
    std::cout << "  bam_file   Input BAM file (required)" << std::endl;
    std::cout << "  te_fasta   TE library FASTA (optional, enables Gate 1)" << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    std::string bam_path = argv[1];
    std::string te_fasta_path = (argc >= 3) ? argv[2] : "";

    std::cout << "=== PLACER Phase 2 Demo ===" << std::endl;
    std::cout << "Input BAM: " << bam_path << std::endl;
    if (!te_fasta_path.empty()) {
        std::cout << "TE Library: " << te_fasta_path << std::endl;
    } else {
        std::cout << "TE Library: (not specified, Gate 1 disabled)" << std::endl;
    }

    // Configure components
    WindowBuffer::Config wconfig;
    wconfig.window_size = 10000;
    wconfig.window_step = 5000;
    wconfig.max_priority_reads = 50;
    wconfig.max_normal_reads = 200;
    wconfig.min_clip_bp = 50;
    wconfig.min_sa_reads = 3;

    WindowBuffer window_buffer(wconfig);

    // Task queue for component building
    TaskQueue task_queue(2);

    // Gate 1 filter (optional)
    Gate1FilterConfig g1_config;
    g1_config.te_fasta_path = te_fasta_path;
    g1_config.enabled = !te_fasta_path.empty();
    g1_config.max_reads_per_window = 100;
    g1_config.sample_ratio = 1.0;

    Gate1Filter gate1_filter(g1_config);

    // Build TE index if enabled
    if (gate1_filter.is_enabled()) {
        std::cout << "Building TE index from " << te_fasta_path << "..." << std::endl;
        if (!gate1_filter.build_index()) {
            std::cerr << "Error: Failed to build TE index" << std::endl;
            return 1;
        }
        std::cout << "TE index built: " << gate1_filter.is_index_ready()
                  << " (size: " << gate1_filter.config().te_fasta_path << ")" << std::endl;
    }

    // Window processor with Gate 1 integration
    WindowProcessor::Config wp_config;
    wp_config.task_queue = &task_queue;
    wp_config.gate1_filter = gate1_filter.is_enabled() ? &gate1_filter : nullptr;
    WindowProcessor window_processor(wp_config);

    int64_t total_reads = 0;
    int64_t triggered_count = 0;

    auto start_time = std::chrono::high_resolution_clock::now();

    // Stream BAM
    BamReader reader(bam_path);
    if (!reader.is_valid()) {
        std::cerr << "Error: Cannot open BAM file: " << bam_path << std::endl;
        return 1;
    }

    std::cout << "Streaming BAM..." << std::endl;

    auto callback = [&](const ReadSketch& read) {
        total_reads++;
        window_buffer.add_read(read);

        // Progress report
        if (total_reads % 10000 == 0) {
            std::cout << "\r  Processed " << total_reads << " reads..."
                      << std::flush;
        }
    };

    int64_t processed = reader.stream(callback);

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    std::cout << "\r  Processed " << total_reads << " reads in "
              << duration.count() << "ms" << std::endl;

    // Summary
    std::cout << "\n=== Summary ===" << std::endl;
    std::cout << "Total reads processed: " << processed << std::endl;

    // Collect stats from all windows
    auto windows = window_buffer.list_windows();
    uint64_t total_clip_bp = 0;
    uint64_t total_sa_reads = 0;
    uint64_t total_ins_events = 0;

    for (auto id : windows) {
        const Window* win = window_buffer.get_window(id);
        if (win) {
            total_clip_bp += win->stats.clip_bp;
            total_sa_reads += win->stats.sa_reads;
            total_ins_events += win->stats.ins_events;
            if (win->triggered) triggered_count++;
        }
    }

    std::cout << "Active windows: " << windows.size() << std::endl;
    std::cout << "Triggered windows: " << triggered_count << std::endl;
    std::cout << "Total clip_bp: " << total_clip_bp << std::endl;
    std::cout << "Total sa_reads: " << total_sa_reads << std::endl;

    // Process triggered windows
    std::cout << "\nProcessing triggered windows..." << std::endl;
    auto sealed = window_buffer.seal_and_flush(1000000);

    // Use WindowProcessor to handle triggered windows with optional Gate 1 filtering
    window_processor.process_triggered_windows(sealed);

    // Close task queue and wait for completion
    task_queue.close();
    task_queue.wait();

    std::cout << "Memory reclamation: " << sealed.size() << " windows sealed" << std::endl;

    // Gate 1 statistics
    if (gate1_filter.is_enabled()) {
        std::cout << "\n=== Gate 1 Statistics ===" << std::endl;
        std::cout << "Reads processed: " << gate1_filter.get_total_reads_processed() << std::endl;
        std::cout << "Reads passed: " << gate1_filter.get_total_reads_passed() << std::endl;
        std::cout << "Reads filtered: " << gate1_filter.get_total_reads_filtered() << std::endl;

        double pass_rate = (gate1_filter.get_total_reads_processed() > 0)
            ? (100.0 * gate1_filter.get_total_reads_passed() / gate1_filter.get_total_reads_processed())
            : 0.0;
        std::cout << "Pass rate: " << pass_rate << "%" << std::endl;
    }

    // Window processor statistics
    std::cout << "\n=== Window Processor Statistics ===" << std::endl;
    std::cout << "Windows processed: " << window_processor.get_total_windows_processed() << std::endl;
    std::cout << "Reads submitted: " << window_processor.get_total_reads_submitted() << std::endl;
    std::cout << "Reads filtered: " << window_processor.get_total_reads_filtered() << std::endl;

    std::cout << "Throughput: " << (processed * 1000.0 / duration.count()) << " reads/sec" << std::endl;

    return 0;
}
