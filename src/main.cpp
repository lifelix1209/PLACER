/**
 * PLACER - Phase 1: Stream Layer Demo
 *
 * Demonstrates single-pass BAM streaming, window buffering, and triggering.
 */

#include <iostream>
#include <chrono>
#include "bam_reader.h"
#include "window_buffer.h"
#include "window_stats.h"

using namespace placer;

void print_usage(const char* prog) {
    std::cout << "Usage: " << prog << " <bam_file>" << std::endl;
    std::cout << std::endl;
    std::cout << "PLACER Phase 1 - Stream Layer Demo" << std::endl;
    std::cout << "  Streams BAM, buffers windows, triggers on anomalies" << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    std::string bam_path = argv[1];

    std::cout << "=== PLACER Phase 1 Demo ===" << std::endl;
    std::cout << "Input: " << bam_path << std::endl;

    // Configure components
    WindowBuffer::Config wconfig;
    wconfig.window_size = 10000;
    wconfig.window_step = 5000;
    wconfig.max_priority_reads = 50;
    wconfig.max_normal_reads = 200;
    wconfig.min_clip_bp = 50;
    wconfig.min_sa_reads = 3;

    WindowBuffer window_buffer(wconfig);
    WindowStatsCollector stats_collector;

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
    std::cout << "Throughput: " << (processed * 1000.0 / duration.count()) << " reads/sec" << std::endl;

    // Memory reclamation test
    auto sealed = window_buffer.seal_and_flush(1000000);
    std::cout << "\nMemory reclamation: " << sealed.size() << " windows sealed" << std::endl;

    return 0;
}
