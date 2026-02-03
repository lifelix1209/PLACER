#include <iostream>
#include <cassert>
#include "bam_reader.h"
#include "window_buffer.h"
#include "window_stats.h"
#include "trigger.h"
#include "task_queue.h"

using namespace placer;

void test_bam_reader() {
    std::cout << "Testing BamReader..." << std::endl;

    BamReader reader("/mnt/home1/miska/hl725/projects/tldr_optimized/test/test.bam");
    assert(reader.is_valid());

    int64_t count = 0;
    auto callback = [&count](const ReadSketch& read) {
        count++;
        assert(!read.qname.empty());
        assert(read.tid >= 0);
        assert(read.pos >= 0);
        assert(read.end_pos > read.pos);
        assert(read.mapq >= 0);
        // Verify sequence extraction (required for Gate 1)
        assert(!read.sequence.empty() || read.pos == 0);
    };

    int64_t result = reader.stream(callback);
    std::cout << "  Processed " << result << " records" << std::endl;
    assert(result > 0);

    std::cout << "  BamReader tests passed!" << std::endl;
}

void test_window_buffer() {
    std::cout << "Testing WindowBuffer..." << std::endl;

    WindowBuffer::Config config;
    config.window_size = 10000;
    config.window_step = 5000;
    config.max_priority_reads = 10;
    config.max_normal_reads = 50;
    config.min_clip_bp = 50;
    config.min_sa_reads = 3;

    WindowBuffer buffer(config);

    // Test read spanning multiple windows
    ReadSketch read1;
    read1.qname = "read1";
    read1.tid = 0;
    read1.pos = 1000;
    read1.end_pos = 25000;  // Spans multiple 10kb windows
    read1.mapq = 60;
    read1.cigar_ops = {{'M', 24000}};

    buffer.add_read(read1);

    auto windows = buffer.list_windows();
    std::cout << "  Windows created: " << windows.size() << std::endl;
    assert(windows.size() >= 2);  // Should span at least 2 windows

    // Check that a window contains this read
    const Window* win = buffer.get_window(windows[0]);
    assert(win != nullptr);
    std::cout << "  Window stats: bases_covered=" << win->stats.bases_covered << std::endl;

    // Test seal_and_flush (memory reclamation)
    auto sealed = buffer.seal_and_flush(0, 50000);  // Safe pos past all windows
    std::cout << "  Sealed windows: " << sealed.size() << std::endl;

    std::cout << "  WindowBuffer tests passed!" << std::endl;
}

void test_window_stats() {
    std::cout << "Testing WindowStats..." << std::endl;

    // Test basic WindowStats struct
    WindowStats stats;
    assert(stats.clip_bp == 0);
    assert(stats.sa_reads == 0);
    assert(stats.ins_events == 0);
    assert(stats.q99_clip == 0.0);
    assert(stats.q99_sa == 0.0);
    assert(stats.q99_ins == 0.0);

    // Test QuantileTracker
    QuantileTracker tracker(0.99);
    assert(tracker.count() == 0);

    for (int i = 0; i < 100; i++) {
        tracker.add(i * 10);
    }

    double estimate = tracker.estimate();
    std::cout << "  Quantile estimate: " << estimate << std::endl;
    std::cout << "  Observations: " << tracker.count() << std::endl;

    assert(tracker.count() == 100);
    assert(estimate >= 0 && estimate <= 1000);
    assert(tracker.min() == 0);
    assert(tracker.max() == 990);
    assert(tracker.mean() == 495);

    // Test WindowStatsCollector
    WindowStatsCollector collector;
    collector.add_clip_bp("chr1", 100);
    collector.add_clip_bp("chr1", 200);
    collector.add_clip_bp("chr1", 300);

    double q99 = collector.get_q99_clip("chr1");
    std::cout << "  Q99 clip for chr1: " << q99 << std::endl;
    assert(q99 > 100);

    std::cout << "  WindowStats tests passed!" << std::endl;
}

void test_trigger() {
    std::cout << "Testing Trigger..." << std::endl;

    TriggerConfig config;
    Trigger trigger(config);

    WindowStats stats;
    stats.clip_bp = 500;
    stats.sa_reads = 5;
    stats.ins_events = 3;
    stats.q99_clip = 100.0;   // Set Q99 thresholds
    stats.q99_sa = 2.0;
    stats.q99_ins = 1.0;

    TriggerResult result = trigger.evaluate("test_window", stats);
    std::cout << "  Trigger score: " << result.score << std::endl;
    std::cout << "  Triggered: " << (result.triggered ? "yes" : "no") << std::endl;
    std::cout << "  Clip trig: " << (result.clip_trig ? "yes" : "no") << std::endl;

    assert(result.score > 0);
    assert(result.clip_trig);  // 500 > (100 * 2.0, 50) = 200, so should trigger

    // Test baseline protection (Q99 = 0 case)
    WindowStats zero_stats;
    zero_stats.clip_bp = 10;   // 10bp clip
    zero_stats.sa_reads = 0;
    zero_stats.ins_events = 0;
    zero_stats.q99_clip = 0.0;  // Cold start case
    zero_stats.q99_sa = 0.0;
    zero_stats.q99_ins = 0.0;

    TriggerResult zero_result = trigger.evaluate("zero_window", zero_stats);
    // Should NOT trigger because min_baseline_clip = 50, and 10 < 50
    assert(!zero_result.triggered);

    std::cout << "  Trigger tests passed!" << std::endl;
}

void test_task_queue() {
    std::cout << "Testing TaskQueue..." << std::endl;

    TaskQueue queue(2);

    for (int i = 0; i < 10; i++) {
        auto task = TaskFactory::create_component_build("window_" + std::to_string(i));
        queue.submit(std::move(task));
    }

    queue.wait();
    std::cout << "  All tasks processed" << std::endl;
    assert(queue.size() == 0);

    queue.close();
    assert(!queue.submit(TaskFactory::create_local_align("test")));

    std::cout << "  TaskQueue tests passed!" << std::endl;
}

void test_integration() {
    std::cout << "Testing integration (full pipeline)..." << std::endl;

    WindowBuffer::Config wconfig;
    wconfig.window_size = 10000;
    wconfig.max_priority_reads = 50;
    wconfig.max_normal_reads = 200;

    WindowBuffer window_buffer(wconfig);

    BamReader reader("/mnt/home1/miska/hl725/projects/tldr_optimized/test/test.bam");
    assert(reader.is_valid());

    int64_t count = 0;
    int64_t triggered_count = 0;

    auto callback = [&window_buffer, &triggered_count, &count](const ReadSketch& read) {
        count++;
        window_buffer.add_read(read);

        // Check if this read caused any window to trigger
        auto windows = window_buffer.list_windows();
        for (auto id : windows) {
            const Window* win = window_buffer.get_window(id);
            if (win && win->triggered) {
                triggered_count++;
            }
        }
    };

    reader.stream(callback);

    std::cout << "  Total reads: " << count << std::endl;

    // Test memory reclamation
    auto sealed = window_buffer.seal_and_flush(0, 100000);
    std::cout << "  Triggered windows sealed: " << sealed.size() << std::endl;

    std::cout << "  Integration tests passed!" << std::endl;
}

int main() {
    std::cout << "=== PLACER Phase 1 Tests ===" << std::endl;

    try {
        test_window_stats();
        test_bam_reader();
        test_window_buffer();
        test_trigger();
        test_task_queue();
        test_integration();

        std::cout << "\n=== All tests passed! ===" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Test failed: " << e.what() << std::endl;
        return 1;
    }
}
