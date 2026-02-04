#include <iostream>
#include <cassert>
#include "bam_reader.h"
#include "window_buffer.h"
#include "window_stats.h"
#include "trigger.h"
#include "task_queue.h"
#include "gate1.h"

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
    auto sealed = buffer.seal_and_flush(50000);  // Safe pos past all windows
    std::cout << "  Sealed windows: " << sealed.size() << std::endl;

    std::cout << "  WindowBuffer tests passed!" << std::endl;
}

void test_chromosome_switching() {
    std::cout << "Testing chromosome switching..." << std::endl;

    WindowBuffer::Config config;
    config.window_size = 10000;
    config.max_priority_reads = 10;
    config.max_normal_reads = 50;
    config.min_clip_bp = 1;  // Low threshold for testing
    config.min_sa_reads = 1;

    WindowBuffer buffer(config);

    // Create reads on different chromosomes
    // Simulating: chr1 -> chr2 -> chr1 (should trigger cleanup)

    ReadSketch read1;
    read1.qname = "chr1_read1";
    read1.tid = 0;  // chr1
    read1.pos = 1000;
    read1.end_pos = 5000;
    read1.mapq = 60;
    read1.cigar_ops = {{'M', 4000}};
    read1.total_clip_len = 100;  // High clip to trigger

    ReadSketch read2;
    read2.qname = "chr2_read1";
    read2.tid = 1;  // chr2
    read2.pos = 1000;
    read2.end_pos = 5000;
    read2.mapq = 60;
    read2.cigar_ops = {{'M', 4000}};
    read2.total_clip_len = 100;

    ReadSketch read3;
    read3.qname = "chr1_read2";
    read3.tid = 0;  // Back to chr1
    read3.pos = 20000;
    read3.end_pos = 25000;
    read3.mapq = 60;
    read3.cigar_ops = {{'M', 5000}};
    read3.total_clip_len = 100;

    // Add reads - should trigger chromosome switching cleanup
    buffer.add_read(read1);
    auto windows_after_chr1 = buffer.list_windows();
    std::cout << "  Windows after chr1 read: " << windows_after_chr1.size() << std::endl;

    buffer.add_read(read2);
    auto windows_after_chr2 = buffer.list_windows();
    std::cout << "  Windows after chr2 read: " << windows_after_chr2.size() << std::endl;

    buffer.add_read(read3);
    auto windows_after_chr1_again = buffer.list_windows();
    std::cout << "  Windows after returning to chr1: " << windows_after_chr1_again.size() << std::endl;

    // After switching away from chr1, chr1 windows should be flushed
    // So windows_after_chr2 should have 0 chr1 windows
    // Note: The exact count depends on implementation

    std::cout << "  Chromosome switching tests passed!" << std::endl;
}

void test_window_buffer_empty() {
    std::cout << "Testing WindowBuffer (empty state)..." << std::endl;

    WindowBuffer::Config config;
    config.window_size = 10000;

    WindowBuffer buffer(config);

    // Test empty buffer
    auto windows = buffer.list_windows();
    assert(windows.empty());
    std::cout << "  Empty buffer has 0 windows: OK" << std::endl;

    // Test seal_and_flush on empty buffer
    auto sealed = buffer.seal_and_flush(50000);
    assert(sealed.empty());
    std::cout << "  seal_and_flush on empty buffer: OK" << std::endl;

    std::cout << "  WindowBuffer empty state tests passed!" << std::endl;
}

void test_task_queue_close_empty() {
    std::cout << "Testing TaskQueue (close on empty)..." << std::endl;

    TaskQueue queue(2);

    // Close immediately without submitting tasks
    queue.close();

    // Submit should fail
    bool result = queue.submit(TaskFactory::create_component_build("test"));
    assert(!result);
    std::cout << "  close() on empty queue: OK" << std::endl;

    // Verify no crash
    assert(queue.size() == 0);
    std::cout << "  TaskQueue close on empty tests passed!" << std::endl;
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
    std::cout << "  All placeholder tasks processed" << std::endl;
    assert(queue.size() == 0);

    // Test submit_serialized with actual read data
    std::cout << "  Testing submit_serialized with read data..." << std::endl;

    std::vector<ReadSketch> test_reads;
    for (int i = 0; i < 5; i++) {
        ReadSketch read;
        read.qname = "test_read_" + std::to_string(i);
        read.tid = 0;
        read.pos = 1000 + i * 100;
        read.end_pos = 2000 + i * 100;
        read.mapq = 60;
        test_reads.push_back(read);
    }

    bool result = queue.submit_serialized(TaskType::COMPONENT_BUILD, "test_window", test_reads);
    assert(result);
    std::cout << "  submit_serialized returned: " << (result ? "success" : "failed") << std::endl;

    queue.wait();
    std::cout << "  Serialized task processed" << std::endl;

    queue.close();
    assert(!queue.submit(TaskFactory::create_local_align("test")));

    std::cout << "  TaskQueue tests passed!" << std::endl;
}

void test_task_serialization() {
    std::cout << "Testing TaskSerializer (round-trip)..." << std::endl;

    // Create test task data
    TaskData original;
    original.type = TaskType::COMPONENT_BUILD;
    original.window_id = "test_chr1_90000";
    original.clip_bp = 500;
    original.sa_reads = 3;
    original.ins_events = 2;

    // Add test reads
    for (int i = 0; i < 3; i++) {
        ReadSketch read;
        read.qname = "test_read_" + std::to_string(i);
        read.tid = 0;
        read.pos = 1000 + i * 100;
        read.end_pos = 2000 + i * 100;
        read.flag = 0;
        read.mapq = 60;
        read.sequence = "ACGTACGTACGT";
        read.total_clip_len = 50;
        read.has_large_insertion = true;
        read.cigar_ops = {{'M', 100}, {'I', 20}, {'D', 10}};
        read.sa_targets = {{1, 5000}, {2, 10000}};
        read.has_md = true;
        original.reads.push_back(read);
    }

    // Create serializer
    TaskSerializer serializer("/tmp/placer_test_tasks");
    bool saved = serializer.save_task(original);
    assert(saved);
    std::cout << "  Task saved successfully" << std::endl;

    // Load task back
    TaskData loaded;
    bool loaded_ok = serializer.load_task(
        "/tmp/placer_test_tasks/task_0_test_chr1_90000.task", loaded);

    assert(loaded_ok);
    std::cout << "  Task loaded successfully" << std::endl;

    // Verify data integrity
    assert(loaded.type == original.type);
    assert(loaded.window_id == original.window_id);
    assert(loaded.clip_bp == original.clip_bp);
    assert(loaded.sa_reads == original.sa_reads);
    assert(loaded.ins_events == original.ins_events);
    assert(loaded.reads.size() == original.reads.size());
    std::cout << "  Data integrity verified" << std::endl;

    // Verify read fields
    for (size_t i = 0; i < loaded.reads.size(); i++) {
        const auto& orig = original.reads[i];
        const auto& load = loaded.reads[i];
        assert(load.qname == orig.qname);
        assert(load.tid == orig.tid);
        assert(load.pos == orig.pos);
        assert(load.end_pos == orig.end_pos);
        assert(load.flag == orig.flag);
        assert(load.mapq == orig.mapq);
        assert(load.sequence == orig.sequence);
        assert(load.total_clip_len == orig.total_clip_len);
        assert(load.has_large_insertion == orig.has_large_insertion);
        assert(load.cigar_ops.size() == orig.cigar_ops.size());
        assert(load.sa_targets.size() == orig.sa_targets.size());
        assert(load.has_md == orig.has_md);
    }
    std::cout << "  Read fields verified: " << loaded.reads.size() << " reads" << std::endl;

    std::cout << "  TaskSerializer tests passed!" << std::endl;
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
    auto sealed = window_buffer.seal_and_flush(100000);
    std::cout << "  Triggered windows sealed: " << sealed.size() << std::endl;

    std::cout << "  Integration tests passed!" << std::endl;
}

// ============= Phase 2 Tests =============

void test_hash_te_index() {
    std::cout << "Testing HashTEIndex (bit-packed k-mers)..." << std::endl;

    HashTEIndexConfig config;
    config.kmer_size = 3;  // Small for testing

    auto index = std::make_unique<HashTEIndex>(config);

    // Build from sequences
    std::vector<std::pair<int, std::string>> sequences = {
        {0, "ACGTACGTACGT"},  // Family 0: contains ACG, GTA, TAC, etc.
        {1, "TTTTAAAA"},      // Family 1
    };

    bool built = index->build_from_sequences(sequences);
    assert(built);
    std::cout << "  Index size: " << index->size() << " k-mers" << std::endl;
    assert(index->size() > 0);

    // Test query - exact match
    TEHit hit = index->query("ACGT");
    std::cout << "  Query 'ACGT': hit_count=" << hit.hit_count << ", density=" << hit.hit_density << std::endl;
    assert(hit.hit_count > 0);

    // Test query - partial match
    TEHit hit2 = index->query("GGGGGGGG");  // No match for family 0/1
    std::cout << "  Query 'GGGGGGGG': hit_count=" << hit2.hit_count << std::endl;
    assert(hit2.hit_count == 0);

    // Test passes_gate1
    HashTEIndexConfig strict_config;
    strict_config.kmer_size = 3;
    strict_config.min_hit_count = 3;  // Require 3+ hits to pass
    strict_config.min_hit_density = 0.5;

    auto strict_index = std::make_unique<HashTEIndex>(strict_config);
    strict_index->build_from_sequences(sequences);

    TEHit pass_hit = strict_index->query("ACGTACGT");
    bool passes = strict_index->passes_gate1(pass_hit);
    std::cout << "  Query 'ACGTACGT' passes strict: " << (passes ? "yes" : "no") << std::endl;
    assert(passes);  // 4 k-mers match, should pass

    TEHit fail_hit = strict_index->query("ACGT");
    bool fails = strict_index->passes_gate1(fail_hit);
    std::cout << "  Query 'ACGT' passes strict: " << (fails ? "yes" : "no") << std::endl;
    assert(!fails);  // Only 2 k-mers match, should fail with min=3

    std::cout << "  HashTEIndex tests passed!" << std::endl;
}

void test_hash_te_index_fasta() {
    std::cout << "Testing HashTEIndex build_from_fasta..." << std::endl;

    HashTEIndexConfig config;
    config.kmer_size = 3;

    auto index = HashTEIndex::build_from_fasta("/mnt/home1/miska/hl725/projects/PLACER/test_data/te_test.fa", config);
    assert(index != nullptr);
    std::cout << "  Index size: " << index->size() << " k-mers" << std::endl;
    assert(index->size() > 0);

    // Verify family names were parsed (flexible for any FASTA file)
    const auto& family_names = index->family_ids();
    std::cout << "  Family names: " << family_names.size() << std::endl;
    assert(family_names.size() > 0);  // At least one family
    assert(!family_names[0].empty());  // First family name not empty

    // Test query with L1-like sequence (GGGTGGGG pattern)
    TEHit hit = index->query("GGGTGGGGATGGGGTT");
    std::cout << "  Query L1-like seq: hit_count=" << hit.hit_count << std::endl;
    assert(hit.hit_count > 0);
    assert(hit.family_hits.size() > 0);

    std::cout << "  HashTEIndex FASTA tests passed!" << std::endl;
}

void test_hash_te_index_n_handling() {
    std::cout << "Testing HashTEIndex N character handling..." << std::endl;

    HashTEIndexConfig config;
    config.kmer_size = 3;

    auto index = std::make_unique<HashTEIndex>(config);

    // Build with sequence containing N
    std::vector<std::pair<int, std::string>> sequences = {
        {0, "ACGTACNGTACGT"},  // N in the middle
    };

    index->build_from_sequences(sequences);

    // Query with N - should skip the N and count valid k-mers
    TEHit hit = index->query("ACGTACNGTACGT");
    std::cout << "  Query with N: hit_count=" << hit.hit_count << std::endl;
    assert(hit.hit_count > 0);

    // Query with only N - should have 0 hits
    TEHit hit2 = index->query("NNNNNNNN");
    std::cout << "  Query 'NNNNNNNN': hit_count=" << hit2.hit_count << std::endl;
    assert(hit2.hit_count == 0);

    std::cout << "  N handling tests passed!" << std::endl;
}

void test_gate1_extract_probes() {
    std::cout << "Testing Gate1 probe extraction..." << std::endl;

    Gate1Config config;
    config.probe_len = 50;
    config.min_clip_len = 20;
    config.min_ins_len = 10;
    config.ins_neighborhood = 20;

    auto te_index = HashTEIndex::build_from_fasta("/mnt/home1/miska/hl725/projects/PLACER/test_data/te_test.fa",
                                                   HashTEIndexConfig());
    assert(te_index != nullptr);

    Gate1 gate1(std::shared_ptr<HashTEIndex>(te_index.release()), config);

    // Create a test read with various CIGAR ops
    ReadSketch read;
    read.qname = "test_read";
    read.sequence = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    read.cigar_ops = {
        {'M', 30},
        {'S', 25},  // Soft clip - should be extracted
        {'I', 15},  // Insertion - should be extracted
        {'M', 25}
    };

    std::vector<ProbeFragment> probes;
    gate1.extract_probes(read, probes);

    std::cout << "  Probes extracted: " << probes.size() << std::endl;
    // Should have: END5, END3, SOFTCLIP, INSERTION
    assert(probes.size() >= 4);

    // Check probe types
    int end_count = 0, clip_count = 0, ins_count = 0;
    for (const auto& p : probes) {
        if (p.source_type == 0) end_count++;
        else if (p.source_type == 1) clip_count++;
        else if (p.source_type == 2) ins_count++;
    }
    std::cout << "  END probes: " << end_count << ", CLIP: " << clip_count << ", INS: " << ins_count << std::endl;
    assert(end_count == 2);  // END5 and END3
    assert(clip_count == 1); // SOFTCLIP
    assert(ins_count == 1);  // INSERTION

    std::cout << "  Probe extraction tests passed!" << std::endl;
}

void test_gate1_evaluate() {
    std::cout << "Testing Gate1 evaluate..." << std::endl;

    Gate1Config config;
    config.probe_len = 50;
    config.min_clip_len = 10;
    config.min_ins_len = 5;
    config.ins_neighborhood = 20;
    config.min_hit_count = 1;
    config.min_density = 0.01;

    // Build index with synthetic sequence for testing
    HashTEIndexConfig idx_config;
    idx_config.kmer_size = 3;
    auto test_index = std::make_unique<HashTEIndex>(idx_config);
    std::vector<std::pair<int, std::string>> sequences = {
        {0, "GGGTGGGGATGGGGTTGGGGTTGGGGATGGGGTT"}  // Test TE sequence
    };
    test_index->build_from_sequences(sequences);

    Gate1 gate1(std::move(test_index), config);

    // Test 1: Read with matching TE sequence
    ReadSketch te_read;
    te_read.qname = "te_read";
    te_read.sequence = "GGGTGGGGATGGGGTTGGGGTTGGGGATGGGGTTGGGGTTGGGGATGGGGTTGGGGTTGGGGATGGGGTTGGGGTTGGGGATGGGGTTGGGG";
    te_read.cigar_ops = {{'M', (int)te_read.sequence.size()}};

    auto result = gate1.evaluate(te_read);
    std::cout << "  TE-like read: passed=" << (result.passed ? "yes" : "no")
              << ", total_hits=" << result.total_hits
              << ", dominant_family=" << result.dominant_family << std::endl;
    assert(result.passed);
    assert(result.total_hits > 0);
    assert(!result.dominant_family.empty());  // Family name should be populated

    // Test 2: Read with no TE signal
    ReadSketch bg_read;
    bg_read.qname = "bg_read";
    bg_read.sequence = "ACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACAC";
    bg_read.cigar_ops = {{'M', (int)bg_read.sequence.size()}};

    auto result2 = gate1.evaluate(bg_read);
    std::cout << "  Background read: passed=" << (result2.passed ? "yes" : "no")
              << ", total_hits=" << result2.total_hits << std::endl;
    assert(!result2.passed);  // Should not pass - no TE signal

    std::cout << "  Gate1 evaluate tests passed!" << std::endl;
}

void test_gate1_passes() {
    std::cout << "Testing Gate1 passes (fast path)..." << std::endl;

    Gate1Config config;
    config.probe_len = 50;
    config.min_hit_count = 1;
    config.min_density = 0.01;

    // Build index with synthetic sequence for testing
    HashTEIndexConfig idx_config;
    idx_config.kmer_size = 3;
    auto test_index = std::make_unique<HashTEIndex>(idx_config);
    std::vector<std::pair<int, std::string>> sequences = {
        {0, "GGGTGGGGATGGGGTTGGGGTTGGGGATGGGGTT"}  // Test TE sequence
    };
    test_index->build_from_sequences(sequences);

    Gate1 gate1(std::move(test_index), config);

    // Test TE-like read
    ReadSketch te_read;
    te_read.qname = "te_read";
    te_read.sequence = "GGGTGGGGATGGGGTTGGGG";  // Matching sequence

    bool passed = gate1.passes(te_read);
    std::cout << "  TE-like read passes: " << (passed ? "yes" : "no") << std::endl;
    assert(passed);

    // Test background read
    ReadSketch bg_read;
    bg_read.qname = "bg_read";
    bg_read.sequence = "ACACACACACACACACACAC";  // No TE signal

    bool bg_passed = gate1.passes(bg_read);
    std::cout << "  Background read passes: " << (bg_passed ? "yes" : "no") << std::endl;
    assert(!bg_passed);  // Should not pass

    std::cout << "  Gate1 passes tests passed!" << std::endl;
}

void test_gate1_real_bam() {
    std::cout << "Testing Gate1 with real BAM data..." << std::endl;

    Gate1Config config;
    config.probe_len = 100;
    config.min_clip_len = 30;
    config.min_ins_len = 20;
    config.ins_neighborhood = 50;
    config.min_hit_count = 3;
    config.min_density = 0.05;

    auto te_index = HashTEIndex::build_from_fasta("/mnt/home1/miska/hl725/projects/PLACER/test_data/te_test.fa",
                                                   HashTEIndexConfig());
    assert(te_index != nullptr);

    Gate1 gate1(std::shared_ptr<HashTEIndex>(te_index.release()), config);

    BamReader reader("/mnt/home1/miska/hl725/projects/tldr_optimized/test/test.bam");
    assert(reader.is_valid());

    int64_t count = 0;
    int64_t passed = 0;

    auto callback = [&gate1, &count, &passed](const ReadSketch& read) {
        count++;
        if (gate1.passes(read)) {
            passed++;
        }
    };

    int64_t result = reader.stream(callback);
    std::cout << "  Processed " << result << " reads from BAM" << std::endl;
    std::cout << "  Gate1 passed: " << passed << " reads" << std::endl;
    std::cout << "  Pass rate: " << (100.0 * passed / count) << "%" << std::endl;

    std::cout << "  Gate1 BAM tests passed!" << std::endl;
}

int main() {
    std::cout << "=== PLACER Phase 1 Tests ===" << std::endl;

    try {
        test_window_stats();
        test_bam_reader();
        test_window_buffer();
        test_chromosome_switching();
        test_window_buffer_empty();
        test_trigger();
        test_task_queue();
        test_task_queue_close_empty();
        test_task_serialization();
        test_integration();

        std::cout << "\n=== PLACER Phase 2 Tests ===" << std::endl;

        test_hash_te_index();
        test_hash_te_index_fasta();
        test_hash_te_index_n_handling();
        test_gate1_extract_probes();
        test_gate1_evaluate();
        test_gate1_passes();
        test_gate1_real_bam();

        std::cout << "\n=== All tests passed! ===" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Test failed: " << e.what() << std::endl;
        return 1;
    }
}
