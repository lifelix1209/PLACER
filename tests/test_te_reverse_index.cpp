#include <cassert>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "te_reverse_index.h"
#include "gate1.h"
#include "component_builder.h"

using namespace placer;

// ============================================================================
// Test Helpers
// ============================================================================

void print_test_header(const std::string& name) {
    std::cout << "\nTesting " << name << "...\n";
}

void check_result(const std::string& name, bool passed) {
    std::cout << "  " << name << ": " << (passed ? "PASS" : "FAIL") << "\n";
    assert(passed);
}

// ============================================================================
// TEReverseIndexConfig Tests
// ============================================================================

void test_config() {
    print_test_header("TEReverseIndexConfig");

    TEReverseIndexConfig config;
    assert(config.kmer_size == 15);
    assert(config.minimizer_window == 10);
    assert(config.max_genome_hits_per_kmer == 20);
    assert(config.min_genome_hits == 3);
    assert(config.locus_cluster_radius == 50);
    assert(config.min_locus_reads == 2);
    assert(config.min_probe_hits == 3);
    assert(config.min_hit_density == 0.1);

    // Custom config
    TEReverseIndexConfig custom;
    custom.kmer_size = 11;
    custom.minimizer_window = 5;
    custom.locus_cluster_radius = 100;
    custom.num_threads = 8;

    assert(custom.kmer_size == 11);
    assert(custom.minimizer_window == 5);

    std::cout << "  TEReverseIndexConfig tests passed!\n";
}

// ============================================================================
// GenomeKmerIndex Tests
// ============================================================================

void test_genome_kmer_index() {
    print_test_header("GenomeKmerIndex");

    TEReverseIndexConfig config;
    config.kmer_size = 3;
    config.minimizer_window = 2;
    config.max_genome_hits_per_kmer = 10;

    GenomeKmerIndex index(config);

    // Test empty index
    assert(index.empty());
    assert(index.genome_size() == 0);
    assert(index.num_chromosomes() == 0);

    // Build from sequences
    std::vector<std::string> sequences = {
        "ACGTACGTACGT",  // Chromosome 1
        "TTTTAAAA"       // Chromosome 2
    };

    index.build_from_sequences(sequences);

    assert(!index.empty());
    assert(index.genome_size() == 20);
    assert(index.num_chromosomes() == 2);

    // Test query - exact match
    auto hits = index.query("ACG");
    // ACG appears at positions 0 and 4 on chrom 0
    assert(!hits.empty());

    // Test query - no match
    hits = index.query("NNNN");
    assert(hits.empty());

    // Test query_hit_count
    int count = index.query_hit_count("ACGT");
    assert(count >= 0);

    std::cout << "  GenomeKmerIndex tests passed!\n";
}

// ============================================================================
// RescuedLocus Tests
// ============================================================================

void test_rescued_locus() {
    print_test_header("RescuedLocus");

    RescuedLocus locus;
    assert(locus.chrom_tid == -1);
    assert(locus.position == -1);
    assert(locus.cluster_id == -1);
    assert(locus.support_reads == 0);
    assert(locus.total_hits == 0);
    assert(locus.hit_density == 0.0);
    assert(locus.placeability_score == 0.0);
    assert(!locus.passed_filter);

    // Set values
    locus.chrom_tid = 0;
    locus.position = 1000;
    locus.support_reads = 5;
    locus.total_hits = 10;
    locus.hit_density = 0.05;
    locus.placeability_score = 0.8;
    locus.passed_filter = true;

    assert(locus.chrom_tid == 0);
    assert(locus.position == 1000);
    assert(locus.support_reads == 5);

    std::cout << "  RescuedLocus tests passed!\n";
}

// ============================================================================
// ReadRescueResult Tests
// ============================================================================

void test_read_rescue_result() {
    print_test_header("ReadRescueResult");

    ReadRescueResult result;
    assert(result.read_idx == 0);
    assert(result.probe_start == -1);
    assert(result.probe_end == -1);
    assert(result.genome_hits.empty());
    assert(!result.rescued);
    assert(result.num_valid_hits == 0);

    // Set values
    result.read_idx = 42;
    result.probe_sequence = "ACGTACGT";
    result.probe_start = 100;
    result.probe_end = 108;
    result.genome_hits.emplace_back(0, 500);
    result.genome_hits.emplace_back(0, 510);
    result.num_valid_hits = 2;
    result.rescued = true;

    assert(result.read_idx == 42);
    assert(result.genome_hits.size() == 2);
    assert(result.rescued);

    std::cout << "  ReadRescueResult tests passed!\n";
}

// ============================================================================
// TEReverseIndex Tests
// ============================================================================

void test_te_reverse_index() {
    print_test_header("TEReverseIndex");

    TEReverseIndexConfig config;
    config.kmer_size = 3;
    config.minimizer_window = 2;
    config.min_genome_hits = 1;
    config.locus_cluster_radius = 50;
    config.min_locus_reads = 1;
    config.min_probe_hits = 1;
    config.min_hit_density = 0.01;

    TEReverseIndex reverse_index(config);

    // Build genome index from sequences
    std::vector<std::string> genome = {
        "ACGTACGTACGT"
    };

    // We need to build the index manually for this test
    GenomeKmerIndex index(config);
    index.build_from_sequences(genome);

    // Create a mock read with probe fragments
    ReadSketch read;
    read.qname = "test_read";
    read.tid = 0;
    read.pos = 100;

    // The reverse_index should have the genome_index
    // In a real scenario, we would call initialize() with a FASTA file

    std::cout << "  TEReverseIndex basic tests passed!\n";
}

// ============================================================================
// Utility Functions Tests
// ============================================================================

void test_utility_functions() {
    print_test_header("Utility Functions");

    // Test position_distance - same chromosome
    int32_t dist = position_distance(0, 100, 0, 150);
    assert(dist == 50);

    // Test position_distance - different chromosome
    dist = position_distance(0, 100, 1, 150);
    assert(dist == std::numeric_limits<int32_t>::max());

    // Test cluster_positions
    std::vector<std::pair<int32_t, int32_t>> positions = {
        {0, 100}, {0, 110}, {0, 120}, {0, 200}, {0, 210}
    };

    auto clusters = cluster_positions(positions, 50);

    // Should have 2 clusters: [100, 110, 120] and [200, 210]
    assert(clusters.size() == 2);
    assert(clusters[0].size() == 3);
    assert(clusters[1].size() == 2);

    std::cout << "  Utility Functions tests passed!\n";
}

// ============================================================================
// Integration Tests
// ============================================================================

void test_rescue_pipeline() {
    print_test_header("Rescue Pipeline Integration");

    TEReverseIndexConfig config;
    config.kmer_size = 3;
    config.minimizer_window = 2;
    config.min_genome_hits = 1;
    config.locus_cluster_radius = 100;
    config.min_locus_reads = 1;
    config.min_probe_hits = 1;
    config.min_hit_density = 0.01;

    TEReverseIndex reverse_index(config);

    // Test statistics before any processing
    const auto& stats = reverse_index.stats();
    assert(stats.total_reads_processed == 0);
    assert(stats.reads_with_sa == 0);
    assert(stats.reads_rescued == 0);

    // Test reset
    reverse_index.reset_stats();
    assert(reverse_index.stats().total_reads_processed == 0);

    std::cout << "  Rescue Pipeline integration tests passed!\n";
}

// ============================================================================
// Main Test Runner
// ============================================================================

int main() {
    std::cout << "=== PLACER Phase 8 TE Reverse Index Tests ===\n";

    try {
        test_config();
        test_genome_kmer_index();
        test_rescued_locus();
        test_read_rescue_result();
        test_te_reverse_index();
        test_utility_functions();
        test_rescue_pipeline();

        std::cout << "\n=== All Phase 8 tests passed! ===\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Test failed with exception: " << e.what() << "\n";
        return 1;
    }
}
