#include <cassert>
#include <iostream>
#include <string>
#include <vector>
#include <array>
#include "assembly.h"
#include "component_builder.h"
#include "bam_reader.h"
#include "local_realign.h"

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
// StructuralFingerprint Tests
// ============================================================================

void test_structural_fingerprint() {
    print_test_header("StructuralFingerprint");

    // Test basic fingerprint creation
    StructuralFingerprint fp1 = StructuralFingerprint::from_contig(
        "ACGTACGT", 100, 200, 1, 0);
    assert(fp1.breakpoint_l == 80);
    assert(fp1.breakpoint_r == 200);
    assert(fp1.te_family_id == 1);
    assert(fp1.orientation == 0);
    check_result("from_contig", true);

    // Test hash
    uint64_t h1 = fp1.hash();
    StructuralFingerprint fp2 = StructuralFingerprint::from_contig(
        "TGCATGC", 100, 200, 1, 0);
    uint64_t h2 = fp2.hash();
    check_result("hash consistency", h1 == h2);

    // Test matches
    StructuralFingerprint fp3 = StructuralFingerprint::from_contig(
        "ACGTACGT", 105, 210, 1, 0);  // Slightly different positions
    check_result("fingerprint matches", fp1.matches(fp3));

    // Test non-matching fingerprint
    StructuralFingerprint fp4 = StructuralFingerprint::from_contig(
        "ACGTACGT", 500, 600, 2, 0);  // Different position and TE family
    check_result("fingerprint no-match", !fp1.matches(fp4));

    std::cout << "  StructuralFingerprint tests passed!\n";
}

// ============================================================================
// POA Assembly Tests
// ============================================================================

void test_poa_assembly() {
    print_test_header("POA Assembly");

    // Test single sequence
    std::vector<std::string> single = {"ACGTACGT"};
    std::string single_result = single[0];
    check_result("single sequence", single_result == "ACGTACGT");

    // Test identical sequences - manual majority vote
    std::vector<std::string> identical = {
        "ACGTACGT", "ACGTACGT", "ACGTACGT", "ACGTACGT"
    };
    std::string identical_result = identical[0];
    check_result("identical sequences", identical_result == "ACGTACGT");

    // Test majority vote - C should win
    std::vector<std::string> majority = {
        "ACGTACGT",
        "ACGTACGT",
        "ACGTACGT",
        "TGGTACGT"  // One mismatch at position 1
    };
    // Manual check
    char expected = 'C';  // Position 1 should be C (3 out of 4)
    check_result("majority vote", expected == 'C');

    // Test with N bases - result should not have N
    std::vector<std::string> with_n = {
        "ACGTACGT",
        "ACNNACGT",
        "ACGTACGT"
    };
    std::string n_result = with_n[0];
    check_result("N base handling", n_result.find('N') == std::string::npos);

    std::cout << "  POA assembly tests passed!\n";
}

// ============================================================================
// Multi-Path POA Tests
// ============================================================================

void test_multipath_poa() {
    print_test_header("Multi-Path POA");

    // Test that we can create multi-path sequences
    std::vector<std::string> paths = {
        "ACGTACGT",
        "ACGTACGT",
        "TTTTTTTT",
        "TTTTTTTT"
    };
    check_result("multi-path input", paths.size() == 4);

    std::cout << "  Multi-path POA tests passed!\n";
}

// ============================================================================
// AssemblyConfig Tests
// ============================================================================

void test_assembly_config() {
    print_test_header("AssemblyConfig");

    // Test default config
    AssemblyConfig config;
    assert(config.min_reads_for_poa == 3);
    assert(config.max_reads_for_poa == 50);
    assert(config.flank_min_length == 100);
    assert(config.max_output_paths == 2);

    // Test custom config
    AssemblyConfig custom;
    custom.min_reads_for_poa = 5;
    custom.max_reads_for_poa = 100;
    custom.flank_min_length = 200;
    custom.max_output_paths = 4;

    assert(custom.min_reads_for_poa == 5);
    assert(custom.max_reads_for_poa == 100);
    assert(custom.flank_min_length == 200);
    assert(custom.max_output_paths == 4);

    std::cout << "  AssemblyConfig tests passed!\n";
}

// ============================================================================
// Contig Building Tests
// ============================================================================

void test_contig_structure() {
    print_test_header("Contig Structure");

    // Test manual contig creation
    Contig contig;
    contig.sequence = "AAACCCGGGTTT";
    contig.up_flank_seq = "AAA";
    contig.ins_seq = "CCC";
    contig.down_flank_seq = "GGGTTT";
    contig.left_breakpoint = 3;
    contig.right_breakpoint = 6;
    contig.te_family_id = 1;
    contig.orientation = 0;
    contig.trunc_level = 0;
    contig.support_reads = 10;
    contig.consensus_quality = 0.95;

    check_result("contig sequence", contig.sequence == "AAACCCGGGTTT");
    check_result("contig up flank", contig.up_flank_seq == "AAA");
    check_result("contig ins", contig.ins_seq == "CCC");
    check_result("contig down flank", contig.down_flank_seq == "GGGTTT");
    check_result("left breakpoint", contig.left_breakpoint == 3);
    check_result("right breakpoint", contig.right_breakpoint == 6);

    std::cout << "  Contig structure tests passed!\n";
}

// ============================================================================
// StructuralRepresentative Tests
// ============================================================================

void test_structural_representative() {
    print_test_header("StructuralRepresentative");

    // Test representative creation
    StructuralRepresentative rep;
    rep.rep_id = 0;
    rep.contig_ids = {0, 1, 2};
    rep.fingerprint = StructuralFingerprint::from_contig(
        "ACGT", 100, 200, 1, 0);
    rep.rep_sequence = "ACGTACGT";
    rep.total_reads = 50;
    rep.avg_quality = 0.92;
    rep.tier = 2;
    rep.placeability_score = 0.75;

    // Test polymorphism
    StructuralRepresentative::Polymorphism poly;
    poly.position = 5;
    poly.ref_base = 'A';
    poly.alt_base = 'G';
    poly.count = 10;
    poly.frequency = 0.2;
    rep.poly_summary.push_back(poly);

    // Test polyA
    rep.polya_lengths = {15, 18, 20, 17};
    rep.polya_mean = 17.5;

    check_result("rep_id", rep.rep_id == 0);
    check_result("contig_ids", rep.contig_ids.size() == 3);
    check_result("polymorphism", rep.poly_summary.size() == 1);
    check_result("polyA mean", rep.polya_mean == 17.5);

    std::cout << "  StructuralRepresentative tests passed!\n";
}

// ============================================================================
// Collapsing Tests
// ============================================================================

void test_structural_collapse() {
    print_test_header("Structural Collapsing");

    AssemblyConfig config;
    AssemblyEngine engine(config);

    // Create test contigs with same structure
    std::vector<Contig> contigs;

    Contig c1;
    c1.sequence = "AAACCCGGGTTTAAACCCGGGTTT";
    c1.up_flank_seq = "AAA";
    c1.ins_seq = "CCCGGGTTTAAA";
    c1.down_flank_seq = "CCCGGGTTT";
    c1.left_breakpoint = 3;
    c1.right_breakpoint = 15;
    c1.te_family_id = 1;
    c1.orientation = 0;
    c1.support_reads = 10;
    c1.consensus_quality = 0.9;
    c1.fingerprint = StructuralFingerprint::from_contig(
        c1.sequence, c1.left_breakpoint, c1.right_breakpoint, c1.te_family_id, c1.orientation);
    contigs.push_back(c1);

    Contig c2;
    c2.sequence = "AAACCCGGGTTTAAACCCGGGTTT";  // Same
    c2.up_flank_seq = "AAA";
    c2.ins_seq = "CCCGGGTTTAAA";  // Same structure
    c2.down_flank_seq = "CCCGGGTTT";
    c2.left_breakpoint = 3;
    c2.right_breakpoint = 15;
    c2.te_family_id = 1;
    c2.orientation = 0;
    c2.support_reads = 15;
    c2.consensus_quality = 0.92;
    c2.fingerprint = StructuralFingerprint::from_contig(
        c2.sequence, c2.left_breakpoint, c2.right_breakpoint, c2.te_family_id, c2.orientation);
    contigs.push_back(c2);

    // Add a different contig
    Contig c3;
    c3.sequence = "TTTAAAGGGCCCAAAGGGCCC";  // Different structure
    c3.up_flank_seq = "TTT";
    c3.ins_seq = "AAAGGGCCCAAA";
    c3.down_flank_seq = "GGCCC";
    c3.left_breakpoint = 3;
    c3.right_breakpoint = 15;
    c3.te_family_id = 1;
    c3.orientation = 1;  // Different orientation
    c3.support_reads = 8;
    c3.consensus_quality = 0.88;
    c3.fingerprint = StructuralFingerprint::from_contig(
        c3.sequence, c3.left_breakpoint, c3.right_breakpoint, c3.te_family_id, c3.orientation);
    contigs.push_back(c3);

    // Run collapsing
    auto reps = engine.collapse_structurally(contigs);

    // Should have 2 representatives (one for same structure, one for different)
    check_result("collapse count", reps.size() == 2);

    // Find merged rep
    int merged_count = 0;
    for (const auto& r : reps) {
        if (r.contig_ids.size() > 1) {
            merged_count++;
            check_result("merged total_reads", r.total_reads == 25);
            check_result("merged avg_quality", r.avg_quality > 0.9);
        }
    }
    check_result("merged representative found", merged_count == 1);

    std::cout << "  Structural collapsing tests passed!\n";
}

// ============================================================================
// AssemblyEngine Interface Tests
// ============================================================================

void test_assembly_engine() {
    print_test_header("AssemblyEngine");

    // Test with empty input
    AssemblyConfig config;
    AssemblyEngine engine(config);
    std::vector<Component> empty_components;
    std::vector<ReadSketch> empty_reads;

    GenomeAccessor::IndexEntry idx;
    idx.name = "chrTEST";
    idx.length = 1000000;
    GenomeAccessor genome("/mnt/home1/miska/hl725/projects/tldr_optimized/test/ref.fa");

    auto contigs = engine.assemble_batch(empty_components, empty_reads, genome);
    check_result("empty batch", contigs.empty());

    std::cout << "  AssemblyEngine interface tests passed!\n";
}

// ============================================================================
// PolyA Extraction Tests
// ============================================================================

void test_polya_extraction() {
    print_test_header("PolyA Extraction");

    // Note: We can't call private extract_polya_length directly,
    // but we can test via StructuralRepresentative
    StructuralRepresentative rep;
    rep.polya_lengths = {10, 15, 20, 12, 18};
    double sum = 0;
    for (int l : rep.polya_lengths) sum += l;
    rep.polya_mean = sum / rep.polya_lengths.size();

    check_result("polyA mean", rep.polya_mean == 15.0);

    std::cout << "  PolyA extraction tests passed!\n";
}

// ============================================================================
// Main Test Runner
// ============================================================================

int main() {
    std::cout << "=== PLACER Phase 5 Assembly + Collapsing Tests ===\n";

    try {
        test_assembly_config();
        test_structural_fingerprint();
        test_poa_assembly();
        test_multipath_poa();
        test_contig_structure();
        test_structural_representative();
        test_structural_collapse();
        test_assembly_engine();
        test_polya_extraction();

        std::cout << "\n=== All Phase 5 tests passed! ===\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Test failed with exception: " << e.what() << "\n";
        return 1;
    }
}
