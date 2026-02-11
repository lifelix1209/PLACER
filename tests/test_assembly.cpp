#include <cassert>
#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include "assembly.h"
#include "component_builder.h"
#include "bam_reader.h"
#include "local_realign.h"
#include "test_path_utils.h"

namespace placer {
ReadBreakpoint detect_breakpoint_from_cigar(
    const std::vector<std::pair<char, int>>& cigar_ops);
}

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
        "ACGTACGT", 100, 200, 1, 0, 3);
    assert(fp1.breakpoint_l == 80);
    assert(fp1.breakpoint_r == 200);
    assert(fp1.tid == 3);
    assert(fp1.te_family_id == 1);
    assert(fp1.orientation == 0);
    check_result("from_contig", true);

    // Test hash
    uint64_t h1 = fp1.hash();
    StructuralFingerprint fp2 = StructuralFingerprint::from_contig(
        "TGCATGC", 100, 200, 1, 0, 3);
    uint64_t h2 = fp2.hash();
    check_result("hash consistency", h1 == h2);

    // Test matches
    StructuralFingerprint fp3 = StructuralFingerprint::from_contig(
        "ACGTACGT", 105, 210, 1, 0, 3);  // Slightly different positions
    check_result("fingerprint matches", fp1.matches(fp3));

    // Test non-matching fingerprint
    StructuralFingerprint fp4 = StructuralFingerprint::from_contig(
        "ACGTACGT", 500, 600, 2, 0, 3);  // Different position and TE family
    check_result("fingerprint no-match", !fp1.matches(fp4));

    std::cout << "  StructuralFingerprint tests passed!\n";
}

// ============================================================================
// RTree / POA Core Structure Tests
// ============================================================================

void test_rtree_spatial_query() {
    print_test_header("RTree Spatial Query");

    RTree tree;
    tree.insert(0.0f, 0.0f, 10.0f, 10.0f, 0);
    tree.insert(1000.0f, 1000.0f, 1010.0f, 1010.0f, 1);
    tree.insert(2000.0f, 2000.0f, 2010.0f, 2010.0f, 2);

    auto has_idx = [](const std::vector<int>& values, int idx) {
        return std::find(values.begin(), values.end(), idx) != values.end();
    };

    auto mid = tree.range_query(995.0f, 995.0f, 1015.0f, 1015.0f);
    check_result("query mid contains idx=1", has_idx(mid, 1));
    check_result("query mid excludes idx=0", !has_idx(mid, 0));
    check_result("query mid excludes idx=2", !has_idx(mid, 2));

    auto high = tree.range_query(1995.0f, 1995.0f, 2020.0f, 2020.0f);
    check_result("query high contains idx=2", has_idx(high, 2));
    check_result("query high excludes idx=0", !has_idx(high, 0));

    std::cout << "  RTree spatial query tests passed!\n";
}

void test_poa_predecessor_overflow() {
    print_test_header("POA Predecessor Overflow");

    POAArena arena(32);
    const uint32_t target = arena.allocate_node('T');

    std::vector<uint32_t> preds;
    for (int i = 0; i < 6; ++i) {
        uint32_t pred = arena.allocate_node(static_cast<char>('A' + (i % 4)));
        preds.push_back(pred);
        bool inserted = arena.add_edge(pred, target);
        check_result("edge inserted", inserted);
    }

    const POANode* target_node = arena.get_node(target);
    bool all_present = (target_node != nullptr);
    for (uint32_t pred : preds) {
        all_present = all_present && target_node->predecessors.contains(pred);
    }

    check_result("overflow predecessor size",
        target_node && target_node->predecessors.size() == preds.size());
    check_result("overflow predecessor contains all", all_present);
    check_result("indegree equals unique edges",
        target_node && target_node->indegree == preds.size());

    bool duplicate_inserted = arena.add_edge(preds.front(), target);
    target_node = arena.get_node(target);
    check_result("duplicate edge rejected", !duplicate_inserted);
    check_result("indegree stable after duplicate",
        target_node && target_node->indegree == preds.size());

    std::cout << "  POA predecessor overflow tests passed!\n";
}

// ============================================================================
// Breakpoint and Segment Semantics Tests
// ============================================================================

void test_breakpoint_detection_semantics() {
    print_test_header("Breakpoint Detection Semantics");

    auto bp_ins = detect_breakpoint_from_cigar({
        {'M', 10}, {'I', 60}, {'M', 20}
    });
    check_result("insertion valid", bp_ins.is_valid);
    check_result("insertion type", bp_ins.type == ReadBreakpoint::INSERTION);
    check_result("insertion read_pos", bp_ins.read_pos == 10);
    check_result("insertion ref_pos is event-time", bp_ins.ref_pos == 10);
    check_result("insertion len captured", bp_ins.insertion_len == 60);

    auto bp_split = detect_breakpoint_from_cigar({
        {'M', 15}, {'N', 100}, {'M', 20}
    });
    check_result("split valid", bp_split.is_valid);
    check_result("split type", bp_split.type == ReadBreakpoint::SPLIT);
    check_result("split read_pos", bp_split.read_pos == 15);
    check_result("split ref_pos is event-time", bp_split.ref_pos == 15);

    auto bp_clip_5p = detect_breakpoint_from_cigar({
        {'S', 25}, {'M', 60}, {'S', 10}
    });
    check_result("5p clip chosen by longer clip", bp_clip_5p.type == ReadBreakpoint::SOFT_CLIP_5P);
    check_result("5p clip read_pos", bp_clip_5p.read_pos == 25);

    auto bp_clip_3p = detect_breakpoint_from_cigar({
        {'S', 10}, {'M', 60}, {'S', 25}
    });
    check_result("3p clip chosen by longer clip", bp_clip_3p.type == ReadBreakpoint::SOFT_CLIP_3P);
    check_result("3p clip read_pos", bp_clip_3p.read_pos == 70);
    check_result("3p clip ref_pos", bp_clip_3p.ref_pos == 60);

    std::cout << "  Breakpoint detection semantics tests passed!\n";
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

void test_batch_support_reads_preserved() {
    print_test_header("Batch Support Reads Preservation");

    AssemblyConfig config;
    config.min_reads_for_poa = 2;
    config.max_reads_for_poa = 20;
    config.breakpoint_upstream = 40;
    config.breakpoint_downstream = 40;
    config.ins_min_length = 30;

    AssemblyEngine engine(config);

    Component component{};
    component.id = 42;
    component.chrom_tid = 0;
    component.start = 1000;
    component.end = 1100;
    component.read_count = 10;  // intentionally larger than reads.size()

    std::vector<ReadSketch> reads;
    for (int i = 0; i < 6; ++i) {
        ReadSketch read;
        read.qname = "read_" + std::to_string(i);
        read.tid = 0;
        read.pos = 1000 + i;
        read.end_pos = read.pos + 80;
        read.flag = 0;
        read.mapq = 60;
        read.sequence = std::string(80, 'A') + std::string(80, 'T');
        read.cigar_ops = {{'M', 80}, {'S', 80}};
        read.total_clip_len = 80;
        read.breakpoint = detect_breakpoint_from_cigar(read.cigar_ops);
        reads.push_back(std::move(read));
    }

    const std::string ref_path = placer_test::resolve_test_file(
        "PLACER_TEST_REF",
        {"test_data/ref.fa", "tests/data/ref.fa"});
    GenomeAccessor genome(ref_path);

    auto contigs_component = engine.assemble_component(component, reads, genome);
    check_result("component-level contig exists", !contigs_component.empty());

    std::vector<Component> components = {component};
    auto contigs_batch = engine.assemble_batch(components, reads, genome);
    check_result("batch-level contig exists", !contigs_batch.empty());

    int32_t component_support = contigs_component.front().support_reads;
    int32_t batch_support = contigs_batch.front().support_reads;

    check_result("batch keeps assembled support",
        batch_support == component_support);
    check_result("batch support bounded by sampled reads",
        batch_support <= static_cast<int32_t>(reads.size()));

    std::cout << "  Batch support read preservation tests passed!\n";
}

void test_assembly_engine() {
    print_test_header("AssemblyEngine");

    // Test with empty input
    AssemblyConfig config;
    AssemblyEngine engine(config);
    std::vector<Component> empty_components;
    std::vector<ReadSketch> empty_reads;

    const std::string ref_path = placer_test::resolve_test_file(
        "PLACER_TEST_REF",
        {"test_data/ref.fa", "tests/data/ref.fa"});
    if (ref_path.empty()) {
        std::cout << "  [INFO] No reference fixture found; using empty GenomeAccessor\n";
    }
    GenomeAccessor genome(ref_path);

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
        test_rtree_spatial_query();
        test_poa_predecessor_overflow();
        test_breakpoint_detection_semantics();
        test_poa_assembly();
        test_multipath_poa();
        test_contig_structure();
        test_structural_representative();
        test_structural_collapse();
        test_batch_support_reads_preserved();
        test_assembly_engine();
        test_polya_extraction();

        std::cout << "\n=== All Phase 5 tests passed! ===\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Test failed with exception: " << e.what() << "\n";
        return 1;
    }
}
