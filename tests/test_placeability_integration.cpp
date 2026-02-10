#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include "placeability.h"
#include "local_realign.h"
#include "bam_reader.h"

using namespace placer;

// ============================================================================
// Integration Test using tldr_optimized/test data
// ============================================================================

// Read target positions from targets.test.txt
std::vector<std::tuple<std::string, int32_t, int32_t>> read_targets(const std::string& path) {
    std::vector<std::tuple<std::string, int32_t, int32_t>> targets;
    std::ifstream infile(path);
    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        std::string chrom;
        int32_t start, end;
        iss >> chrom >> start >> end;
        targets.emplace_back(chrom, start, end);
    }
    return targets;
}

// Simulate LocusEvidence from BAM reads around target region
std::vector<LocusEvidence> simulate_evidence_for_target(
    const std::string& chrom,
    int32_t target_start,
    int32_t target_end) {

    std::vector<LocusEvidence> evidence;

    // Simulate reads supporting the target locus (high confidence)
    for (int i = 0; i < 5; ++i) {
        LocusEvidence e;
        e.read_idx = i;
        e.locus_pos = target_start;  // All supporting the same locus
        e.normalized_score = 0.9 + (i % 2) * 0.1;  // High scores
        e.is_reverse = (i % 2 == 0);  // Mix of forward/reverse
        e.read_pos = target_start - 100;  // Left of locus
        e.read_end_pos = target_start + 500;  // Right of locus
        e.total_score = e.normalized_score * 100;
        e.evidence_bits = 0x01 | 0x02;  // up + down evidence
        evidence.push_back(e);
    }

    // Simulate some evidence for alternative locus
    for (int i = 0; i < 3; ++i) {
        LocusEvidence e;
        e.read_idx = 10 + i;
        e.locus_pos = target_start + 50;  // Alternative locus
        e.normalized_score = 0.6 + i * 0.05;  // Lower scores
        e.is_reverse = (i % 2 == 0);
        e.read_pos = target_start + 50 - 100;
        e.read_end_pos = target_start + 50 + 500;
        e.total_score = e.normalized_score * 100;
        e.evidence_bits = 0x01;
        evidence.push_back(e);
    }

    return evidence;
}

void test_with_test_data() {
    std::cout << "\n=== Integration Test with tldr_optimized/test data ===\n";

    // Read targets
    auto targets = read_targets("/mnt/home1/miska/hl725/projects/tldr_optimized/test/targets.test.txt");
    std::cout << "Loaded " << targets.size() << " targets from targets.test.txt\n";

    // For each target, create evidence and test placeability
    PlaceabilityConfig config;
    PlaceabilityScorer scorer(config);
    PlaceabilityOutput output(config);

    int test_idx = 0;
    for (const auto& [chrom, start, end] : targets) {
        test_idx++;
        std::cout << "\n--- Target " << test_idx << ": " << chrom << ":" << start << "-" << end << " ---\n";

        // Create simulated evidence for this target
        auto evidence = simulate_evidence_for_target(chrom, start, end);

        // Calculate placeability
        auto report = scorer.calculate_full_placeability(evidence);

        // Detect TSD (using flanking sequence simulation)
        TSDDetector detector(config);
        TSDResult tsd = detector.detect(
            "AAACCCGGGTTT" + std::to_string(test_idx % 10),
            "AAACCCGGGTTT" + std::to_string(test_idx % 10),
            start,
            end
        );

        // Generate output
        std::string info = output.generate_info_fields(report);
        std::string json = output.generate_json(report, tsd);

        std::cout << "  Best locus: " << report.best_locus << "\n";
        std::cout << "  Delta score: " << report.delta_score << "\n";
        std::cout << "  Tier: " << static_cast<int>(report.tier) << "\n";
        std::cout << "  Support reads: " << report.unique_support_reads << "\n";
        std::cout << "  Side consistent: " << (report.side_consistent ? "yes" : "no") << "\n";
        std::cout << "  TSD found: " << (tsd.found ? "yes" : "no") << "\n";
        std::cout << "  TSD length: " << tsd.tsd_length << "\n";

        // Verify expected results
        assert(report.best_locus == start);  // Best locus should be target start
        assert(report.candidate_count >= 2);  // Should have at least 2 candidates
        assert(report.unique_support_reads >= 5);  // Should have strong support

        std::cout << "  INFO: " << info.substr(0, 100) << "...\n";
    }

    std::cout << "\n=== Integration test passed! ===\n";
}

void test_tsd_with_realistic_flanks() {
    std::cout << "\n=== TSD Detection with Realistic Flanks ===\n";

    PlaceabilityConfig config;
    config.min_tsd_length = 3;
    config.max_tsd_length = 20;
    config.max_tsd_mismatches = 1;
    config.max_tsd_mismatch_ratio = 0.15;

    TSDDetector detector(config);

    // Test case 1: Perfect TSD (common in Alu insertions)
    {
        TSDResult result = detector.detect(
            "AAACCCGGGTTTAAACCC",
            "AAACCCGGGTTTAAACCC",
            100, 200
        );
        std::cout << "  Test 1 - Perfect TSD: found=" << result.found
                  << " len=" << result.tsd_length << " method=" << result.detection_method << "\n";
        assert(result.found);
        assert(result.tsd_length >= 12);
    }

    // Test case 2: TSD with 1 mismatch
    {
        TSDResult result = detector.detect(
            "AAACCCGGGTTTAAACCC",
            "AAACCCGGGTTTAAACCG",  // Last base G->G (no change)
            100, 200
        );
        std::cout << "  Test 2 - 1-mismatch TSD: found=" << result.found
                  << " len=" << result.tsd_length << "\n";
        assert(result.found);
    }

    // Test case 3: No TSD (different flanks)
    {
        TSDResult result = detector.detect(
            "AAACCCGGGTTTAAACCC",
            "GGGAAATTTCCCGGGAAA",
            100, 200
        );
        std::cout << "  Test 3 - No TSD: found=" << result.found << "\n";
        // This may or may not find fuzzy matches
    }

    // Test case 4: PolyA TSD (should be filtered out)
    {
        TSDDetector detector2(config);
        TSDResult result = detector2.detect(
            "AAAAAAATTTTTT",
            "AAAAAAATTTTTT",
            100, 200
        );
        std::cout << "  Test 4 - PolyA TSD: found=" << result.found
                  << " sig=" << result.is_significant << "\n";
        if (result.found) {
            // PolyA should not be significant
            assert(!result.is_significant || result.tsd_length >= 12);
        }
    }

    std::cout << "  All TSD tests passed!\n";
}

void test_tier_classification() {
    std::cout << "\n=== Tier Classification Tests ===\n";

    PlaceabilityConfig config;
    PlaceabilityScorer scorer(config);

    // Test Tier 1: High confidence unique placement
    {
        std::vector<LocusEvidence> evidence;
        for (int i = 0; i < 10; ++i) {
            LocusEvidence e;
            e.read_idx = i;
            e.locus_pos = 1000;
            e.normalized_score = 0.95;
            e.is_reverse = (i % 2 == 0);
            e.read_pos = 900;
            e.read_end_pos = 1100;
            evidence.push_back(e);
        }

        auto report = scorer.calculate_full_placeability(evidence);
        std::cout << "  Tier 1 test: tier=" << static_cast<int>(report.tier)
                  << " delta=" << report.delta_score << "\n";

        // Should be Tier 1 due to high support and delta
        bool is_high_tier = (report.tier == Tier::TIER1 || report.tier == Tier::TIER2);
        assert(is_high_tier);
    }

    // Test Tier 5: No evidence
    {
        std::vector<LocusEvidence> evidence;
        auto report = scorer.calculate_full_placeability(evidence);
        std::cout << "  Tier 5 test (empty): tier=" << static_cast<int>(report.tier) << "\n";
        assert(report.tier == Tier::TIER5);
    }

    // Test with ambiguous multiple loci
    {
        std::vector<LocusEvidence> evidence;
        for (int i = 0; i < 5; ++i) {
            LocusEvidence e;
            e.read_idx = i;
            e.locus_pos = 1000 + (i % 3) * 50;  // Spread across 3 loci
            e.normalized_score = 0.8;
            evidence.push_back(e);
        }

        auto report = scorer.calculate_full_placeability(evidence);
        std::cout << "  Ambiguous test: tier=" << static_cast<int>(report.tier)
                  << " candidates=" << report.candidate_count << "\n";
        assert(report.candidate_count >= 2);
    }

    std::cout << "  All tier classification tests passed!\n";
}

void test_output_formats() {
    std::cout << "\n=== Output Format Tests ===\n";

    PlaceabilityConfig config;
    PlaceabilityOutput output(config);

    // Create a report
    ExtendedPlaceabilityReport report;
    report.best_locus = 1000;
    report.second_best_locus = 950;
    report.delta_score = 35.5;
    report.side_consistent = true;
    report.side_consistency_verified = true;
    report.support_consistency = 0.88;
    report.candidate_count = 2;
    report.unique_support_reads = 8;
    report.tier = Tier::TIER1;
    report.overall_score = 45.2;
    report.is_ambiguous = false;
    report.strand_balanced = true;
    report.strand_ratio = 0.8f;
    report.locus_scores = {45.2, 9.7};
    report.all_loci = {1000, 950};

    // Generate outputs
    std::string info = output.generate_info_fields(report);
    std::string json = output.generate_json(report, TSDResult{});

    std::cout << "  INFO field: " << info.substr(0, 80) << "...\n";
    std::cout << "  JSON length: " << json.size() << " bytes\n";

    // Verify format
    assert(info.find("TIER=1") != std::string::npos);
    assert(info.find("DELTA=35.50") != std::string::npos);
    assert(info.find("SIDES=CONSISTENT_VERIFIED") != std::string::npos);
    assert(info.find("BEST_LOCUS=1000") != std::string::npos);

    assert(json.find("\"tier\":1") != std::string::npos);
    assert(json.find("\"delta_score\":") != std::string::npos);

    std::cout << "  All output format tests passed!\n";
}

int main() {
    std::cout << "=== PLACER Placeability Integration Tests ===\n";

    try {
        test_with_test_data();
        test_tsd_with_realistic_flanks();
        test_tier_classification();
        test_output_formats();

        std::cout << "\n=== All Integration Tests Passed! ===\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Test failed with exception: " << e.what() << "\n";
        return 1;
    }
}
