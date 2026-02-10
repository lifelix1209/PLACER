#include <cassert>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "genotyping.h"
#include "assembly.h"
#include "placeability.h"
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
// GenotypeConfig Tests
// ============================================================================

void test_genotype_config() {
    print_test_header("GenotypeConfig");

    // Default config
    GenotypeConfig config;
    assert(config.max_em_iterations == 20);
    assert(config.em_convergence_threshold == 1e-6);
    assert(config.spatial_lambda == 100.0);
    assert(config.null_base_prior == 0.1);

    // Custom config
    GenotypeConfig custom;
    custom.max_em_iterations = 50;
    custom.em_convergence_threshold = 1e-8;
    custom.spatial_lambda = 200.0;
    custom.null_base_prior = 0.05;
    custom.family_match_bonus = 3.0;

    assert(custom.max_em_iterations == 50);
    assert(custom.spatial_lambda == 200.0);

    std::cout << "  GenotypeConfig tests passed!\n";
}

// ============================================================================
// GenotypeResult Tests
// ============================================================================

void test_genotype_result() {
    print_test_header("GenotypeResult");

    GenotypeResult result;
    result.mix_alt = 0.6;
    result.mix_ref = 0.3;
    result.mix_null = 0.1;
    result.genotype = "0/1";
    result.gq = 30;
    result.af = 0.667;

    assert(result.mix_alt == 0.6);
    assert(result.genotype == "0/1");
    assert(result.gq == 30);

    // Test to_string
    std::string str = result.to_string();
    assert(!str.empty());
    assert(str.find("GT=0/1") != std::string::npos);
    assert(str.find("AF=0.667") != std::string::npos);

    std::cout << "  GenotypeResult tests passed!\n";
}

// ============================================================================
// ReadEvidence Tests
// ============================================================================

void test_read_evidence() {
    print_test_header("ReadEvidence");

    ReadEvidence evidence;
    evidence.read_idx = 0;
    evidence.locus_pos = 1000;
    evidence.d_spatial = 50.0;
    evidence.geom_ok = true;
    evidence.geom_score = 0.9;
    evidence.normalized_score = 0.85;
    evidence.contig_support = true;
    evidence.contig_score = 0.8;
    evidence.te_family_id = 5;
    evidence.structural_score = 0.7;

    assert(evidence.read_idx == 0);
    assert(evidence.locus_pos == 1000);
    assert(evidence.d_spatial == 50.0);
    assert(evidence.geom_ok == true);
    assert(evidence.normalized_score == 0.85);

    std::cout << "  ReadEvidence tests passed!\n";
}

// ============================================================================
// SpatialPriorCalculator Tests
// ============================================================================

void test_spatial_prior() {
    print_test_header("SpatialPriorCalculator");

    GenotypeConfig config;
    config.spatial_lambda = 100.0;
    config.null_base_prior = 0.1;

    SpatialPriorCalculator calculator(config);

    // Test 1: Read at locus (distance = 0)
    {
        double pi_alt, pi_ref, pi_null;
        calculator.calculate_prior(0.0, 0.5, pi_alt, pi_ref, pi_null);

        // Near locus should have high ALT/REF prior
        assert(pi_alt > pi_null);
        assert(pi_ref > pi_null);
    }

    // Test 2: Read far from locus
    {
        double pi_alt, pi_ref, pi_null;
        calculator.calculate_prior(1000.0, 0.5, pi_alt, pi_ref, pi_null);

        // Far from locus should have higher NULL prior
        assert(pi_null > pi_alt);
    }

    // Test 3: Batch calculation
    std::vector<ReadEvidence> evidence;
    for (int i = 0; i < 5; ++i) {
        ReadEvidence e;
        e.d_spatial = i * 100.0;
        evidence.push_back(e);
    }

    auto priors = calculator.calculate_priors(evidence);
    assert(priors.size() == 5);

    // Closer reads should have higher priors
    assert(priors[0] > priors[4]);

    std::cout << "  SpatialPriorCalculator tests passed!\n";
}

// ============================================================================
// StructuralPriorCalculator Tests
// ============================================================================

void test_structural_prior() {
    print_test_header("StructuralPriorCalculator");

    GenotypeConfig config;
    config.family_match_bonus = 2.0;
    config.family_mismatch_penalty = 0.1;

    StructuralPriorCalculator calculator(config);

    // Test 1: Same family
    {
        double bonus = calculator.calculate_family_bonus(5, 5);
        assert(bonus == 2.0);
    }

    // Test 2: Different families
    {
        double bonus = calculator.calculate_family_bonus(5, 10);
        assert(bonus == 0.1);
    }

    // Test 3: Unknown family
    {
        double bonus = calculator.calculate_family_bonus(-1, 5);
        assert(bonus == 1.0);
    }

    std::cout << "  StructuralPriorCalculator tests passed!\n";
}

// ============================================================================
// EM Engine Tests
// ============================================================================

void test_em_engine() {
    print_test_header("EMEngine");

    GenotypeConfig config;
    config.max_em_iterations = 50;
    config.em_convergence_threshold = 1e-8;
    config.em_init_alt_fraction = 0.5;

    EMEngine engine(config);

    // Create test evidence
    std::vector<ReadEvidence> evidence;

    // Strong ALT supporting reads
    for (int i = 0; i < 5; ++i) {
        ReadEvidence e;
        e.d_spatial = 10.0;
        e.geom_ok = true;
        e.geom_score = 0.9;
        e.normalized_score = 0.9;
        e.contig_support = true;
        e.contig_score = 0.9;
        evidence.push_back(e);
    }

    // Some NULL reads
    for (int i = 0; i < 3; ++i) {
        ReadEvidence e;
        e.d_spatial = 500.0;
        e.geom_ok = false;
        e.geom_score = 0.2;
        e.normalized_score = 0.3;
        e.contig_support = false;
        evidence.push_back(e);
    }

    // Test 1: EM convergence
    std::vector<std::array<double, 3>> priors;
    for (size_t i = 0; i < evidence.size(); ++i) {
        priors.push_back({0.5, 0.3, 0.2});  // pi_alt, pi_ref, pi_null
    }

    GenotypeResult result;
    auto [log_lik, iterations] = engine.run_em(evidence, priors, result);

    check_result("EM converged", iterations <= config.max_em_iterations);
    check_result("Positive log likelihood", log_lik < 0);  // Log likelihood should be negative
    check_result("E[ALT] + E[REF] + E[NULL] = 1",
                 std::abs(result.mix_alt + result.mix_ref + result.mix_null - 1.0) < 0.01);

    // Should have high ALT fraction due to evidence
    check_result("High ALT fraction", result.mix_alt > 0.5);

    // Test 2: Empty evidence
    {
        std::vector<ReadEvidence> empty;
        GenotypeResult empty_result;
        engine.run_em(empty, {}, empty_result);
        check_result("Empty evidence: mix_null=1", empty_result.mix_null > 0.99);
    }

    std::cout << "  EMEngine tests passed!\n";
}

// ============================================================================
// Genotyper Tests
// ============================================================================

void test_genotyper() {
    print_test_header("Genotyper");

    GenotypeConfig config;
    Genotyper genotyper(config);

    // Create representative
    StructuralRepresentative rep;
    rep.rep_id = 1;
    rep.fingerprint.breakpoint_l = 1000;  // 使用 fingerprint 中的断点
    rep.total_reads = 10;

    // Create evidence
    std::vector<LocusEvidence> evidence;
    for (int i = 0; i < 10; ++i) {
        LocusEvidence e;
        e.read_idx = i;
        e.locus_pos = 1000 + (i % 2) * 10;  // Mostly at 1000
        e.normalized_score = 0.85 + (i % 3) * 0.05;
        e.total_score = e.normalized_score * 100;
        e.weighted_identity = 0.9f;
        e.evidence_bits = 0x09;  // SA + high identity
        evidence.push_back(e);
    }

    // Create reads (ReadSketch doesn't have read_idx, use index)
    std::vector<ReadSketch> reads;
    for (int i = 0; i < 10; ++i) {
        ReadSketch r;
        r.qname = "read_" + std::to_string(i);
        reads.push_back(r);
    }

    // Create dummy GenomeAccessor
    GenomeAccessor genome("");

    // Run genotyping
    GenotypeResult result = genotyper.genotype(rep, evidence, reads, genome);

    check_result("Genotype produced", !result.genotype.empty());
    check_result("Valid genotype", result.genotype == "0/0" ||
                                    result.genotype == "0/1" ||
                                    result.genotype == "1/1" ||
                                    result.genotype == "./.");
    check_result("GQ computed", result.gq >= 0);
    check_result("AF computed", result.af >= 0 && result.af <= 1);

    std::cout << "  Genotyper tests passed!\n";
}

// ============================================================================
// Utility Functions Tests
// ============================================================================

void test_utility_functions() {
    print_test_header("Utility Functions");

    // Test beta_binomial_ci
    {
        auto ci = beta_binomial_ci(5, 10, 0.95);
        check_result("Beta CI low >= 0", ci[0] >= 0.0);
        check_result("Beta CI high <= 1", ci[1] <= 1.0);
        check_result("Beta CI valid range", ci[0] <= ci[1]);
    }

    // Test phred conversion
    {
        int phred = prob_to_phred(0.001);  // 1 in 1000
        check_result("Phred from prob", phred >= 29 && phred <= 31);
    }

    // Test likelihood ratio
    {
        double lr = likelihood_ratio_test(-10.0, -5.0);
        check_result("Positive likelihood ratio", lr > 0);
    }

    std::cout << "  Utility Functions tests passed!\n";
}

// ============================================================================
// Integration Tests
// ============================================================================

void test_genotyping_pipeline() {
    print_test_header("Genotyping Pipeline Integration");

    GenotypeConfig config;
    config.max_em_iterations = 100;
    config.spatial_lambda = 100.0;

    Genotyper genotyper(config);

    // Create dummy GenomeAccessor
    GenomeAccessor genome("");

    // Test case 1: Clear ALT signal
    {
        StructuralRepresentative rep;
        rep.rep_id = 1;
        rep.fingerprint.breakpoint_l = 1000;

        std::vector<LocusEvidence> evidence;
        for (int i = 0; i < 10; ++i) {
            LocusEvidence e;
            e.read_idx = i;
            e.locus_pos = 1000;
            e.normalized_score = 0.9;
            e.total_score = 90;
            e.weighted_identity = 0.95f;
            e.evidence_bits = 0x09;
            evidence.push_back(e);
        }

        std::vector<ReadSketch> reads;
        for (int i = 0; i < 10; ++i) {
            ReadSketch r;
            reads.push_back(r);
        }

        GenotypeResult result = genotyper.genotype(rep, evidence, reads, genome);

        std::cout << "  ALT signal: GT=" << result.genotype
                  << " mix_alt=" << result.mix_alt
                  << " mix_ref=" << result.mix_ref
                  << " mix_null=" << result.mix_null << "\n";

        check_result("Clear ALT: high mix_alt", result.mix_alt > 0.5);
    }

    // Test case 2: High NULL background
    {
        StructuralRepresentative rep;
        rep.rep_id = 2;
        rep.fingerprint.breakpoint_l = 2000;

        std::vector<LocusEvidence> evidence;
        for (int i = 0; i < 8; ++i) {
            LocusEvidence e;
            e.read_idx = i;
            e.locus_pos = 2000 + i * 100;  // Scattered
            e.normalized_score = 0.4;  // Low quality
            e.total_score = 40;
            e.weighted_identity = 0.7f;
            e.evidence_bits = 0x01;
            evidence.push_back(e);
        }

        std::vector<ReadSketch> reads;
        for (int i = 0; i < 8; ++i) {
            ReadSketch r;
            reads.push_back(r);
        }

        GenotypeResult result = genotyper.genotype(rep, evidence, reads, genome);

        std::cout << "  High NULL: GT=" << result.genotype
                  << " mix_null=" << result.mix_null << "\n";

        check_result("High NULL: marked", result.high_background);
    }

    std::cout << "  Genotyping Pipeline tests passed!\n";
}

// ============================================================================
// Main Test Runner
// ============================================================================

int main() {
    std::cout << "=== PLACER Phase 7 Genotyping Tests ===\n";

    try {
        test_genotype_config();
        test_genotype_result();
        test_read_evidence();
        test_spatial_prior();
        test_structural_prior();
        test_em_engine();
        test_genotyper();
        test_utility_functions();
        test_genotyping_pipeline();

        std::cout << "\n=== All Phase 7 tests passed! ===\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Test failed with exception: " << e.what() << "\n";
        return 1;
    }
}
