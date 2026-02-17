#include "pipeline.h"

#include <algorithm>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>

namespace placer {
namespace {

bool env_flag_enabled(const char* key) {
    const char* value = std::getenv(key);
    if (!value) {
        return false;
    }
    const std::string v(value);
    return v == "1" || v == "true" || v == "TRUE" || v == "on" || v == "ON";
}

bool env_try_double(const char* key, double& out) {
    const char* value = std::getenv(key);
    if (!value || !*value) {
        return false;
    }
    char* end = nullptr;
    const double parsed = std::strtod(value, &end);
    if (end == value || (end && *end != '\0')) {
        return false;
    }
    out = parsed;
    return true;
}

bool env_try_int32(const char* key, int32_t& out) {
    const char* value = std::getenv(key);
    if (!value || !*value) {
        return false;
    }
    char* end = nullptr;
    const long parsed = std::strtol(value, &end, 10);
    if (end == value || (end && *end != '\0')) {
        return false;
    }
    out = static_cast<int32_t>(parsed);
    return true;
}

void write_scientific_txt(const PipelineResult& result, const std::string& output_path) {
    std::ofstream out(output_path);
    if (!out.is_open()) {
        throw std::runtime_error("Failed to open output file: " + output_path);
    }

    out << "#PLACER streaming pipeline summary\n";
    out << "total_reads\t" << result.total_reads << "\n";
    out << "gate1_passed\t" << result.gate1_passed << "\n";
    out << "processed_bins\t" << result.processed_bins << "\n";
    out << "components\t" << result.built_components << "\n";
    out << "evidence_rows\t" << result.evidence_rows << "\n";
    out << "assemblies\t" << result.assembled_calls << "\n";
    out << "placeability_calls\t" << result.placeability_calls << "\n";
    out << "genotype_calls\t" << result.genotype_calls << "\n";

    out << "\n#chrom\ttid\tpos\twindow_start\twindow_end\tte\tte_vote_frac\tte_median_ident\tte_fragments"
        << "\tte_theta\tte_mad_fwd\tte_mad_rev\tte_bp_core\tte_bp_win_start\tte_bp_win_end"
        << "\tte_core_candidates\tte_core_set\tsplit_sa_core_frac\tte_ref_junc_min\tte_ref_junc_max\tte_qc"
        << "\ttier\tsupport_reads\tgt\taf\tgq\n";
    for (const auto& call : result.final_calls) {
        out << call.chrom << "\t"
            << call.tid << "\t"
            << call.pos << "\t"
            << call.window_start << "\t"
            << call.window_end << "\t"
            << (call.te_name.empty() ? "NA" : call.te_name) << "\t"
            << call.te_vote_fraction << "\t"
            << call.te_median_identity << "\t"
            << call.te_fragment_count << "\t"
            << call.te_theta << "\t"
            << call.te_mad_fwd << "\t"
            << call.te_mad_rev << "\t"
            << call.te_breakpoint_core << "\t"
            << call.te_breakpoint_window_start << "\t"
            << call.te_breakpoint_window_end << "\t"
            << call.te_core_candidates << "\t"
            << call.te_core_set << "\t"
            << call.te_split_sa_core_frac << "\t"
            << call.te_ref_junc_pos_min << "\t"
            << call.te_ref_junc_pos_max << "\t"
            << call.te_qc << "\t"
            << call.tier << "\t"
            << call.support_reads << "\t"
            << call.genotype << "\t"
            << call.af << "\t"
            << call.gq << "\n";
    }
}

}  // namespace
}  // namespace placer

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: placer <input.bam> <ref.fa> [te.fa]" << std::endl;
        return 1;
    }

    placer::PipelineConfig config;
    config.bam_path = argv[1];
    config.reference_fasta_path = argv[2];
    if (argc >= 4) {
        config.te_fasta_path = argv[3];
    }

    config.enable_parallel = placer::env_flag_enabled("PLACER_PARALLEL");
    {
        double v = 0.0;
        if (placer::env_try_double("PLACER_TE_MEDIAN_IDENTITY_MIN", v)) {
            config.te_median_identity_min = std::clamp(v, 0.0, 1.0);
        }
        if (placer::env_try_double("PLACER_TE_VOTE_FRACTION_MIN", v)) {
            config.te_vote_fraction_min = std::clamp(v, 0.0, 1.0);
        }
        int32_t i = 0;
        if (placer::env_try_int32("PLACER_TE_MIN_FRAGMENTS_FOR_VOTE", i)) {
            config.te_min_fragments_for_vote = std::max(1, i);
        }
    }

    try {
        auto pipeline = placer::build_default_pipeline(config);
        placer::PipelineResult result = pipeline->run();

        std::cerr << "[PLACER] pipeline finished\n"
                  << "  mode=" << (config.enable_parallel ? "parallel" : "streaming") << "\n"
                  << "  total_reads=" << result.total_reads << "\n"
                  << "  gate1_passed=" << result.gate1_passed << "\n"
                  << "  processed_bins=" << result.processed_bins << "\n"
                  << "  components=" << result.built_components << "\n"
                  << "  assemblies=" << result.assembled_calls << "\n"
                  << "  genotype_calls=" << result.genotype_calls << std::endl;

        placer::write_scientific_txt(result, "scientific.txt");
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "[PLACER] fatal: " << ex.what() << std::endl;
        return 2;
    }
}
