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

bool env_try_bool(const char* key, bool& out) {
    const char* value = std::getenv(key);
    if (!value) {
        return false;
    }
    const std::string v(value);
    if (v == "1" || v == "true" || v == "TRUE" || v == "on" || v == "ON") {
        out = true;
        return true;
    }
    if (v == "0" || v == "false" || v == "FALSE" || v == "off" || v == "OFF") {
        out = false;
        return true;
    }
    return false;
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

bool env_try_string(const char* key, std::string& out) {
    const char* value = std::getenv(key);
    if (!value) {
        return false;
    }
    out = value;
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
    out << "final_te_certain\t" << result.final_te_certain << "\n";
    out << "final_te_uncertain\t" << result.final_te_uncertain << "\n";
    out << "final_non_te\t" << result.final_non_te << "\n";
    out << "final_high_confidence\t" << result.final_high_confidence << "\n";
    out << "final_low_confidence\t" << result.final_low_confidence << "\n";
    out << "bootstrap_exported_calls\t" << result.bootstrap_exported_calls << "\n";

    out << "\n#chrom\ttid\tpos\twindow_start\twindow_end\tte\tte_vote_frac\tte_median_ident\tte_fragments"
        << "\tte_theta\tte_mad_fwd\tte_mad_rev\tte_bp_core\tte_bp_win_start\tte_bp_win_end"
        << "\tte_core_candidates\tte_core_set\tsplit_sa_core_frac\tte_ref_junc_min\tte_ref_junc_max\tinsertion_qc\tte_qc"
        << "\tte_status\tte_top1_name\tte_top2_name\tte_post_top1\tte_post_top2\tte_post_margin\tte_conf_prob\tconfidence"
        << "\ttier\tsupport_reads\tgt\taf\tgq\tasm_mode\tasm_input_fragments\tasm_used_fragments"
        << "\tasm_consensus_len\tasm_identity_est\tasm_qc\n";
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
            << call.insertion_qc << "\t"
            << call.te_qc << "\t"
            << call.te_status << "\t"
            << (call.te_top1_name.empty() ? "NA" : call.te_top1_name) << "\t"
            << (call.te_top2_name.empty() ? "NA" : call.te_top2_name) << "\t"
            << call.te_posterior_top1 << "\t"
            << call.te_posterior_top2 << "\t"
            << call.te_posterior_margin << "\t"
            << call.te_confidence_prob << "\t"
            << call.confidence << "\t"
            << call.tier << "\t"
            << call.support_reads << "\t"
            << call.genotype << "\t"
            << call.af << "\t"
            << call.gq << "\t"
            << call.asm_mode << "\t"
            << call.asm_input_fragments << "\t"
            << call.asm_used_fragments << "\t"
            << call.asm_consensus_len << "\t"
            << call.asm_identity_est << "\t"
            << call.asm_qc << "\n";
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
    config.emit_low_confidence_calls = placer::env_flag_enabled("PLACER_EMIT_LOW_CONFIDENCE_CALLS");
    config.bootstrap_export_enable = placer::env_flag_enabled("PLACER_BOOTSTRAP_EXPORT");
    {
        double v = 0.0;
        bool b = false;
        if (placer::env_try_bool("PLACER_BOOTSTRAP_EXPORT_INCLUDE_NON_TE", b)) {
            config.bootstrap_export_include_non_te = b;
        }
        if (placer::env_try_double("PLACER_TE_MEDIAN_IDENTITY_MIN", v)) {
            config.te_median_identity_min = std::clamp(v, 0.0, 1.0);
        }
        if (placer::env_try_double("PLACER_TE_VOTE_FRACTION_MIN", v)) {
            config.te_vote_fraction_min = std::clamp(v, 0.0, 1.0);
        }
        if (placer::env_try_double("PLACER_TE_RESCUE_VOTE_FRACTION_MIN", v)) {
            config.te_rescue_vote_fraction_min = std::clamp(v, 0.0, 1.0);
        }
        if (placer::env_try_double("PLACER_TE_RESCUE_MEDIAN_IDENTITY_MIN", v)) {
            config.te_rescue_median_identity_min = std::clamp(v, 0.0, 1.0);
        }
        int32_t i = 0;
        if (placer::env_try_int32("PLACER_TE_KMER_SIZE", i)) {
            config.te_kmer_size = std::max(7, i);
        }
        if (placer::env_try_int32("PLACER_BAM_THREADS", i)) {
            config.bam_threads = std::max(1, i);
        }
        if (placer::env_try_int32("PLACER_TE_MIN_FRAGMENTS_FOR_VOTE", i)) {
            config.te_min_fragments_for_vote = std::max(1, i);
        }
        if (placer::env_try_int32("PLACER_PARALLEL_WORKERS", i)) {
            config.parallel_workers = std::max(1, i);
        }
        if (placer::env_try_int32("PLACER_TE_SOFTCLIP_LOW_COMPLEXITY_HOMOPOLYMER_MIN", i)) {
            config.te_softclip_low_complexity_homopolymer_min = std::max(1, i);
        }
        if (placer::env_try_int32("PLACER_PURE_SOFTCLIP_MIN_READS", i)) {
            config.te_pure_softclip_min_reads = std::max(1, i);
        }
        if (placer::env_try_int32("PLACER_PURE_SOFTCLIP_MIN_FRAGMENTS", i)) {
            config.te_pure_softclip_min_fragments = std::max(1, i);
        }
        if (placer::env_try_int32("PLACER_LOW_CONF_MIN_SUPPORT_READS", i)) {
            config.low_conf_min_support_reads = std::max(1, i);
        }
        if (placer::env_try_int32("PLACER_LOW_CONF_MAX_TIER", i)) {
            config.low_conf_max_tier = std::clamp(i, 1, 3);
        }
        if (placer::env_try_int32("PLACER_BOOTSTRAP_MIN_CONSENSUS_LEN", i)) {
            config.bootstrap_export_min_consensus_len = std::max(20, i);
        }

        if (placer::env_try_double("PLACER_TE_SOFTCLIP_LOW_COMPLEXITY_AT_FRAC_MIN", v)) {
            config.te_softclip_low_complexity_at_frac_min = std::clamp(v, 0.0, 1.0);
        }
        if (placer::env_try_double("PLACER_PURE_SOFTCLIP_MIN_IDENTITY", v)) {
            config.te_pure_softclip_min_identity = std::clamp(v, 0.0, 1.0);
        }
        if (placer::env_try_double("PLACER_TE_PROXY_POST_TOP1_MIN", v)) {
            config.te_proxy_posterior_top1_min = std::clamp(v, 0.0, 1.0);
        }
        if (placer::env_try_double("PLACER_TE_PROXY_POST_MARGIN_MIN", v)) {
            config.te_proxy_posterior_margin_min = std::clamp(v, 0.0, 1.0);
        }
        if (placer::env_try_double("PLACER_TE_PROXY_MIN_IDENTITY", v)) {
            config.te_proxy_identity_min = std::clamp(v, 0.0, 1.0);
        }
        if (placer::env_try_double("PLACER_TE_CONF_BIAS", v)) {
            config.te_confidence_bias = v;
        }
        if (placer::env_try_double("PLACER_TE_CONF_W_TOP1", v)) {
            config.te_confidence_w_top1 = v;
        }
        if (placer::env_try_double("PLACER_TE_CONF_W_MARGIN", v)) {
            config.te_confidence_w_margin = v;
        }
        if (placer::env_try_double("PLACER_TE_CONF_W_ASM_IDENTITY", v)) {
            config.te_confidence_w_asm_identity = v;
        }
        if (placer::env_try_double("PLACER_TE_CONF_W_SUPPORT", v)) {
            config.te_confidence_w_support = v;
        }
        if (placer::env_try_double("PLACER_TE_CONF_CERTAIN_MIN", v)) {
            config.te_confidence_prob_certain_min = std::clamp(v, 0.0, 1.0);
        }
        if (placer::env_try_double("PLACER_TE_CONF_UNCERTAIN_MIN", v)) {
            config.te_confidence_prob_uncertain_min = std::clamp(v, 0.0, 1.0);
        }
        std::string s;
        if (placer::env_try_string("PLACER_BOOTSTRAP_FASTA_PATH", s)) {
            config.bootstrap_consensus_fasta_path = s;
        }
        if (placer::env_try_string("PLACER_BOOTSTRAP_TSV_PATH", s)) {
            config.bootstrap_metadata_tsv_path = s;
        }
        if (placer::env_try_double("PLACER_EVIDENCE_MIN_SUPPORT_ALPHA", v)) {
            config.evidence_min_support_alpha = std::clamp(v, 0.0, 1.0);
        }
        if (placer::env_try_double("PLACER_EVIDENCE_MIN_SUPPORT_LAMBDA", v)) {
            config.evidence_min_support_lambda = std::clamp(v, 0.0, 1.0);
        }
        if (placer::env_try_double("PLACER_EVIDENCE_BREAKPOINT_MAD_MAX", v)) {
            config.evidence_breakpoint_mad_max = std::max(0.0, v);
        }
        if (placer::env_try_double("PLACER_EVIDENCE_LOW_COMPLEX_SOFTCLIP_FRAC_MAX", v)) {
            config.evidence_low_complex_softclip_frac_max = std::clamp(v, 0.0, 1.0);
        }
        if (placer::env_try_double("PLACER_EVIDENCE_TIER1_PROB", v)) {
            config.evidence_tier1_prob = std::clamp(v, 0.0, 1.0);
        }
        if (placer::env_try_double("PLACER_EVIDENCE_TIER2_PROB", v)) {
            config.evidence_tier2_prob = std::clamp(v, 0.0, 1.0);
        }
        if (placer::env_try_double("PLACER_EVIDENCE_LOGIT_BIAS", v)) {
            config.evidence_logit_bias = std::max(0.0, v);
        }
        if (config.evidence_tier2_prob > config.evidence_tier1_prob) {
            config.evidence_tier2_prob = config.evidence_tier1_prob;
        }
        if (placer::env_try_int32("PLACER_GENOTYPE_MIN_DEPTH", i)) {
            config.genotype_min_depth = std::max(1, i);
        }
        if (placer::env_try_double("PLACER_GENOTYPE_ERROR_RATE", v)) {
            config.genotype_error_rate = std::clamp(v, 1e-4, 0.25);
        }
        if (placer::env_try_int32("PLACER_ASSEMBLY_POA_MIN_READS", i)) {
            config.assembly_poa_min_reads = std::max(2, i);
        }
        if (placer::env_try_int32("PLACER_ASSEMBLY_POA_MAX_READS", i)) {
            config.assembly_poa_max_reads = std::max(2, i);
        }
        if (placer::env_try_int32("PLACER_ASSEMBLY_MIN_FRAGMENT_LEN", i)) {
            config.assembly_min_fragment_len = std::max(20, i);
        }
        if (placer::env_try_int32("PLACER_ASSEMBLY_MIN_CONSENSUS_LEN", i)) {
            config.assembly_min_consensus_len = std::max(20, i);
        }
        if (placer::env_try_int32("PLACER_ASSEMBLY_KMER_SIZE", i)) {
            config.assembly_kmer_size = std::max(5, i);
        }
        if (placer::env_try_double("PLACER_ASSEMBLY_MIN_IDENTITY_EST", v)) {
            config.assembly_min_identity_est = std::clamp(v, 0.0, 1.0);
        }
        config.te_rescue_vote_fraction_min = std::min(
            config.te_rescue_vote_fraction_min,
            config.te_vote_fraction_min);
        config.te_rescue_median_identity_min = std::min(
            config.te_rescue_median_identity_min,
            config.te_median_identity_min);
        if (config.te_confidence_prob_uncertain_min > config.te_confidence_prob_certain_min) {
            config.te_confidence_prob_uncertain_min = config.te_confidence_prob_certain_min;
        }
        config.assembly_poa_min_reads = std::min(
            config.assembly_poa_min_reads,
            config.assembly_poa_max_reads);
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
                  << "  genotype_calls=" << result.genotype_calls << "\n"
                  << "  final_te_certain=" << result.final_te_certain << "\n"
                  << "  final_te_uncertain=" << result.final_te_uncertain << "\n"
                  << "  final_non_te=" << result.final_non_te << "\n"
                  << "  final_high_confidence=" << result.final_high_confidence << "\n"
                  << "  final_low_confidence=" << result.final_low_confidence << "\n"
                  << "  bootstrap_exported_calls=" << result.bootstrap_exported_calls << std::endl;

        placer::write_scientific_txt(result, "scientific.txt");
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "[PLACER] fatal: " << ex.what() << std::endl;
        return 2;
    }
}
