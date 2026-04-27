#include "denovo.h"
#include "pipeline.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
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

std::string safe_current_path() {
    try {
        return std::filesystem::current_path().string();
    } catch (const std::exception&) {
        return ".";
    }
}

std::string safe_absolute_path(const std::string& path) {
    try {
        return std::filesystem::absolute(path).string();
    } catch (const std::exception&) {
        return path;
    }
}

void print_usage() {
    std::cerr << "Usage: placer [--contig <chrom> | --region <chrom:start-end>] "
              << "<input.bam> <ref.fa> <te.fa>" << std::endl;
}

bool parse_positive_i64(const std::string& text, int64_t& out) {
    if (text.empty()) {
        return false;
    }
    for (char ch : text) {
        if (!std::isdigit(static_cast<unsigned char>(ch))) {
            return false;
        }
    }
    try {
        out = std::stoll(text);
    } catch (const std::exception&) {
        return false;
    }
    return out > 0;
}

BamRegionScope parse_region_scope_or_throw(const std::string& region) {
    if (region.empty()) {
        throw std::runtime_error("empty region");
    }

    BamRegionScope scope;
    scope.enabled = true;

    const size_t colon = region.rfind(':');
    if (colon == std::string::npos) {
        scope.chrom = region;
        scope.start = 0;
        scope.end = -1;
        return scope;
    }

    scope.chrom = region.substr(0, colon);
    const std::string range = region.substr(colon + 1);
    const size_t dash = range.find('-');
    if (scope.chrom.empty() || dash == std::string::npos) {
        throw std::runtime_error("invalid region: " + region);
    }

    int64_t start_1based = 0;
    int64_t end_1based = 0;
    if (!parse_positive_i64(range.substr(0, dash), start_1based) ||
        !parse_positive_i64(range.substr(dash + 1), end_1based) ||
        end_1based < start_1based ||
        end_1based > std::numeric_limits<int32_t>::max()) {
        throw std::runtime_error("invalid region: " + region);
    }

    scope.start = static_cast<int32_t>(start_1based - 1);
    scope.end = static_cast<int32_t>(end_1based);
    return scope;
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
    out << "event_consensus_calls\t" << result.event_consensus_calls << "\n";
    out << "genotype_calls\t" << result.genotype_calls << "\n";
    out << "final_pass_calls\t" << result.final_pass_calls << "\n";
    out << "schema_version\t1.0.0\n";

    out << "\n#chrom\tpos\tbp_left\tbp_right\tte\tfamily\tsubfamily\tstrand\tinsert_len"
        << "\tsupport_reads\talt_struct_reads\tref_span_reads\tlow_mapq_ref_span_reads"
        << "\tgt\taf\tgq"
        << "\tbest_te_identity\tbest_te_query_coverage\tcross_family_margin"
        << "\ttsd_type\ttsd_len"
        << "\tleft_flank_align_len\tright_flank_align_len\tconsensus_len"
        << "\tqc\n";
    for (const auto& call : result.final_calls) {
        out << call.chrom << "\t"
            << call.pos << "\t"
            << call.bp_left << "\t"
            << call.bp_right << "\t"
            << (call.te_name.empty() ? "NA" : call.te_name) << "\t"
            << call.family << "\t"
            << call.subfamily << "\t"
            << call.strand << "\t"
            << call.insert_len << "\t"
            << call.support_reads << "\t"
            << call.alt_struct_reads << "\t"
            << call.ref_span_reads << "\t"
            << call.low_mapq_ref_span_reads << "\t"
            << call.genotype << "\t"
            << call.af << "\t"
            << call.gq << "\t"
            << call.best_te_identity << "\t"
            << call.best_te_query_coverage << "\t"
            << call.cross_family_margin << "\t"
            << call.tsd_type << "\t"
            << call.tsd_len << "\t"
            << call.left_flank_align_len << "\t"
            << call.right_flank_align_len << "\t"
            << call.event_consensus_len << "\t"
            << call.final_qc << "\n";
    }
}

}  // namespace
}  // namespace placer

int main(int argc, char** argv) {
    if (argc >= 2 && std::string(argv[1]) == "denovo") {
        return placer::run_denovo_cli(argc - 1, argv + 1);
    }

    placer::BamRegionScope bam_region_scope;
    int argi = 1;
    while (argi < argc) {
        const std::string arg = argv[argi];
        if (arg == "--region") {
            if (bam_region_scope.enabled || argi + 1 >= argc) {
                placer::print_usage();
                return 1;
            }
            try {
                bam_region_scope = placer::parse_region_scope_or_throw(argv[argi + 1]);
            } catch (const std::exception& ex) {
                std::cerr << "[PLACER] invalid --region: " << ex.what() << std::endl;
                return 1;
            }
            argi += 2;
            continue;
        }
        if (arg == "--contig") {
            if (bam_region_scope.enabled || argi + 1 >= argc) {
                placer::print_usage();
                return 1;
            }
            bam_region_scope.enabled = true;
            bam_region_scope.chrom = argv[argi + 1];
            bam_region_scope.start = 0;
            bam_region_scope.end = -1;
            argi += 2;
            continue;
        }
        if (arg.rfind("--", 0) == 0) {
            placer::print_usage();
            return 1;
        }
        break;
    }

    if (argc - argi != 3) {
        placer::print_usage();
        return 1;
    }

    placer::PipelineConfig config;
    config.bam_path = argv[argi];
    config.reference_fasta_path = argv[argi + 1];
    config.te_fasta_path = argv[argi + 2];
    config.bam_region_scope = bam_region_scope;

    config.enable_parallel = placer::env_flag_enabled("PLACER_PARALLEL");
    {
        double v = 0.0;
        bool b = false;
        int32_t i = 0;
        if (placer::env_try_int32("PLACER_PROGRESS_INTERVAL", i)) {
            config.progress_interval = std::max<int64_t>(0, static_cast<int64_t>(i));
        }
        if (placer::env_try_bool("PLACER_LOG_STAGE_BINS", b)) {
            config.log_stage_bins = b;
        }
        if (placer::env_try_bool("PLACER_LOG_STAGE_COMPONENTS", b)) {
            config.log_stage_components = b;
        }
        std::string s;
        if (placer::env_try_string("PLACER_INS_FRAGMENTS_FASTA_PATH", s)) {
            config.ins_fragments_fasta_path = s;
        }
        if (placer::env_try_string("PLACER_INS_FRAGMENT_HITS_TSV_PATH", s)) {
            config.ins_fragment_hits_tsv_path = s;
        }
        if (placer::env_try_int32("PLACER_TE_KMER_SIZE", i)) {
            config.te_kmer_size = std::max(7, i);
        }
        if (placer::env_try_string("PLACER_TE_KMER_SIZES", s)) {
            config.te_kmer_sizes_csv = s;
        }
        if (placer::env_try_int32("PLACER_TE_FAMILY_TOPN", i)) {
            config.te_family_topn = std::max(1, i);
        }
        if (placer::env_try_int32("PLACER_TE_FAMILY_REPRESENTATIVES", i)) {
            config.te_family_representatives = std::max(1, i);
        }
        if (placer::env_try_int32("PLACER_TE_TEMPLATE_REFINE_TOPN", i)) {
            config.te_template_refine_topn = std::max(1, i);
        }
        if (placer::env_try_int32("PLACER_TE_EXACT_ALIGN_TOPN", i)) {
            config.te_exact_align_topn = std::max(1, i);
        }
        if (placer::env_try_double("PLACER_TE_FAMILY_MARGIN_MIN", v)) {
            config.te_family_margin_min = std::clamp(v, 0.0, 1.0);
        }
        if (placer::env_try_double("PLACER_TE_SUBFAMILY_MARGIN_MIN", v)) {
            config.te_subfamily_margin_min = std::clamp(v, 0.0, 1.0);
        }
        if (placer::env_try_int32("PLACER_BAM_THREADS", i)) {
            config.bam_threads = std::max(1, i);
        }
        if (placer::env_try_int32("PLACER_PARALLEL_WORKERS", i)) {
            config.parallel_workers = std::max(1, i);
        }
        if (placer::env_try_int32("PLACER_PARALLEL_QUEUE_MAX_TASKS", i)) {
            config.parallel_queue_max_tasks = i;
        }
        if (placer::env_try_int32("PLACER_PARALLEL_RESULT_BUFFER_MAX", i)) {
            config.parallel_result_buffer_max = i;
        }
        if (placer::env_try_bool("PLACER_LOG_PARALLEL_PROGRESS", b)) {
            config.log_parallel_progress = b;
        }
        if (placer::env_try_bool("PLACER_TSD_ENABLE", b)) {
            config.tsd_enable = b;
        }
        if (placer::env_try_int32("PLACER_TSD_MIN_LEN", i)) {
            config.tsd_min_len = std::max(1, i);
        }
        if (placer::env_try_int32("PLACER_TSD_MAX_LEN", i)) {
            config.tsd_max_len = std::max(1, i);
        }
        if (placer::env_try_int32("PLACER_TSD_FLANK_WINDOW", i)) {
            config.tsd_flank_window = std::max(10, i);
        }

        if (placer::env_try_double("PLACER_TSD_BG_P_MAX", v)) {
            config.tsd_bg_p_max = std::clamp(v, 0.0, 1.0);
        }
        if (placer::env_try_int32("PLACER_GENOTYPE_MIN_DEPTH", i)) {
            config.genotype_min_depth = std::max(1, i);
        }
        if (placer::env_try_double("PLACER_GENOTYPE_ERROR_RATE", v)) {
            config.genotype_error_rate = std::clamp(v, 1e-4, 0.25);
        }
        if (placer::env_try_int32("PLACER_EVENT_CONSENSUS_POA_MIN_READS", i)) {
            config.event_consensus_poa_min_reads = std::max(2, i);
        }
        if (placer::env_try_int32("PLACER_EVENT_CONSENSUS_POA_MAX_READS", i)) {
            config.event_consensus_poa_max_reads = std::max(2, i);
        }
        config.tsd_max_len = std::max(config.tsd_min_len, config.tsd_max_len);
        config.event_consensus_poa_min_reads = std::min(
            config.event_consensus_poa_min_reads,
            config.event_consensus_poa_max_reads);
    }

    try {
        const std::string scientific_path = placer::safe_absolute_path("scientific.txt");
        std::cerr << "[PLACER] run started\n"
                  << "  cwd=" << placer::safe_current_path() << "\n"
                  << "  bam=" << config.bam_path << "\n"
                  << "  bam_region=" << (config.bam_region_scope.enabled
                                           ? config.bam_region_scope.chrom + ":" +
                                                 std::to_string(config.bam_region_scope.start + 1) + "-" +
                                                 (config.bam_region_scope.end > 0
                                                      ? std::to_string(config.bam_region_scope.end)
                                                      : std::string("end"))
                                           : std::string("ALL")) << "\n"
                  << "  reference=" << config.reference_fasta_path << "\n"
                  << "  te_fasta=" << (config.te_fasta_path.empty() ? "NA" : config.te_fasta_path) << "\n"
                  << "  mode=" << (config.enable_parallel ? "parallel" : "streaming") << "\n"
                  << "  bam_threads=" << config.bam_threads << "\n"
                  << "  parallel_workers=" << config.parallel_workers << "\n"
                  << "  parallel_queue_max_tasks=" << config.parallel_queue_max_tasks << "\n"
                  << "  parallel_result_buffer_max=" << config.parallel_result_buffer_max << "\n"
                  << "  log_parallel_progress=" << (config.log_parallel_progress ? 1 : 0) << "\n"
                  << "  progress_interval=" << config.progress_interval << "\n"
                  << "  log_stage_bins=" << (config.log_stage_bins ? 1 : 0) << "\n"
                  << "  log_stage_components=" << (config.log_stage_components ? 1 : 0) << "\n"
                  << "  scientific_txt=" << scientific_path << std::endl;

        auto pipeline = placer::build_default_pipeline(config);
        placer::PipelineResult result = pipeline->run();

        std::cerr << "[PLACER] pipeline finished\n"
                  << "  mode=" << (config.enable_parallel ? "parallel" : "streaming") << "\n"
                  << "  total_reads=" << result.total_reads << "\n"
                  << "  gate1_passed=" << result.gate1_passed << "\n"
                  << "  processed_bins=" << result.processed_bins << "\n"
                  << "  components=" << result.built_components << "\n"
                  << "  event_consensus_calls=" << result.event_consensus_calls << "\n"
                  << "  genotype_calls=" << result.genotype_calls << "\n"
                  << "  final_pass_calls=" << result.final_pass_calls << std::endl;

        placer::write_scientific_txt(result, "scientific.txt");
        std::cerr << "[PLACER] wrote scientific.txt path=" << scientific_path << std::endl;
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "[PLACER] fatal: " << ex.what() << std::endl;
        return 2;
    }
}
