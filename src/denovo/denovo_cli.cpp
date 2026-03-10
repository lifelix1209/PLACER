#include "denovo.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace placer {
namespace {

int32_t env_int_or_default(const char* key, int32_t fallback) {
    const char* raw = std::getenv(key);
    if (raw == nullptr) {
        return fallback;
    }
    try {
        return std::stoi(raw);
    } catch (...) {
        return fallback;
    }
}

std::string trim_copy(std::string value) {
    const auto first = std::find_if_not(value.begin(), value.end(), [](unsigned char ch) {
        return std::isspace(ch) != 0;
    });
    const auto last = std::find_if_not(value.rbegin(), value.rend(), [](unsigned char ch) {
        return std::isspace(ch) != 0;
    }).base();
    if (first >= last) {
        return "";
    }
    return std::string(first, last);
}

std::string to_upper_copy(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::toupper(ch));
    });
    return value;
}

bool parse_int32_string(const std::string& text, int32_t& out) {
    if (text.empty() || text == "NA") {
        return false;
    }
    try {
        const int parsed = std::stoi(text);
        out = static_cast<int32_t>(parsed);
        return true;
    } catch (...) {
        return false;
    }
}

std::vector<std::string> split_tab_line(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream ss(line);
    std::string field;
    while (std::getline(ss, field, '\t')) {
        fields.push_back(field);
    }
    if (!line.empty() && line.back() == '\t') {
        fields.emplace_back();
    }
    return fields;
}

std::string get_field(
    const std::vector<std::string>& fields,
    const std::unordered_map<std::string, size_t>& columns,
    const std::string& name,
    const std::string& fallback = "") {
    const auto it = columns.find(name);
    if (it == columns.end() || it->second >= fields.size()) {
        return fallback;
    }
    return fields[it->second];
}

int32_t get_int32_field(
    const std::vector<std::string>& fields,
    const std::unordered_map<std::string, size_t>& columns,
    const std::string& name,
    int32_t fallback = -1) {
    int32_t out = fallback;
    const std::string raw = get_field(fields, columns, name);
    if (parse_int32_string(raw, out)) {
        return out;
    }
    return fallback;
}

std::string join_strings(const std::vector<std::string>& values, const std::string& delim) {
    std::ostringstream out;
    for (size_t i = 0; i < values.size(); ++i) {
        if (i > 0) {
            out << delim;
        }
        out << values[i];
    }
    return out.str();
}

void write_denovo_usage(std::ostream& out) {
    out << "Usage:\n"
        << "  placer denovo --child-scientific <child.scientific.txt>"
        << " --parent-bam-list <parents.fofn> --ref <ref.fa> --te <te.fa>"
        << " [--out-prefix trio_denovo]\n\n"
        << "Options:\n"
        << "  --child-min-support-reads <int>\n"
        << "  --child-max-tier <int>\n"
        << "  --child-confidence <HIGH|UNCERTAIN|...>\n"
        << "  --child-allow-uncertain\n"
        << "  --fetch-window <int>\n"
        << "  --default-match-window <int>\n"
        << "  --max-match-window <int>\n"
        << "  --min-softclip-len <int>\n"
        << "  --min-insertion-len <int>\n"
        << "  --parent-mapq-min <int>\n"
        << "  --no-family-match-veto\n"
        << "  --no-emit-review-status\n"
        << "  --dry-run\n";
}

bool file_exists(const std::string& path) {
    std::ifstream in(path);
    return in.good();
}

std::vector<std::string> load_path_list_file(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("Failed to open list file: " + path);
    }

    std::vector<std::string> paths;
    std::string line;
    while (std::getline(in, line)) {
        line = trim_copy(line);
        if (line.empty() || line[0] == '#') {
            continue;
        }
        paths.push_back(line);
    }
    return paths;
}

std::unordered_map<std::string, size_t> build_column_map(const std::vector<std::string>& header) {
    std::unordered_map<std::string, size_t> columns;
    columns.reserve(header.size());
    for (size_t i = 0; i < header.size(); ++i) {
        std::string key = header[i];
        if (!key.empty() && key.front() == '#') {
            key.erase(key.begin());
        }
        columns.emplace(std::move(key), i);
    }
    return columns;
}

void require_columns(
    const std::unordered_map<std::string, size_t>& columns,
    const std::vector<std::string>& required) {
    for (const auto& name : required) {
        if (columns.find(name) == columns.end()) {
            throw std::runtime_error("scientific.txt missing required column: " + name);
        }
    }
}

DenovoChildCandidate parse_child_candidate_row(
    const std::vector<std::string>& fields,
    const std::unordered_map<std::string, size_t>& columns,
    const DenovoConfig& config) {
    DenovoChildCandidate candidate;
    candidate.chrom = get_field(fields, columns, "chrom");
    candidate.pos = get_int32_field(fields, columns, "pos");
    candidate.te_name = get_field(fields, columns, "te", "NA");
    candidate.te_status = get_field(fields, columns, "te_status", "NON_TE");
    candidate.confidence = get_field(fields, columns, "confidence", "NA");
    candidate.tier = get_int32_field(fields, columns, "tier", 3);
    candidate.support_reads = get_int32_field(fields, columns, "support_reads", 0);

    const int32_t ref_junc_min = get_int32_field(fields, columns, "te_ref_junc_min", -1);
    const int32_t ref_junc_max = get_int32_field(fields, columns, "te_ref_junc_max", -1);
    const int32_t bp_win_start = get_int32_field(fields, columns, "te_bp_win_start", -1);
    const int32_t bp_win_end = get_int32_field(fields, columns, "te_bp_win_end", -1);

    if (ref_junc_min >= 0 && ref_junc_max >= ref_junc_min) {
        candidate.event_start = ref_junc_min;
        candidate.event_end = ref_junc_max;
    } else if (bp_win_start >= 0 && bp_win_end >= bp_win_start) {
        candidate.event_start = bp_win_start;
        candidate.event_end = bp_win_end;
    } else {
        candidate.event_start = std::max(0, candidate.pos - config.default_match_window);
        candidate.event_end = candidate.pos + config.default_match_window;
    }

    return candidate;
}

bool candidate_passes_child_filters(const DenovoChildCandidate& candidate, const DenovoConfig& config) {
    const std::string te_name = to_upper_copy(candidate.te_name);
    if (te_name.empty() || te_name == "NA") {
        return false;
    }

    if (to_upper_copy(candidate.te_status) == "NON_TE") {
        return false;
    }

    if (candidate.support_reads < config.child_min_support_reads) {
        return false;
    }

    if (candidate.tier > config.child_max_tier) {
        return false;
    }

    const std::string confidence = to_upper_copy(candidate.confidence);
    const std::string required = to_upper_copy(config.child_confidence);
    if (!required.empty()) {
        if (config.child_allow_uncertain && required == "HIGH") {
            if (confidence != "HIGH" && confidence != "UNCERTAIN") {
                return false;
            }
        } else if (confidence != required) {
            return false;
        }
    }

    return true;
}

void write_calls_file(const DenovoConfig& config, const DenovoResult& result) {
    std::ofstream out(config.out_prefix + ".calls.tsv");
    if (!out.is_open()) {
        throw std::runtime_error("Failed to open denovo calls output");
    }

    out << "#chrom\tpos\tte\tchild_support_reads\tchild_confidence\tchild_tier"
        << "\tchild_event_start\tchild_event_end\tparent_total_support_reads"
        << "\tparent_exact_te_reads\tparent_family_te_reads\tparent_ambiguous_signal_reads"
        << "\tparent_support_bams\tscanner_status\tstatus\tde_novo\n";
    for (const auto& call : result.calls) {
        out << call.child.chrom << "\t"
            << call.child.pos << "\t"
            << call.child.te_name << "\t"
            << call.child.support_reads << "\t"
            << call.child.confidence << "\t"
            << call.child.tier << "\t"
            << call.child.event_start << "\t"
            << call.child.event_end << "\t"
            << call.parent_summary.total_support_reads << "\t"
            << call.parent_summary.exact_te_reads << "\t"
            << call.parent_summary.family_te_reads << "\t"
            << call.parent_summary.ambiguous_signal_reads << "\t"
            << (call.parent_summary.support_bams.empty()
                    ? "NA"
                    : join_strings(call.parent_summary.support_bams, ","))
            << "\t"
            << call.parent_summary.scanner_status << "\t"
            << call.status << "\t"
            << call.de_novo << "\n";
    }
}

void write_veto_reads_file(const DenovoConfig& config, const DenovoResult& result) {
    std::ofstream out(config.out_prefix + ".parent_veto_reads.tsv");
    if (!out.is_open()) {
        throw std::runtime_error("Failed to open denovo parent veto output");
    }

    out << "#chrom\tchild_pos\tchild_te\tparent_bam\tread_name\tevidence_type"
        << "\tbreakpoint_pos\tfragment_len\tmatched_te\tmatched_family\tveto_reason\n";
    for (const auto& call : result.calls) {
        for (const auto& evidence : call.parent_summary.veto_reads) {
            out << call.child.chrom << "\t"
                << call.child.pos << "\t"
                << call.child.te_name << "\t"
                << evidence.parent_bam_path << "\t"
                << evidence.read_name << "\t"
                << evidence.evidence_type << "\t"
                << evidence.breakpoint_pos << "\t"
                << evidence.fragment_len << "\t"
                << evidence.matched_te << "\t"
                << evidence.matched_family << "\t"
                << evidence.veto_reason << "\n";
        }
    }
}

void write_summary_file(const DenovoConfig& config, const DenovoResult& result) {
    std::ofstream out(config.out_prefix + ".summary.txt");
    if (!out.is_open()) {
        throw std::runtime_error("Failed to open denovo summary output");
    }

    out << "implementation_status\t" << result.implementation_status << "\n";
    out << "child_rows_total\t" << result.child_rows_total << "\n";
    out << "child_candidates_considered\t" << result.child_candidates_considered << "\n";
    out << "parent_bams\t" << result.parent_bams << "\n";
    out << "calls_written\t" << result.calls_written << "\n";
    out << "inherited_prefilter_hits\t" << result.inherited_prefilter_hits << "\n";
    out << "parent_veto_calls\t" << result.parent_veto_calls << "\n";
    out << "review_calls\t" << result.review_calls << "\n";
    out << "denovo_pass_calls\t" << result.denovo_pass_calls << "\n";
}

}  // namespace

bool parse_denovo_cli_args(
    int argc,
    char** argv,
    DenovoConfig& config,
    std::string& error_message) {
    config = DenovoConfig{};
    config.bam_threads = std::max<int32_t>(1, env_int_or_default("PLACER_BAM_THREADS", 1));

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto require_value = [&](const char* flag) -> std::string {
            if (i + 1 >= argc) {
                throw std::runtime_error(std::string("Missing value for ") + flag);
            }
            ++i;
            return argv[i];
        };

        try {
            if (arg == "--child-scientific") {
                config.child_scientific_path = require_value("--child-scientific");
            } else if (arg == "--parent-bam-list") {
                config.parent_bam_list_path = require_value("--parent-bam-list");
            } else if (arg == "--ref") {
                config.reference_fasta_path = require_value("--ref");
            } else if (arg == "--te") {
                config.te_fasta_path = require_value("--te");
            } else if (arg == "--out-prefix") {
                config.out_prefix = require_value("--out-prefix");
            } else if (arg == "--child-min-support-reads") {
                config.child_min_support_reads =
                    std::max<int32_t>(1, std::stoi(require_value("--child-min-support-reads")));
            } else if (arg == "--child-max-tier") {
                config.child_max_tier = std::max<int32_t>(1, std::stoi(require_value("--child-max-tier")));
            } else if (arg == "--child-confidence") {
                config.child_confidence = require_value("--child-confidence");
            } else if (arg == "--child-allow-uncertain") {
                config.child_allow_uncertain = true;
            } else if (arg == "--fetch-window") {
                config.fetch_window = std::max<int32_t>(1, std::stoi(require_value("--fetch-window")));
            } else if (arg == "--default-match-window") {
                config.default_match_window =
                    std::max<int32_t>(1, std::stoi(require_value("--default-match-window")));
            } else if (arg == "--max-match-window") {
                config.max_match_window =
                    std::max<int32_t>(config.default_match_window, std::stoi(require_value("--max-match-window")));
            } else if (arg == "--min-softclip-len") {
                config.min_softclip_len = std::max<int32_t>(1, std::stoi(require_value("--min-softclip-len")));
            } else if (arg == "--min-insertion-len") {
                config.min_insertion_len = std::max<int32_t>(1, std::stoi(require_value("--min-insertion-len")));
            } else if (arg == "--parent-mapq-min") {
                config.parent_mapq_min = std::max<int32_t>(0, std::stoi(require_value("--parent-mapq-min")));
            } else if (arg == "--no-family-match-veto") {
                config.family_match_veto = false;
            } else if (arg == "--no-emit-review-status") {
                config.emit_review_status = false;
            } else if (arg == "--dry-run") {
                config.dry_run = true;
            } else {
                error_message = "Unknown argument: " + arg;
                return false;
            }
        } catch (const std::exception& ex) {
            error_message = ex.what();
            return false;
        }
    }

    if (config.child_scientific_path.empty()) {
        error_message = "Missing required argument: --child-scientific";
        return false;
    }
    if (config.parent_bam_list_path.empty()) {
        error_message = "Missing required argument: --parent-bam-list";
        return false;
    }
    if (config.reference_fasta_path.empty()) {
        error_message = "Missing required argument: --ref";
        return false;
    }
    if (config.te_fasta_path.empty()) {
        error_message = "Missing required argument: --te";
        return false;
    }

    config.max_match_window = std::max(config.default_match_window, config.max_match_window);
    return true;
}

std::vector<DenovoChildCandidate> load_denovo_child_candidates(
    const DenovoConfig& config,
    int64_t* total_rows) {
    if (total_rows != nullptr) {
        *total_rows = 0;
    }
    if (!file_exists(config.child_scientific_path)) {
        throw std::runtime_error("child scientific.txt not found: " + config.child_scientific_path);
    }

    std::ifstream in(config.child_scientific_path);
    if (!in.is_open()) {
        throw std::runtime_error("Failed to open child scientific.txt: " + config.child_scientific_path);
    }

    std::string line;
    std::vector<std::string> header;
    std::unordered_map<std::string, size_t> columns;
    std::vector<DenovoChildCandidate> candidates;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (line.rfind("#chrom\t", 0) == 0) {
            header = split_tab_line(line);
            columns = build_column_map(header);
            require_columns(
                columns,
                {"chrom", "pos", "te", "te_status", "confidence", "tier", "support_reads",
                 "te_bp_win_start", "te_bp_win_end", "te_ref_junc_min", "te_ref_junc_max"});
            continue;
        }
        if (line[0] == '#') {
            continue;
        }
        if (header.empty()) {
            continue;
        }

        std::vector<std::string> fields = split_tab_line(line);
        if (total_rows != nullptr) {
            *total_rows += 1;
        }

        DenovoChildCandidate candidate = parse_child_candidate_row(fields, columns, config);
        if (!candidate_passes_child_filters(candidate, config)) {
            continue;
        }
        candidates.push_back(std::move(candidate));
    }

    if (header.empty()) {
        throw std::runtime_error("Failed to find #chrom header in child scientific.txt");
    }

    return candidates;
}

DenovoResult run_denovo(const DenovoConfig& config) {
    if (!file_exists(config.child_scientific_path)) {
        throw std::runtime_error("child scientific.txt not found: " + config.child_scientific_path);
    }
    if (!file_exists(config.parent_bam_list_path)) {
        throw std::runtime_error("parent BAM list not found: " + config.parent_bam_list_path);
    }
    if (!file_exists(config.reference_fasta_path)) {
        throw std::runtime_error("reference FASTA not found: " + config.reference_fasta_path);
    }
    if (!file_exists(config.te_fasta_path)) {
        throw std::runtime_error("TE FASTA not found: " + config.te_fasta_path);
    }

    DenovoConfig resolved = config;
    if (resolved.parent_bam_paths.empty()) {
        resolved.parent_bam_paths = load_path_list_file(resolved.parent_bam_list_path);
    }
    if (resolved.parent_bam_paths.empty()) {
        throw std::runtime_error("No parent BAMs found in list: " + resolved.parent_bam_list_path);
    }

    DenovoResult result;
    result.parent_bams = static_cast<int64_t>(resolved.parent_bam_paths.size());

    int64_t total_rows = 0;
    std::vector<DenovoChildCandidate> candidates = load_denovo_child_candidates(resolved, &total_rows);
    result.child_rows_total = total_rows;
    result.child_candidates_considered = static_cast<int64_t>(candidates.size());

    std::unique_ptr<ParentPoolScanner> scanner;
    if (resolved.dry_run) {
        result.implementation_status = "DRY_RUN";
    } else {
        scanner = std::make_unique<ParentPoolScanner>(resolved);
        result.implementation_status = scanner->implementation_status();
    }

    result.calls.reserve(candidates.size());
    for (const auto& candidate : candidates) {
        DenovoCall call;
        call.child = candidate;
        if (resolved.dry_run) {
            call.parent_summary.scanner_status = "DRY_RUN";
            call.status = "DRY_RUN";
            call.de_novo = "DRY_RUN";
        } else {
            call.parent_summary = scanner->scan_candidate(candidate);
            const bool parent_veto =
                call.parent_summary.exact_te_reads > 0 ||
                (resolved.family_match_veto && call.parent_summary.family_te_reads > 0);
            const bool review_signal =
                call.parent_summary.ambiguous_signal_reads > 0 ||
                (!resolved.family_match_veto && call.parent_summary.family_te_reads > 0);

            if (parent_veto) {
                call.status = "PARENT_VETO";
                call.de_novo = "0";
                result.parent_veto_calls += 1;
            } else if (review_signal && resolved.emit_review_status) {
                call.status = "PARENT_AMBIGUOUS_SIGNAL";
                call.de_novo = "REVIEW";
                result.review_calls += 1;
            } else {
                call.status = "DENOVO_PASS";
                call.de_novo = "1";
                result.denovo_pass_calls += 1;
            }
        }
        result.calls.push_back(std::move(call));
    }

    result.calls_written = static_cast<int64_t>(result.calls.size());
    return result;
}

void write_denovo_outputs(const DenovoConfig& config, const DenovoResult& result) {
    write_calls_file(config, result);
    write_veto_reads_file(config, result);
    write_summary_file(config, result);
}

int run_denovo_cli(int argc, char** argv) {
    if (argc <= 1) {
        write_denovo_usage(std::cerr);
        return 1;
    }

    const std::string first_arg = argv[1];
    if (first_arg == "--help" || first_arg == "-h") {
        write_denovo_usage(std::cout);
        return 0;
    }

    DenovoConfig config;
    std::string error_message;
    if (!parse_denovo_cli_args(argc, argv, config, error_message)) {
        std::cerr << "[DENOVO] " << error_message << "\n\n";
        write_denovo_usage(std::cerr);
        return 1;
    }

    try {
        DenovoResult result = run_denovo(config);
        write_denovo_outputs(config, result);

        std::cerr << "[DENOVO] skeleton finished\n"
                  << "  implementation_status=" << result.implementation_status << "\n"
                  << "  child_rows_total=" << result.child_rows_total << "\n"
                  << "  child_candidates_considered=" << result.child_candidates_considered << "\n"
                  << "  parent_bams=" << result.parent_bams << "\n"
                  << "  calls_written=" << result.calls_written << std::endl;
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "[DENOVO] fatal: " << ex.what() << std::endl;
        return 2;
    }
}

}  // namespace placer
