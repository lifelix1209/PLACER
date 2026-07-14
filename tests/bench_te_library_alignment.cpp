#include <chrono>
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unistd.h>
#include <vector>

#define private public
#include "pipeline.h"
#undef private

namespace {

std::string make_temp_path(const std::string& stem) {
    return "/tmp/" + stem + "_" + std::to_string(static_cast<long long>(::getpid())) + ".fa";
}

std::string build_sequence(size_t len, uint32_t seed) {
    static const char kBases[] = {'A', 'C', 'G', 'T'};
    uint32_t state = seed;
    std::string out;
    out.reserve(len);
    for (size_t i = 0; i < len; ++i) {
        state = (state * 1664525u) + 1013904223u;
        out.push_back(kBases[(state >> 24) & 3u]);
    }
    return out;
}

std::string mutate_positions(
    std::string seq,
    const std::vector<size_t>& positions) {
    for (size_t pos : positions) {
        if (pos >= seq.size()) {
            continue;
        }
        switch (seq[pos]) {
            case 'A': seq[pos] = 'C'; break;
            case 'C': seq[pos] = 'G'; break;
            case 'G': seq[pos] = 'T'; break;
            case 'T': seq[pos] = 'A'; break;
            default: seq[pos] = 'A'; break;
        }
    }
    return seq;
}

void write_te_fasta(
    const std::string& path,
    const std::vector<std::pair<std::string, std::string>>& entries) {
    std::ofstream out(path);
    for (const auto& entry : entries) {
        out << ">" << entry.first << "\n";
        out << entry.second << "\n";
    }
}

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> fields;
    size_t start = 0;
    while (start <= line.size()) {
        const size_t tab = line.find('\t', start);
        if (tab == std::string::npos) {
            fields.push_back(line.substr(start));
            break;
        }
        fields.push_back(line.substr(start, tab - start));
        start = tab + 1;
    }
    return fields;
}

char uppercase_base(char c) {
    c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    switch (c) {
        case 'A':
        case 'C':
        case 'G':
        case 'T':
        case 'N':
            return c;
        default:
            return '\0';
    }
}

std::string extract_blue_insert_sequence(const std::string& consensus) {
    std::string out;
    std::string ansi_code = "0";
    for (size_t i = 0; i < consensus.size();) {
        if (consensus[i] == '\x1b' && (i + 1) < consensus.size() &&
            consensus[i + 1] == '[') {
            const size_t end = consensus.find('m', i + 2);
            if (end == std::string::npos) {
                break;
            }
            ansi_code = consensus.substr(i + 2, end - (i + 2));
            i = end + 1;
            continue;
        }
        if (ansi_code == "34") {
            const char base = uppercase_base(consensus[i]);
            if (base != '\0') {
                out.push_back(base);
            }
        }
        ++i;
    }
    return out;
}

std::string display_or_na(const std::string& value) {
    return value.empty() ? "NA" : value;
}

bool subfamily_matches(const std::string& truth, const std::string& predicted) {
    if (truth == predicted) {
        return true;
    }
    if (predicted.empty() || predicted == "NA") {
        return false;
    }
    return false;
}

int run_truth_table_benchmark(
    const std::string& te_fasta,
    const std::string& truth_table,
    int32_t limit,
    int32_t seed) {
    placer::PipelineConfig cfg;
    cfg.te_fasta_path = te_fasta;
    cfg.te_kmer_size = 9;
    cfg.te_kmer_sizes_csv = "7,9,11";
    cfg.te_family_topn = 8;
    cfg.te_template_refine_topn = 8;
    cfg.te_exact_align_topn = 8;
    cfg.te_family_representatives = 4;
    placer::TEKmerQuickClassifierModule classifier(cfg);
    if (!classifier.is_enabled()) {
        std::cerr << "failed_to_load_te_library=" << te_fasta << "\n";
        return 2;
    }

    std::ifstream in(truth_table);
    if (!in.is_open()) {
        std::cerr << "failed_to_open_truth_table=" << truth_table << "\n";
        return 2;
    }

    std::string header_line;
    if (!std::getline(in, header_line)) {
        std::cerr << "empty_truth_table=" << truth_table << "\n";
        return 2;
    }
    const auto header = split_tab(header_line);
    std::unordered_map<std::string, size_t> column;
    for (size_t i = 0; i < header.size(); ++i) {
        column.emplace(header[i], i);
    }
    const auto require_col = [&](const std::string& name) -> size_t {
        const auto it = column.find(name);
        if (it == column.end()) {
            throw std::runtime_error("missing column: " + name);
        }
        return it->second;
    };

    const size_t uuid_col = require_col("UUID");
    const size_t family_col = require_col("Family");
    const size_t subfamily_col = require_col("Subfamily");
    const size_t filter_col = require_col("Filter");
    const size_t consensus_col = require_col("Consensus");

    struct TruthRow {
        std::string uuid;
        std::string family;
        std::string subfamily;
        std::string insert_sequence;
    };
    std::vector<TruthRow> rows;
    std::string line;
    while (std::getline(in, line)) {
        const auto fields = split_tab(line);
        if (fields.size() <= std::max({uuid_col, family_col, subfamily_col, filter_col, consensus_col})) {
            continue;
        }
        if (fields[filter_col] != "PASS" ||
            fields[family_col] == "NA" ||
            fields[family_col] == "Unknown" ||
            fields[subfamily_col] == "NA") {
            continue;
        }
        std::string insert_sequence =
            extract_blue_insert_sequence(fields[consensus_col]);
        if (insert_sequence.size() < 80) {
            continue;
        }
        rows.push_back({
            fields[uuid_col],
            fields[family_col],
            fields[subfamily_col],
            std::move(insert_sequence),
        });
    }
    if (rows.empty()) {
        std::cerr << "no_eligible_truth_rows\n";
        return 2;
    }

    uint32_t state = static_cast<uint32_t>(seed);
    for (size_t i = rows.size(); i > 1; --i) {
        state = (state * 1664525u) + 1013904223u;
        const size_t j = static_cast<size_t>(state % static_cast<uint32_t>(i));
        std::swap(rows[i - 1], rows[j]);
    }
    const int32_t n = std::min<int32_t>(
        std::max(1, limit),
        static_cast<int32_t>(rows.size()));

    int32_t pass_count = 0;
    int32_t family_matches = 0;
    int32_t subfamily_matches_count = 0;
    std::cout
        << "uuid\ttruth_family\ttruth_subfamily\tinsert_len\tpred_family"
        << "\tpred_subfamily\tidentity\tquery_coverage\tsecond_family"
        << "\tcross_family_margin\tcoarse_prefilter\tcoarse_chain"
        << "\tconfidence\tqc\tintervals\n";
    for (int32_t i = 0; i < n; ++i) {
        const TruthRow& row = rows[static_cast<size_t>(i)];
        const placer::TEAlignmentEvidence evidence =
            classifier.align_insert_sequence(row.insert_sequence);
        if (evidence.pass) {
            ++pass_count;
        }
        if (evidence.best_family == row.family) {
            ++family_matches;
        }
        if (subfamily_matches(row.subfamily, evidence.best_subfamily)) {
            ++subfamily_matches_count;
        }
        std::cout
            << row.uuid << "\t"
            << row.family << "\t"
            << row.subfamily << "\t"
            << row.insert_sequence.size() << "\t"
            << display_or_na(evidence.best_family) << "\t"
            << display_or_na(evidence.best_subfamily) << "\t"
            << evidence.best_identity << "\t"
            << evidence.best_query_coverage << "\t"
            << evidence.second_family << "\t"
            << evidence.cross_family_margin << "\t"
            << evidence.coarse_prefilter_score << "\t"
            << evidence.coarse_chain_coverage << "\t"
            << evidence.annotation_confidence << "\t"
            << evidence.qc_reason << "\t"
            << evidence.annotation_intervals << "\n";
    }

    std::cout
        << "summary_rows\t" << n << "\n"
        << "summary_pass\t" << pass_count << "\n"
        << "summary_family_matches\t" << family_matches << "\n"
        << "summary_subfamily_matches\t" << subfamily_matches_count << "\n"
        << "summary_family_accuracy\t"
        << (static_cast<double>(family_matches) / static_cast<double>(n)) << "\n"
        << "summary_subfamily_accuracy\t"
        << (static_cast<double>(subfamily_matches_count) / static_cast<double>(n)) << "\n";
    return 0;
}

int run_query_tsv_benchmark(
    const std::string& te_fasta,
    const std::string& query_tsv,
    int32_t limit) {
    placer::PipelineConfig cfg;
    cfg.te_fasta_path = te_fasta;
    cfg.te_kmer_size = 9;
    cfg.te_kmer_sizes_csv = "7,9,11";
    cfg.te_family_topn = 8;
    cfg.te_template_refine_topn = 8;
    cfg.te_exact_align_topn = 8;
    cfg.te_family_representatives = 4;
    placer::TEKmerQuickClassifierModule classifier(cfg);
    if (!classifier.is_enabled()) {
        std::cerr << "failed_to_load_te_library=" << te_fasta << "\n";
        return 2;
    }

    std::ifstream in(query_tsv);
    if (!in.is_open()) {
        std::cerr << "failed_to_open_query_tsv=" << query_tsv << "\n";
        return 2;
    }
    std::string line;
    if (!std::getline(in, line)) {
        std::cerr << "empty_query_tsv=" << query_tsv << "\n";
        return 2;
    }

    struct QueryRow {
        std::string id;
        std::string family;
        std::string sequence;
    };
    std::vector<QueryRow> rows;
    while (std::getline(in, line)) {
        const auto fields = split_tab(line);
        if (fields.size() < 3) {
            continue;
        }
        rows.push_back({fields[0], fields[1], fields[2]});
    }
    const int32_t n = std::min<int32_t>(
        std::max(1, limit),
        static_cast<int32_t>(rows.size()));
    if (n <= 0) {
        std::cerr << "no_queries\n";
        return 2;
    }

    int32_t pass_count = 0;
    int32_t family_matches = 0;
    std::cout
        << "id\ttruth_family\tinsert_len\tpred_family\tpred_subfamily"
        << "\tidentity\tquery_coverage\tconfidence\tqc\tintervals\n";
    for (int32_t i = 0; i < n; ++i) {
        const QueryRow& row = rows[static_cast<size_t>(i)];
        const placer::TEAlignmentEvidence evidence =
            classifier.align_insert_sequence(row.sequence);
        if (evidence.pass) {
            ++pass_count;
        }
        if (evidence.best_family == row.family) {
            ++family_matches;
        }
        std::cout
            << row.id << "\t"
            << row.family << "\t"
            << row.sequence.size() << "\t"
            << display_or_na(evidence.best_family) << "\t"
            << display_or_na(evidence.best_subfamily) << "\t"
            << evidence.best_identity << "\t"
            << evidence.best_query_coverage << "\t"
            << evidence.annotation_confidence << "\t"
            << evidence.qc_reason << "\t"
            << evidence.annotation_intervals << "\n";
    }
    std::cout
        << "summary_rows\t" << n << "\n"
        << "summary_pass\t" << pass_count << "\n"
        << "summary_family_matches\t" << family_matches << "\n"
        << "summary_family_accuracy\t"
        << (static_cast<double>(family_matches) / static_cast<double>(n)) << "\n";
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    using namespace placer;

    if (argc > 1) {
        std::string truth_table;
        std::string query_tsv;
        std::string te_fasta;
        int32_t limit = 20;
        int32_t seed = 23;
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            const auto require_value = [&](const std::string& flag) -> std::string {
                if ((i + 1) >= argc) {
                    throw std::runtime_error("missing value for " + flag);
                }
                ++i;
                return argv[i];
            };
            if (arg == "--truth-table") {
                truth_table = require_value(arg);
            } else if (arg == "--query-tsv") {
                query_tsv = require_value(arg);
            } else if (arg == "--te") {
                te_fasta = require_value(arg);
            } else if (arg == "--limit") {
                limit = std::max(1, std::atoi(require_value(arg).c_str()));
            } else if (arg == "--seed") {
                seed = std::atoi(require_value(arg).c_str());
            } else {
                std::cerr
                    << "usage: " << argv[0]
                    << " --truth-table <table.tsv> --te <te.fa>"
                    << " OR --query-tsv <queries.tsv> --te <te.fa>"
                    << " [--limit 20] [--seed 23]\n";
                return 2;
            }
        }
        if (truth_table.empty() || te_fasta.empty()) {
            if (!query_tsv.empty() && !te_fasta.empty()) {
                return run_query_tsv_benchmark(te_fasta, query_tsv, limit);
            }
            std::cerr
                << "usage: " << argv[0]
                << " --truth-table <table.tsv> --te <te.fa>"
                << " OR --query-tsv <queries.tsv> --te <te.fa>"
                << " [--limit 20] [--seed 23]\n";
            return 2;
        }
        return run_truth_table_benchmark(te_fasta, truth_table, limit, seed);
    }

    const std::string fasta_path = make_temp_path("placer_te_alignment_bench");
    std::vector<std::pair<std::string, std::string>> entries;
    entries.reserve(40);
    std::string query;
    const std::string left_noise = build_sequence(20, 4001u);
    const std::string right_noise = build_sequence(16, 4003u);
    for (int family_index = 0; family_index < 20; ++family_index) {
        const std::string family = "Fam" + std::to_string(family_index);
        const std::string seq_a = build_sequence(
            280,
            static_cast<uint32_t>(19 + (family_index * 97)));
        const std::string seq_b = mutate_positions(
            seq_a,
            {
                12,
                29,
                47,
                static_cast<size_t>(65 + family_index),
                static_cast<size_t>(91 + family_index),
                static_cast<size_t>(119 + family_index),
            });
        entries.push_back({family + ":" + family + "-A", seq_a});
        entries.push_back({family + ":" + family + "-B", seq_b});
        if (family_index == 0) {
            query = left_noise + seq_a + right_noise;
        }
    }
    write_te_fasta(fasta_path, entries);

    PipelineConfig cfg;
    cfg.te_fasta_path = fasta_path;
    cfg.te_kmer_size = 9;
    cfg.te_kmer_sizes_csv = "9,11";
    cfg.te_exact_align_topn = 1;
    TEKmerQuickClassifierModule classifier(cfg);

    constexpr int32_t kWarmupIters = 20;
    constexpr int32_t kMeasureIters = 200;
    for (int32_t i = 0; i < kWarmupIters; ++i) {
        (void)classifier.align_insert_sequence(query);
    }

    TEAlignmentEvidence evidence;
    const auto started = std::chrono::steady_clock::now();
    for (int32_t i = 0; i < kMeasureIters; ++i) {
        evidence = classifier.align_insert_sequence(query);
    }
    const double total_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    const double avg_milliseconds =
        (total_seconds * 1000.0) / static_cast<double>(kMeasureIters);

    std::cout
        << "te_entries=" << entries.size() << "\n"
        << "measure_iterations=" << kMeasureIters << "\n"
        << "best_family=" << evidence.best_family << "\n"
        << "best_subfamily=" << (evidence.best_subfamily.empty() ? "NA" : evidence.best_subfamily) << "\n"
        << "best_identity=" << evidence.best_identity << "\n"
        << "best_query_coverage=" << evidence.best_query_coverage << "\n"
        << "last_exact_alignments=" << classifier.last_exact_alignments_ << "\n"
        << "elapsed_total_seconds=" << total_seconds << "\n"
        << "elapsed_avg_milliseconds=" << avg_milliseconds << "\n";

    std::remove(fasta_path.c_str());
    return 0;
}
