#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
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

}  // namespace

int main() {
    using namespace placer;

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
