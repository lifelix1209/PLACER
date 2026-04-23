#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <cstdio>
#include <fstream>
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
    assert(out.is_open());
    for (const auto& entry : entries) {
        out << ">" << entry.first << "\n";
        out << entry.second << "\n";
    }
}

}  // namespace

int main() {
    using namespace placer;

    const std::string fasta_path = make_temp_path("placer_te_exact_budget");
    const std::string gypsy_a = build_sequence(180, 17u);
    const std::string gypsy_b = mutate_positions(
        gypsy_a,
        {12, 24, 36, 48, 60, 72});
    const std::string gypsy_c = mutate_positions(
        gypsy_a,
        {15, 30, 45, 75, 90, 105, 120});
    write_te_fasta(
        fasta_path,
        {
            {"GypsyA#LTR/Gypsy", gypsy_a},
            {"GypsyB#LTR/Gypsy", gypsy_b},
            {"GypsyC#LTR/Gypsy", gypsy_c},
        });

    PipelineConfig cfg;
    cfg.te_fasta_path = fasta_path;
    cfg.te_kmer_size = 9;
    cfg.te_kmer_sizes_csv = "9,11";
    cfg.te_exact_align_topn = 1;
    TEKmerQuickClassifierModule classifier(cfg);

    const TEAlignmentEvidence evidence = classifier.align_insert_sequence(gypsy_a);
    assert(evidence.pass);
    assert(evidence.best_family == "Gypsy");
    assert(classifier.last_exact_alignments_ <= 1);

    std::remove(fasta_path.c_str());
    return 0;
}
