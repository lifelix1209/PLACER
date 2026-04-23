#ifdef NDEBUG
#undef NDEBUG
#endif
#include "te_family_alignment.h"

#include <algorithm>
#include <cassert>
#include <string>
#include <vector>

namespace {

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

std::string mutate_every_n(std::string seq, size_t start, size_t step) {
    for (size_t i = start; i < seq.size(); i += step) {
        switch (seq[i]) {
            case 'A': seq[i] = 'C'; break;
            case 'C': seq[i] = 'G'; break;
            case 'G': seq[i] = 'T'; break;
            case 'T': seq[i] = 'A'; break;
            default: seq[i] = 'A'; break;
        }
    }
    return seq;
}

std::string reverse_complement(std::string seq) {
    std::reverse(seq.begin(), seq.end());
    for (char& c : seq) {
        switch (c) {
            case 'A': c = 'T'; break;
            case 'C': c = 'G'; break;
            case 'G': c = 'C'; break;
            case 'T': c = 'A'; break;
            default: c = 'N'; break;
        }
    }
    return seq;
}

}  // namespace

int main() {
    using namespace placer;

    const std::string outlier = build_sequence(180, 501u);
    const std::string medoid = build_sequence(180, 777u);
    const std::string close_a = mutate_every_n(medoid, 5, 23);
    const std::string close_b = mutate_every_n(medoid, 11, 29);
    const std::string copia = build_sequence(180, 977u);

    const std::vector<std::string> names = {
        "GypsyOutlier#LTR/Gypsy",
        "GypsyMedoid#LTR/Gypsy",
        "GypsyCloseA#LTR/Gypsy",
        "GypsyCloseB#LTR/Gypsy",
        "CopiaRef#LTR/Copia",
    };
    const std::vector<std::string> seqs = {
        outlier,
        medoid,
        close_a,
        close_b,
        copia,
    };
    const std::vector<std::string> revs = {
        reverse_complement(outlier),
        reverse_complement(medoid),
        reverse_complement(close_a),
        reverse_complement(close_b),
        reverse_complement(copia),
    };

    const TeFamilyCacheBundle cache = build_te_family_cache(names, seqs, revs, 2);
    assert(cache.alignment_index);
    const auto& group = cache.alignment_index->groups.at(
        cache.alignment_index->family_to_group.at("GYPSY"));
    assert(group.representatives.size() == 2);
    assert(group.representatives.front().exact_name == "GypsyMedoid");

    const std::string noisy_gypsy = mutate_every_n(medoid, 3, 13);
    const auto scores = score_family_prefilter(*cache.alignment_index, noisy_gypsy, 9, 4);
    assert(scores.size() >= 2);
    const auto& best = scores.front();
    const auto& runner = scores.at(1);
    const std::string best_family =
        cache.alignment_index->groups.at(static_cast<size_t>(best.group_index)).family;
    const std::string runner_family =
        cache.alignment_index->groups.at(static_cast<size_t>(runner.group_index)).family;
    assert(best_family == "Gypsy");
    assert(runner_family == "Copia");
    assert(best.best_rep_chain_norm > runner.best_rep_chain_norm);
    assert(best.chained_query_coverage > runner.chained_query_coverage);

    return 0;
}
