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

    const std::string repeat = std::string(90, 'A') + std::string(50, 'T');
    const std::string gypsy = repeat + build_sequence(90, 111u);
    const std::string copia = repeat + build_sequence(90, 777u);
    const std::string query = repeat + build_sequence(90, 111u);

    const std::vector<std::string> names = {
        "GypsyRepeat#LTR/Gypsy",
        "CopiaRepeat#LTR/Copia",
    };
    const std::vector<std::string> seqs = {
        gypsy,
        copia,
    };
    const std::vector<std::string> revs = {
        reverse_complement(gypsy),
        reverse_complement(copia),
    };

    const TeFamilyCacheBundle cache = build_te_family_cache(names, seqs, revs, 1);
    assert(cache.alignment_index);

    const auto scores = score_family_prefilter(*cache.alignment_index, query, 9, 2);
    assert(scores.size() == 2);
    const auto& best_family =
        cache.alignment_index->groups.at(static_cast<size_t>(scores[0].group_index)).family;
    assert(best_family == "Gypsy");
    assert(scores[0].family_prefilter_score > scores[1].family_prefilter_score);
    assert(scores[0].best_rep_chain_norm > 0.20);
    const uint64_t repeat_key = cache.alignment_index->representative_seed_keys.front();
    assert(cache.alignment_index->kmer_family_counts.at(repeat_key) >= 2);
    assert(cache.alignment_index->kmer_family_weights.at(repeat_key) < 0.6);

    return 0;
}
