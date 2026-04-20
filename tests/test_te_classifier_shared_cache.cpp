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
    const long long pid = static_cast<long long>(::getpid());
    return "/tmp/" + stem + "_" + std::to_string(pid) + ".fa";
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

    const std::string fasta_path = make_temp_path("placer_te_classifier_shared_cache");
    write_te_fasta(
        fasta_path,
        {
            {"SubGypsyCache#LTR/Gypsy", build_sequence(160, 11u)},
            {"SubCopiaCache#LTR/Copia", build_sequence(160, 97u)},
        });

    PipelineConfig cfg;
    cfg.te_fasta_path = fasta_path;
    cfg.te_kmer_size = 9;
    cfg.te_kmer_sizes_csv = "9,11";

    TEKmerQuickClassifierModule first(cfg);
    TEKmerQuickClassifierModule second(cfg);

    assert(first.is_enabled());
    assert(second.is_enabled());
    assert(first.primary_index_);
    assert(second.primary_index_);
    assert(first.alignment_shortlist_db_);
    assert(second.alignment_shortlist_db_);
    assert(first.primary_index_.get() == second.primary_index_.get());
    assert(first.alignment_shortlist_db_.get() == second.alignment_shortlist_db_.get());
    assert(first.indices_.size() == second.indices_.size());
    for (size_t i = 0; i < first.indices_.size(); ++i) {
        assert(first.indices_[i].get() == second.indices_[i].get());
    }
    assert(first.te_names_.get() == second.te_names_.get());
    assert(first.te_sequences_.get() == second.te_sequences_.get());
    assert(first.te_reverse_complement_sequences_.get() ==
           second.te_reverse_complement_sequences_.get());

    std::remove(fasta_path.c_str());
    return 0;
}
