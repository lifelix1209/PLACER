#include "pipeline.h"

#include <cassert>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>
#include <unistd.h>

namespace {

std::string make_temp_path(const std::string& stem) {
    const long long pid = static_cast<long long>(::getpid());
    return "/tmp/" + stem + "_" + std::to_string(pid) + ".fa";
}

void write_te_fasta(const std::string& path) {
    std::ofstream out(path);
    assert(out.is_open());
    out << ">TEA\n";
    out << "ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA\n";
    out << ">TEB\n";
    out << "TTGGAACCTTGGAACCTTGGAACCTTGGAACCTTGGAACC\n";
}

}  // namespace

int main() {
    using namespace placer;

    const std::string fasta_path = make_temp_path("placer_te_quick_test");
    write_te_fasta(fasta_path);

    {
        PipelineConfig config;
        config.te_fasta_path = fasta_path;
        config.ins_fragment_hits_tsv_path.clear();
        config.te_kmer_size = 13;
        config.te_kmer_sizes_csv = "9,13";
        config.te_low_kmer_rescue_enable = false;

        TEKmerQuickClassifierModule classifier(config);
        assert(classifier.is_enabled());

        InsertionFragment frag;
        frag.fragment_id = "frag_multik";
        frag.sequence = "ACGTTGCAAC";  // len=10, cannot hit k=13, should hit k=9
        std::vector<InsertionFragment> fragments = {frag};

        const auto hits = classifier.classify(fragments);
        assert(hits.size() == 1);
        assert(hits[0].te_name == "TEA");
        assert(hits[0].multik_support > 0.0);
        assert(hits[0].kmer_support > 0.0);
    }

    {
        PipelineConfig config;
        config.te_fasta_path = fasta_path;
        config.ins_fragment_hits_tsv_path.clear();
        config.te_kmer_size = 13;
        config.te_kmer_sizes_csv = "13";
        config.te_median_identity_min = 0.95;  // force low-k rescue trigger
        config.te_low_kmer_rescue_enable = true;
        config.te_low_kmer_rescue_topn = 2;
        config.te_low_kmer_rescue_min_frag_len = 20;
        config.te_low_kmer_rescue_identity_min = 0.80;
        config.te_low_kmer_rescue_margin_max = 0.20;

        TEKmerQuickClassifierModule classifier(config);
        assert(classifier.is_enabled());

        InsertionFragment frag;
        frag.fragment_id = "frag_rescue";
        frag.sequence = "ACGTTGCAACGTTACAACGTTGCAACGTTG";  // one mismatch from TEA
        std::vector<InsertionFragment> fragments = {frag};

        const auto hits = classifier.classify(fragments);
        assert(hits.size() == 1);
        assert(hits[0].te_name == "TEA");
        assert(hits[0].rescue_used);
        assert(hits[0].kmer_support >= hits[0].multik_support);
    }

    std::remove(fasta_path.c_str());
    return 0;
}
