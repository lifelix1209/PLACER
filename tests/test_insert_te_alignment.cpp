#include <algorithm>
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

char complement_base(char c) {
    switch (c) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return 'N';
    }
}

std::string reverse_complement(std::string seq) {
    std::reverse(seq.begin(), seq.end());
    for (char& c : seq) {
        c = complement_base(c);
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

placer::TEAlignmentEvidence run_alignment(
    const std::string& te_fasta,
    const std::string& insert_seq) {
    placer::PipelineConfig cfg;
    cfg.te_fasta_path = te_fasta;
    cfg.te_kmer_size = 9;
    cfg.te_kmer_sizes_csv = "9,11";
    placer::Pipeline pipeline(cfg, nullptr);

    placer::EventSegmentation segmentation;
    segmentation.insert_seq = insert_seq;
    segmentation.pass = true;
    segmentation.qc_reason = "PASS_EVENT_SEGMENTATION";
    return pipeline.align_insert_seq_to_te(segmentation);
}

}  // namespace

int main() {
    using namespace placer;

    const std::string fasta_path = make_temp_path("placer_insert_te_alignment");

    {
        const std::string gypsy = build_sequence(120, 11u);
        const std::string copia = build_sequence(120, 97u);
        write_te_fasta(
            fasta_path,
            {
                {"SubGypsyA#LTR/Gypsy", gypsy},
                {"SubCopiaA#LTR/Copia", copia},
            });

        const TEAlignmentEvidence evidence = run_alignment(fasta_path, gypsy);
        assert(evidence.pass);
        assert(evidence.qc_reason == "PASS_INSERT_TE_ALIGNMENT");
        assert(evidence.best_family == "Gypsy");
        assert(evidence.best_subfamily == "SubGypsyA");
        assert(evidence.best_identity >= 0.99);
        assert(evidence.best_query_coverage >= 0.99);
        assert(evidence.cross_family_margin >= 0.09);
    }

    {
        const std::string gypsy = build_sequence(120, 23u);
        const std::string copia = build_sequence(120, 131u);
        write_te_fasta(
            fasta_path,
            {
                {"SubGypsyRC#LTR/Gypsy", gypsy},
                {"SubCopiaRC#LTR/Copia", copia},
            });

        const TEAlignmentEvidence evidence = run_alignment(
            fasta_path,
            reverse_complement(gypsy));
        assert(evidence.pass);
        assert(evidence.qc_reason == "PASS_INSERT_TE_ALIGNMENT");
        assert(evidence.best_family == "Gypsy");
        assert(evidence.best_subfamily == "SubGypsyRC");
        assert(evidence.best_identity >= 0.99);
        assert(evidence.best_query_coverage >= 0.99);
        assert(evidence.cross_family_margin >= 0.09);
    }

    {
        const std::string gypsy = build_sequence(120, 37u);
        const std::string copia = mutate_positions(
            gypsy,
            {12, 27, 43, 58, 74, 89});
        write_te_fasta(
            fasta_path,
            {
                {"SubGypsyB#LTR/Gypsy", gypsy},
                {"SubCopiaB#LTR/Copia", copia},
            });

        const TEAlignmentEvidence evidence = run_alignment(fasta_path, gypsy);
        assert(!evidence.pass);
        assert(evidence.best_family == "Gypsy");
        assert(evidence.best_subfamily == "SubGypsyB");
        assert(evidence.second_family == "Copia");
        assert(evidence.qc_reason == "TE_ALIGNMENT_CROSS_FAMILY_AMBIGUOUS");
        assert(evidence.cross_family_margin < 0.10);
    }

    {
        const std::string deu = build_sequence(160, 211u);
        const std::string deu_same_family = mutate_positions(
            deu,
            {11, 24, 39, 52, 67, 81});
        const std::string rex = build_sequence(160, 313u);
        write_te_fasta(
            fasta_path,
            {
                {"5S-Deu-L2:5S-Deu-L2-3", deu},
                {"5S-Deu-L2:5S-Deu-L2-4", deu_same_family},
                {"Rex-Babar:Rex-Babar-18", rex},
            });

        const TEAlignmentEvidence evidence = run_alignment(fasta_path, deu);
        assert(evidence.pass);
        assert(evidence.qc_reason == "PASS_INSERT_TE_ALIGNMENT");
        assert(evidence.best_family == "5S-Deu-L2");
        assert(evidence.best_subfamily == "5S-Deu-L2-3");
        assert(evidence.second_family == "Rex-Babar");
        assert(evidence.cross_family_margin >= 0.09);
    }

    {
        const std::string deu = build_sequence(160, 509u);
        const std::string dna_like = mutate_positions(
            deu,
            {10, 25, 40, 55, 70, 85, 100, 115, 130, 145});
        write_te_fasta(
            fasta_path,
            {
                {"5S-Deu-L2:5S-Deu-L2-3", deu},
                {"DNA:NA-DNA-24", dna_like},
            });

        const TEAlignmentEvidence evidence = run_alignment(fasta_path, deu);
        assert(evidence.pass);
        assert(evidence.qc_reason == "PASS_INSERT_TE_ALIGNMENT");
        assert(evidence.best_family == "5S-Deu-L2");
        assert(evidence.best_subfamily == "5S-Deu-L2-3");
        assert(evidence.second_family == "DNA");
        assert(evidence.cross_family_margin >= 0.09);
    }

    {
        const std::string gypsy = build_sequence(120, 101u);
        write_te_fasta(
            fasta_path,
            {
                {"SubGypsyC#LTR/Gypsy", gypsy},
            });

        const TEAlignmentEvidence evidence = run_alignment(
            fasta_path,
            gypsy.substr(0, 60));
        assert(!evidence.pass);
        assert(evidence.qc_reason == "INSERT_SEQ_TOO_SHORT");
    }

    {
        const std::string gypsy = build_sequence(160, 701u);
        const std::string copia = build_sequence(160, 907u);
        std::vector<size_t> block_mutations;
        for (size_t pos = 30; pos < 70; ++pos) {
            block_mutations.push_back(pos);
        }
        const std::string divergent_insert = mutate_positions(gypsy, block_mutations);
        write_te_fasta(
            fasta_path,
            {
                {"SubGypsyD#LTR/Gypsy", gypsy},
                {"SubCopiaD#LTR/Copia", copia},
            });

        const TEAlignmentEvidence evidence = run_alignment(fasta_path, divergent_insert);
        assert(evidence.pass);
        assert(evidence.qc_reason == "PASS_INSERT_TE_ALIGNMENT");
        assert(evidence.best_family == "Gypsy");
        assert(evidence.best_subfamily == "SubGypsyD");
        assert(evidence.second_family == "Copia");
        assert(evidence.best_identity >= 0.80);
        assert(evidence.best_query_coverage >= 0.80);
    }

    std::remove(fasta_path.c_str());
    return 0;
}
