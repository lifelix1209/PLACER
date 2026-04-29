#include <algorithm>
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <cstdio>
#include <fstream>
#include <string>
#include <unordered_map>
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
    cfg.te_family_margin_min = 0.05;
    cfg.te_subfamily_margin_min = 0.04;
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
        const std::string gypsy = build_sequence(140, 17u);
        const std::string copia = build_sequence(140, 97u);
        write_te_fasta(
            fasta_path,
            {
                {"SubGypsyCached#LTR/Gypsy", gypsy},
                {"SubCopiaCached#LTR/Copia", copia},
            });

        PipelineConfig cfg;
        cfg.te_fasta_path = fasta_path;
        cfg.te_kmer_size = 9;
        cfg.te_kmer_sizes_csv = "9,11";
        cfg.te_family_margin_min = 0.05;
        cfg.te_subfamily_margin_min = 0.04;
        Pipeline pipeline(cfg, nullptr);

        std::unordered_map<std::string, TEAlignmentEvidence> cache;

        EventSegmentation failed_segmentation;
        failed_segmentation.insert_seq = gypsy;
        failed_segmentation.pass = false;
        failed_segmentation.qc_reason = "NO_EVENT_SEGMENTATION";
        const TEAlignmentEvidence failed = pipeline.align_insert_seq_to_te_cached(
            failed_segmentation,
            cache);
        assert(!failed.pass);
        assert(failed.qc_reason == "NO_EVENT_SEGMENTATION_FOR_TE_ALIGNMENT");
        assert(cache.empty());

        EventSegmentation segmentation;
        segmentation.insert_seq = gypsy;
        segmentation.pass = true;
        segmentation.qc_reason = "PASS_EVENT_SEGMENTATION";

        const TEAlignmentEvidence first = pipeline.align_insert_seq_to_te_cached(
            segmentation,
            cache);
        assert(first.pass);
        assert(first.best_family == "Gypsy");
        assert(first.best_subfamily == "SubGypsyCached");
        assert(cache.size() == 1);
        assert(pipeline.te_classifier_module_.last_exact_alignments_ > 0);

        pipeline.te_classifier_module_.last_exact_alignments_ = -123;
        const TEAlignmentEvidence second = pipeline.align_insert_seq_to_te_cached(
            segmentation,
            cache);
        assert(second.pass);
        assert(second.qc_reason == first.qc_reason);
        assert(second.best_family == first.best_family);
        assert(second.best_subfamily == first.best_subfamily);
        assert(second.best_identity == first.best_identity);
        assert(cache.size() == 1);
        assert(pipeline.te_classifier_module_.last_exact_alignments_ == -123);
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
        const std::string gypsy = build_sequence(180, 1201u);
        const std::string copia = build_sequence(180, 1601u);
        const std::string contaminated =
            build_sequence(22, 901u) + gypsy + build_sequence(18, 903u);
        write_te_fasta(
            fasta_path,
            {
                {"SubGypsyOccupancy#LTR/Gypsy", gypsy},
                {"SubCopiaOccupancy#LTR/Copia", copia},
            });

        const TEAlignmentEvidence evidence = run_alignment(fasta_path, contaminated);
        assert(evidence.pass);
        assert(evidence.best_family == "Gypsy");
        assert(evidence.best_subfamily == "SubGypsyOccupancy");
        assert(evidence.best_identity >= 0.95);
        assert(evidence.best_query_coverage > 0.78);
        assert(evidence.best_query_coverage < 0.95);
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
        assert(evidence.pass);
        assert(evidence.best_family == "UNKNOWN");
        assert(evidence.best_subfamily == "UNKNOWN");
        assert(evidence.qc_reason == "PASS_INSERT_TE_ALIGNMENT_UNKNOWN");
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
        const std::string gypsy_a = build_sequence(160, 601u);
        const std::string gypsy_b = mutate_positions(
            gypsy_a,
            {10, 30, 50, 70, 90, 110, 130, 150});
        const std::string query = mutate_positions(
            gypsy_a,
            {10, 50, 90, 130});
        const std::string copia = build_sequence(160, 907u);
        write_te_fasta(
            fasta_path,
            {
                {"SubGypsyFamilyOnlyA#LTR/Gypsy", gypsy_a},
                {"SubGypsyFamilyOnlyB#LTR/Gypsy", gypsy_b},
                {"SubCopiaFamilyOnly#LTR/Copia", copia},
            });

        const TEAlignmentEvidence evidence = run_alignment(fasta_path, query);
        assert(evidence.pass);
        assert(evidence.best_family == "Gypsy");
        assert(evidence.best_subfamily.empty());
        assert(evidence.qc_reason == "PASS_INSERT_TE_ALIGNMENT_FAMILY_ONLY");
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
        assert(evidence.cross_family_margin > 0.05);
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
