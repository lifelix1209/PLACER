#ifdef NDEBUG
#undef NDEBUG
#endif
#include "pipeline.h"

#include <cassert>
#include <cstdio>
#include <fstream>
#include <string>
#include <unistd.h>

namespace {

std::string make_temp_path(const std::string& stem, const std::string& ext) {
    const long long pid = static_cast<long long>(::getpid());
    return "/tmp/" + stem + "_" + std::to_string(pid) + ext;
}

void write_fasta(const std::string& path) {
    std::ofstream out(path);
    assert(out.is_open());
    out << ">chr1\n";
    out << "AAAAATGACCCCCTGACGGGGG\n";
    out << ">chr2\n";
    out << "AAAAATGACGGGGG\n";
    out << ">chrN\n";
    out << "AAACNNNGGG\n";
    out << ">chrGap\n";
    out << std::string(100, 'A') << std::string(100, 'N') << std::string(100, 'A') << "\n";
}

}  // namespace

int main() {
    using namespace placer;

    const std::string fasta_path = make_temp_path("placer_tsd_ref", ".fa");
    write_fasta(fasta_path);

    {
        PipelineConfig cfg;
        cfg.reference_fasta_path = fasta_path;
        cfg.tsd_enable = true;
        cfg.tsd_min_len = 3;
        cfg.tsd_max_len = 8;
        cfg.tsd_flank_window = 40;
        cfg.tsd_bg_p_max = 0.20;

        TSDDetector detector(cfg);
        assert(detector.is_enabled());

        // chr1: AAAAA [TGAC] CCCC [TGAC] GGGGG
        // left_bp=9, right_bp=13 should recover DUP len=4 seq=TGAC
        const TsdDetection dup = detector.detect("chr1", 9, 13);
        assert(dup.type == "DUP");
        assert(dup.length == 4);
        assert(dup.sequence == "TGAC");
        assert(dup.bg_p <= 0.20);

        const TsdDetection dup_alias = detector.detect("1", 9, 13);
        assert(dup_alias.type == "DUP");
        assert(dup_alias.length == 4);
        assert(dup_alias.sequence == "TGAC");
        assert(dup_alias.bg_p <= 0.20);
    }

    {
        PipelineConfig cfg;
        cfg.reference_fasta_path = fasta_path;
        cfg.tsd_enable = true;
        cfg.tsd_min_len = 3;
        cfg.tsd_max_len = 8;
        cfg.tsd_flank_window = 40;
        cfg.tsd_bg_p_max = 0.20;

        TSDDetector detector(cfg);
        assert(detector.is_enabled());

        // chr2: AAAAA [TGAC] GGGGG, treat as target-site deletion candidate.
        const TsdDetection del = detector.detect("chr2", 5, 9);
        assert(del.type == "DEL");
        assert(del.length == 4);
        assert(del.sequence == "TGAC");
        assert(del.bg_p <= 0.20);
    }

    {
        PipelineConfig cfg;
        cfg.reference_fasta_path = fasta_path;

        TSDDetector detector(cfg);

        assert(!detector.reference_position_is_poly_n("chrN", 2));
        assert(!detector.reference_position_is_poly_n("chrN", 3));
        assert(detector.reference_position_is_poly_n("chrN", 4));
        assert(detector.reference_position_is_poly_n("chrN", 5));
        assert(detector.reference_position_is_poly_n("chrN", 6));
        assert(!detector.reference_position_is_poly_n("chrN", 7));
        assert(!detector.reference_position_is_poly_n("chrN", -1));

        PipelineResult result;
        FinalCall keep;
        keep.chrom = "chrN";
        keep.tid = 2;
        keep.pos = 2;
        keep.bp_left = 2;
        keep.bp_right = 2;
        keep.support_reads = 6;
        keep.gq = 99;
        keep.te_name = "TE_KEEP";

        FinalCall reject;
        reject.chrom = "chrN";
        reject.tid = 2;
        reject.pos = 4;
        reject.bp_left = 4;
        reject.bp_right = 4;
        reject.support_reads = 12;
        reject.gq = 99;
        reject.te_name = "TE_REJECT";

        result.final_calls.push_back(reject);
        result.final_calls.push_back(keep);

        finalize_final_calls(result, detector);

        assert(result.final_calls.size() == 1);
        assert(result.final_calls.front().te_name == "TE_KEEP");
        assert(result.final_pass_calls == 1);
    }

    {
        PipelineConfig cfg;
        cfg.reference_fasta_path = fasta_path;

        TSDDetector detector(cfg);

        assert(!detector.reference_position_is_poly_n("chrGap", 90));

        PipelineResult result;
        FinalCall reject_context;
        reject_context.chrom = "chrGap";
        reject_context.tid = 3;
        reject_context.pos = 90;
        reject_context.bp_left = 90;
        reject_context.bp_right = 90;
        reject_context.support_reads = 12;
        reject_context.gq = 99;
        reject_context.te_name = "TE_REJECT_N_CONTEXT";

        result.final_calls.push_back(reject_context);

        finalize_final_calls(result, detector);

        assert(result.final_calls.empty());
        assert(result.final_pass_calls == 0);
    }

    std::remove(fasta_path.c_str());
    std::remove((fasta_path + ".fai").c_str());
    return 0;
}
