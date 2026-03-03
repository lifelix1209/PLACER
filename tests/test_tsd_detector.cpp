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
    out << "AAAAATGACCCCCCTGACGGGGG\n";
    out << ">chr2\n";
    out << "AAAAATGACGGGGG\n";
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

        // chr1: AAAAA [TGAC] CCCCC [TGAC] GGGGG
        // left_bp=9, right_bp=13 should recover DUP len=4 seq=TGAC
        const TsdDetection dup = detector.detect("chr1", 9, 13);
        assert(dup.type == "DUP");
        assert(dup.length == 4);
        assert(dup.sequence == "TGAC");
        assert(dup.bg_p <= 0.20);
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

    std::remove(fasta_path.c_str());
    std::remove((fasta_path + ".fai").c_str());
    return 0;
}
