#include <cassert>
#include <cstdio>
#include <fstream>
#include <memory>
#include <string>
#include <unistd.h>

#define private public
#include "pipeline.h"
#undef private

namespace {

std::string make_temp_path(const std::string& stem, const std::string& ext) {
    const long long pid = static_cast<long long>(::getpid());
    return "/tmp/" + stem + "_" + std::to_string(pid) + ext;
}

std::string build_reference(size_t len) {
    static const char kBases[] = {'A', 'C', 'G', 'T'};
    uint32_t state = 17u;
    std::string out;
    out.reserve(len);
    for (size_t i = 0; i < len; ++i) {
        state = (state * 1103515245u) + 12345u;
        out.push_back(kBases[(state >> 16) & 3u]);
    }
    return out;
}

void write_fasta(const std::string& path, const std::string& seq) {
    std::ofstream out(path);
    assert(out.is_open());
    out << ">chr1\n";
    out << seq << '\n';
}

placer::EventSegmentation run_segmentation(
    const std::string& fasta_path,
    const std::string& consensus_seq,
    int32_t bp_left,
    int32_t bp_right) {
    placer::PipelineConfig cfg;
    cfg.reference_fasta_path = fasta_path;
    placer::Pipeline pipeline(cfg, nullptr);

    placer::ComponentCall component;
    component.chrom = "chr1";
    component.tid = 0;

    placer::EventReadEvidence evidence;
    evidence.bp_left = bp_left;
    evidence.bp_right = bp_right;

    placer::EventConsensus consensus;
    consensus.consensus_seq = consensus_seq;
    consensus.consensus_len = static_cast<int32_t>(consensus_seq.size());
    consensus.input_event_reads = 3;
    consensus.qc_pass = true;
    consensus.qc_reason = "PASS_EVENT_CONSENSUS";

    return pipeline.segment_event_consensus(component, evidence, consensus);
}

}  // namespace

int main() {
    using namespace placer;

    const std::string fasta_path = make_temp_path("placer_event_segmentation", ".fa");

    {
        const std::string reference = build_reference(600);
        write_fasta(fasta_path, reference);

        const int32_t bp = 220;
        const std::string left = reference.substr(150, 70);
        const std::string insert = "TTGGAACCTTGGAACCTTGGAACC";
        const std::string right = reference.substr(220, 70);
        const EventSegmentation segmentation = run_segmentation(
            fasta_path,
            left + insert + right,
            bp,
            bp);

        assert(segmentation.pass);
        assert(segmentation.qc_reason == "PASS_EVENT_SEGMENTATION");
        assert(segmentation.left_flank_seq == left);
        assert(segmentation.insert_seq == insert);
        assert(segmentation.right_flank_seq == right);
        assert(segmentation.left_ref_start == 150);
        assert(segmentation.left_ref_end == 220);
        assert(segmentation.right_ref_start == 220);
        assert(segmentation.right_ref_end == 290);
        assert(segmentation.left_flank_align_len == 70);
        assert(segmentation.right_flank_align_len == 70);
        assert(segmentation.left_flank_identity >= 0.99);
        assert(segmentation.right_flank_identity >= 0.99);
    }

    {
        std::string reference = build_reference(700);
        const std::string repeated_left = reference.substr(120, 70);
        reference.replace(250, 70, repeated_left);
        write_fasta(fasta_path, reference);

        const int32_t bp = 320;
        const std::string insert = "AACCGGTTAACCGGTTAACCGGTT";
        const std::string right = reference.substr(320, 70);
        const EventSegmentation segmentation = run_segmentation(
            fasta_path,
            repeated_left + insert + right,
            bp,
            bp);

        assert(!segmentation.pass);
        assert(segmentation.qc_reason == "LEFT_FLANK_AMBIGUOUS");
    }

    std::remove(fasta_path.c_str());
    std::remove((fasta_path + ".fai").c_str());
    return 0;
}
