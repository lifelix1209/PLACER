#include <cassert>
#include <memory>
#include <string>
#include <vector>

#include <htslib/sam.h>

#define private public
#include "pipeline.h"
#undef private

namespace {

using BamPtr = std::unique_ptr<bam1_t, decltype(&bam_destroy1)>;

BamPtr make_record(
    const std::string& qname,
    int32_t pos0,
    const std::vector<uint32_t>& cigar,
    const std::string& seq) {
    BamPtr record(bam_init1(), &bam_destroy1);
    assert(record != nullptr);

    const std::string qual(seq.size(), 'I');
    const int ret = bam_set1(
        record.get(),
        static_cast<size_t>(qname.size() + 1),
        qname.c_str(),
        0,
        0,
        pos0,
        60,
        static_cast<size_t>(cigar.size()),
        cigar.data(),
        -1,
        -1,
        0,
        static_cast<size_t>(seq.size()),
        seq.c_str(),
        qual.c_str(),
        0);
    assert(ret >= 0);
    return record;
}

placer::InsertionFragment make_fragment(
    placer::InsertionFragmentSource source,
    const std::string& read_id,
    int32_t start,
    int32_t length,
    int32_t read_len,
    int32_t ref_junc_pos) {
    placer::InsertionFragment fragment;
    fragment.source = source;
    fragment.read_id = read_id;
    fragment.start = start;
    fragment.length = length;
    fragment.read_len = read_len;
    fragment.ref_junc_pos = ref_junc_pos;
    return fragment;
}

}  // namespace

int main() {
    using namespace placer;

    PipelineConfig config;
    Pipeline pipeline(config, nullptr);

    const auto full_record = make_record(
        "full_read",
        1000,
        {
            bam_cigar_gen(80, BAM_CMATCH),
            bam_cigar_gen(40, BAM_CINS),
            bam_cigar_gen(80, BAM_CMATCH),
        },
        std::string(200, 'A'));
    const auto left_clip_record = make_record(
        "left_clip_read",
        1000,
        {
            bam_cigar_gen(80, BAM_CSOFT_CLIP),
            bam_cigar_gen(120, BAM_CMATCH),
        },
        std::string(200, 'C'));
    const auto right_clip_record = make_record(
        "right_clip_read",
        920,
        {
            bam_cigar_gen(120, BAM_CMATCH),
            bam_cigar_gen(80, BAM_CSOFT_CLIP),
        },
        std::string(200, 'G'));

    std::vector<const bam1_t*> local_records = {
        full_record.get(),
        left_clip_record.get(),
        right_clip_record.get(),
    };

    std::vector<InsertionFragment> fragments = {
        make_fragment(InsertionFragmentSource::kCigarInsertion, "full_read", 80, 40, 200, 1000),
        make_fragment(InsertionFragmentSource::kClipRefLeft, "left_clip_read", 0, 80, 200, 1005),
        make_fragment(InsertionFragmentSource::kClipRefRight, "right_clip_read", 120, 80, 200, 995),
    };

    EventReadEvidence event_evidence;
    event_evidence.bp_left = 995;
    event_evidence.bp_right = 1005;
    event_evidence.support_qnames = {"full_read", "left_clip_read", "right_clip_read"};

    const auto summary = pipeline.collect_event_consensus_inputs(
        local_records,
        fragments,
        event_evidence);
    const auto counts = pipeline.collect_event_consensus_input_counts(
        local_records,
        fragments,
        event_evidence);

    assert(counts.full_context_input_reads == summary.full_context_input_reads);
    assert(counts.partial_context_input_reads == summary.partial_context_input_reads);
    assert(counts.left_anchor_input_reads == summary.left_anchor_input_reads);
    assert(counts.right_anchor_input_reads == summary.right_anchor_input_reads);
    assert(counts.input_event_reads ==
        static_cast<int32_t>(
            summary.full_event_by_qname.empty()
                ? summary.partial_event_by_qname.size()
                : summary.full_event_by_qname.size()));

    return 0;
}
