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
    int32_t mapq,
    const std::vector<uint32_t>& cigar,
    const std::string& seq,
    uint16_t flag = 0) {
    BamPtr record(bam_init1(), &bam_destroy1);
    assert(record != nullptr);

    const std::string qual(seq.size(), 'I');
    const int ret = bam_set1(
        record.get(),
        static_cast<size_t>(qname.size() + 1),
        qname.c_str(),
        flag,
        0,
        pos0,
        mapq,
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

placer::ReadReferenceSpan make_span(int32_t tid, int32_t start, int32_t end) {
    placer::ReadReferenceSpan span;
    span.valid = true;
    span.tid = tid;
    span.start = start;
    span.end = end;
    return span;
}

}  // namespace

int main() {
    using namespace placer;

    PipelineConfig config;
    Pipeline pipeline(config, nullptr);

    const auto indel_record = make_record(
        "offset_indel_support",
        50,
        60,
        {
            bam_cigar_gen(50, BAM_CMATCH),
            bam_cigar_gen(60, BAM_CINS),
            bam_cigar_gen(50, BAM_CMATCH),
        },
        std::string(160, 'A'));
    const auto left_clip_record = make_record(
        "offset_left_clip",
        175,
        60,
        {
            bam_cigar_gen(30, BAM_CSOFT_CLIP),
            bam_cigar_gen(110, BAM_CMATCH),
        },
        std::string(140, 'C'));
    const auto right_clip_record = make_record(
        "offset_right_clip",
        70,
        60,
        {
            bam_cigar_gen(105, BAM_CMATCH),
            bam_cigar_gen(30, BAM_CSOFT_CLIP),
        },
        std::string(135, 'G'));
    const auto ref_record = make_record(
        "offset_ref",
        0,
        60,
        {
            bam_cigar_gen(250, BAM_CMATCH),
        },
        std::string(250, 'T'));

    std::vector<const bam1_t*> local_records = {
        indel_record.get(),
        left_clip_record.get(),
        right_clip_record.get(),
        ref_record.get(),
    };
    std::vector<ReadReferenceSpan> read_spans = {
        make_span(0, 50, 150),
        make_span(0, 175, 285),
        make_span(0, 70, 175),
        make_span(0, 0, 250),
    };

    ComponentCall component;
    component.chrom = "chr1";
    component.tid = 0;
    component.anchor_pos = 175;
    BreakpointCandidate indel_bp;
    indel_bp.pos = 100;
    BreakpointCandidate left_clip_bp;
    left_clip_bp.pos = 175;
    BreakpointCandidate right_clip_bp;
    right_clip_bp.pos = 175;
    component.breakpoint_candidates = {
        indel_bp,
        left_clip_bp,
        right_clip_bp,
    };

    const auto cache = pipeline.build_component_signal_cache(
        local_records,
        read_spans);
    const EventReadEvidence old_evidence = pipeline.collect_event_read_evidence_for_bounds(
        component,
        local_records,
        read_spans,
        {},
        100,
        100,
        100,
        100);
    const EventReadEvidence new_evidence = pipeline.collect_event_read_evidence_for_bounds_cached(
        component,
        cache,
        {},
        100,
        100,
        100,
        100);

    assert(old_evidence.bp_left == new_evidence.bp_left);
    assert(old_evidence.bp_right == new_evidence.bp_right);
    assert(old_evidence.alt_split_reads == new_evidence.alt_split_reads);
    assert(old_evidence.alt_indel_reads == new_evidence.alt_indel_reads);
    assert(old_evidence.alt_left_clip_reads == new_evidence.alt_left_clip_reads);
    assert(old_evidence.alt_right_clip_reads == new_evidence.alt_right_clip_reads);
    assert(old_evidence.alt_struct_reads == new_evidence.alt_struct_reads);
    assert(old_evidence.ref_span_reads == new_evidence.ref_span_reads);
    assert(old_evidence.low_mapq_ref_span_reads == new_evidence.low_mapq_ref_span_reads);
    assert(old_evidence.support_qnames == new_evidence.support_qnames);
    assert(old_evidence.ref_span_qnames == new_evidence.ref_span_qnames);

    return 0;
}
