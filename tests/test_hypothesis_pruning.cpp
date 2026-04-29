#ifdef NDEBUG
#undef NDEBUG
#endif
#include <algorithm>
#include <cassert>
#include <cstring>
#include <memory>
#include <string>
#include <utility>
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
    const char* sa_tag = nullptr,
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

    if (sa_tag != nullptr) {
        const int aux_ret = bam_aux_append(
            record.get(),
            "SA",
            'Z',
            static_cast<int>(std::strlen(sa_tag) + 1),
            reinterpret_cast<const uint8_t*>(sa_tag));
        assert(aux_ret == 0);
    }

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

placer::Pipeline::HypothesisSummary make_summary(
    size_t original_index,
    int32_t bp_left,
    int32_t bp_right,
    int32_t alt_split_reads,
    int32_t alt_indel_reads,
    int32_t alt_left_clip_reads,
    int32_t alt_right_clip_reads,
    int32_t alt_struct_reads,
    int32_t ref_span_reads,
    int32_t inferred_event_length,
    std::vector<std::string> support_qnames) {
    placer::Pipeline::HypothesisSummary summary;
    summary.original_index = original_index;
    summary.bp_left = bp_left;
    summary.bp_right = bp_right;
    summary.alt_split_reads = alt_split_reads;
    summary.alt_indel_reads = alt_indel_reads;
    summary.alt_left_clip_reads = alt_left_clip_reads;
    summary.alt_right_clip_reads = alt_right_clip_reads;
    summary.alt_struct_reads = alt_struct_reads;
    summary.ref_span_reads = ref_span_reads;
    summary.inferred_event_length = inferred_event_length;
    summary.event_evidence.bp_left = bp_left;
    summary.event_evidence.bp_right = bp_right;
    summary.event_evidence.alt_split_reads = alt_split_reads;
    summary.event_evidence.alt_indel_reads = alt_indel_reads;
    summary.event_evidence.alt_left_clip_reads = alt_left_clip_reads;
    summary.event_evidence.alt_right_clip_reads = alt_right_clip_reads;
    summary.event_evidence.alt_struct_reads = alt_struct_reads;
    summary.event_evidence.ref_span_reads = ref_span_reads;
    summary.event_evidence.support_qnames = support_qnames;
    summary.support_qnames = std::move(support_qnames);
    return summary;
}

}  // namespace

int main() {
    using namespace placer;

    PipelineConfig config;
    Pipeline pipeline(config, nullptr);

    {
        const auto exact = make_summary(
            0,
            100,
            100,
            2,
            0,
            1,
            1,
            4,
            2,
            60,
            {"r1", "r2", "r3", "r4"});
        const auto broad = make_summary(
            1,
            110,
            118,
            2,
            0,
            1,
            1,
            4,
            2,
            62,
            {"r1", "r2", "r3", "r4"});
        const auto competitor = make_summary(
            2,
            220,
            220,
            1,
            0,
            0,
            0,
            1,
            4,
            80,
            {"x1"});

        const auto kept = pipeline.collapse_hypothesis_summaries({exact, broad, competitor});
        assert(kept.size() == 2);
        assert(kept.front().original_index == 0);
        assert(kept.front().bp_left == 100);
        assert(kept.front().bp_right == 100);
        assert(kept.back().original_index == 2);
    }

    {
        const auto indel_record = make_record(
            "indel_read",
            50,
            60,
            {
                bam_cigar_gen(50, BAM_CMATCH),
                bam_cigar_gen(60, BAM_CINS),
                bam_cigar_gen(50, BAM_CMATCH),
            },
            std::string(160, 'A'));
        const auto ref_record = make_record(
            "ref_read",
            40,
            60,
            {bam_cigar_gen(200, BAM_CMATCH)},
            std::string(200, 'G'));

        std::vector<const bam1_t*> local_records = {
            indel_record.get(),
            ref_record.get(),
        };
        std::vector<ReadReferenceSpan> read_spans = {
            make_span(0, 50, 150),
            make_span(0, 40, 240),
        };
        const auto signal_cache = pipeline.build_component_signal_cache(
            local_records,
            read_spans);

        ComponentCall component;
        component.chrom = "chr1";
        component.tid = 0;
        component.anchor_pos = 100;
        BreakpointCandidate indel_bp;
        indel_bp.pos = 100;
        indel_bp.read_id = "indel_read";
        indel_bp.ins_len = 60;
        component.breakpoint_candidates.push_back(indel_bp);

        const auto summary = pipeline.build_hypothesis_summary(
            component,
            signal_cache,
            {},
            0,
            100,
            100);
        assert(summary.original_index == 0);
        assert(summary.bp_left == 100);
        assert(summary.bp_right == 100);
        assert(summary.alt_indel_reads == 1);
        assert(summary.alt_struct_reads == 1);
        assert(summary.ref_span_reads == 1);
        assert(summary.inferred_event_length == 60);
        assert(summary.support_qnames.size() == 1);
        assert(summary.support_qnames.front() == "indel_read");
        assert(summary.event_evidence.bp_left == 100);
        assert(summary.event_evidence.bp_right == 100);
        assert(summary.event_evidence.alt_indel_reads == 1);
        assert(summary.event_evidence.alt_struct_reads == 1);
        assert(summary.event_evidence.ref_span_reads == 1);
        assert(summary.event_evidence.support_qnames.size() == 1);
        assert(summary.event_evidence.support_qnames.front() == "indel_read");
    }

    {
        const auto weak_top = make_summary(
            0,
            100,
            100,
            0,
            0,
            1,
            0,
            1,
            10,
            70,
            {"a1"});
        const auto precise_runner_up = make_summary(
            1,
            150,
            150,
            1,
            0,
            0,
            0,
            1,
            8,
            72,
            {"b1"});
        const auto bilateral_clip = make_summary(
            2,
            200,
            200,
            0,
            0,
            1,
            1,
            2,
            8,
            74,
            {"c1", "c2"});
        const auto weak_bilateral = make_summary(
            3,
            250,
            250,
            0,
            0,
            1,
            1,
            1,
            8,
            74,
            {"d1"});
        const auto unilateral_clip = make_summary(
            4,
            300,
            300,
            0,
            0,
            2,
            0,
            2,
            8,
            74,
            {"e1", "e2"});

        assert(pipeline.should_keep_hypothesis_for_expensive_stage(weak_top, true));
        assert(pipeline.should_keep_hypothesis_for_expensive_stage(precise_runner_up, false));
        assert(pipeline.should_keep_hypothesis_for_expensive_stage(bilateral_clip, false));
        assert(!pipeline.should_keep_hypothesis_for_expensive_stage(weak_bilateral, false));
        assert(!pipeline.should_keep_hypothesis_for_expensive_stage(unilateral_clip, false));
    }

    {
        const auto top_exact = make_summary(
            0,
            100,
            100,
            2,
            0,
            1,
            1,
            4,
            2,
            60,
            {"r1", "r2", "r3", "r4"});
        const auto duplicate_broad = make_summary(
            1,
            110,
            118,
            2,
            0,
            1,
            1,
            4,
            2,
            62,
            {"r1", "r2", "r3", "r4"});
        const auto precise_runner_up = make_summary(
            2,
            180,
            180,
            0,
            1,
            0,
            0,
            1,
            7,
            75,
            {"p1"});
        const auto bilateral_clip = make_summary(
            3,
            260,
            260,
            0,
            0,
            1,
            1,
            2,
            7,
            76,
            {"c1", "c2"});
        const auto unilateral_clip = make_summary(
            4,
            340,
            340,
            0,
            0,
            2,
            0,
            2,
            7,
            76,
            {"u1", "u2"});

        const auto selected = pipeline.select_hypothesis_summaries_for_expensive_stage(
            {top_exact, duplicate_broad, precise_runner_up, bilateral_clip, unilateral_clip});
        assert(selected.size() == 3);
        assert(selected[0].original_index == 0);
        assert(selected[1].original_index == 2);
        assert(selected[2].original_index == 3);
    }

    {
        Pipeline::HypothesisValidatorEvidence lhs;
        lhs.summary.original_index = 1;
        lhs.precise_support = 2;
        lhs.full_context_input_reads = 1;
        lhs.left_anchor_input_reads = 2;
        lhs.right_anchor_input_reads = 2;
        lhs.breakpoint_width = 0;
        lhs.anchor_distance = 3;
        lhs.summary.alt_struct_reads = 4;
        lhs.feasible_for_expensive_stage = true;

        Pipeline::HypothesisValidatorEvidence rhs = lhs;
        rhs.summary.original_index = 0;
        rhs.precise_support = 1;

        assert(pipeline.compare_hypothesis_validator_priority(lhs, rhs));
        assert(!pipeline.compare_hypothesis_validator_priority(rhs, lhs));
    }

    {
        Pipeline::HypothesisValidatorEvidence primary;
        primary.summary.original_index = 0;
        primary.summary.bp_left = 1000;
        primary.summary.bp_right = 1000;
        primary.summary.support_qnames = {"a1", "a2", "a3"};
        primary.precise_support = 2;
        primary.full_context_input_reads = 1;
        primary.left_anchor_input_reads = 2;
        primary.right_anchor_input_reads = 2;
        primary.feasible_for_expensive_stage = true;

        Pipeline::HypothesisValidatorEvidence jitter = primary;
        jitter.summary.original_index = 1;
        jitter.summary.bp_left = 1010;
        jitter.summary.bp_right = 1015;
        jitter.summary.support_qnames = {"a1", "a2", "a3"};

        Pipeline::HypothesisValidatorEvidence challenger = primary;
        challenger.summary.original_index = 2;
        challenger.summary.bp_left = 1200;
        challenger.summary.bp_right = 1200;
        challenger.summary.support_qnames = {"c1", "c2"};
        challenger.precise_support = 1;

        const auto shortlist = pipeline.build_expensive_stage_shortlist(
            {primary, jitter, challenger});
        assert(shortlist.size() == 2);
        assert(shortlist[0].is_primary);
        assert(shortlist[0].validator.summary.original_index == 0);
        assert(!shortlist[1].is_primary);
        assert(shortlist[1].validator.summary.original_index == 2);
    }

    {
        Pipeline::HypothesisValidatorEvidence primary;
        primary.summary.original_index = 0;
        primary.summary.bp_left = 1000;
        primary.summary.bp_right = 1000;
        primary.summary.support_qnames = {"a1", "a2", "a3", "a4"};
        primary.precise_support = 4;
        primary.full_context_input_reads = 2;
        primary.left_anchor_input_reads = 3;
        primary.right_anchor_input_reads = 3;
        primary.feasible_for_expensive_stage = true;

        Pipeline::HypothesisValidatorEvidence first_challenger = primary;
        first_challenger.summary.original_index = 1;
        first_challenger.summary.bp_left = 1400;
        first_challenger.summary.bp_right = 1400;
        first_challenger.summary.support_qnames = {"b1", "b2"};
        first_challenger.precise_support = 2;

        Pipeline::HypothesisValidatorEvidence second_challenger = primary;
        second_challenger.summary.original_index = 2;
        second_challenger.summary.bp_left = 2000;
        second_challenger.summary.bp_right = 2000;
        second_challenger.summary.support_qnames = {"c1", "c2"};
        second_challenger.precise_support = 1;

        const auto shortlist = pipeline.build_expensive_stage_shortlist(
            {primary, first_challenger, second_challenger});
        assert(shortlist.size() == 3);
        assert(shortlist[0].is_primary);
        assert(shortlist[0].validator.summary.original_index == 0);
        assert(!shortlist[1].is_primary);
        assert(shortlist[1].validator.summary.original_index == 1);
        assert(!shortlist[2].is_primary);
        assert(shortlist[2].validator.summary.original_index == 2);
    }

    {
        Pipeline::HypothesisSummary summary;
        summary.original_index = 7;
        summary.bp_left = 2000;
        summary.bp_right = 2010;
        summary.alt_split_reads = 1;
        summary.alt_indel_reads = 2;
        summary.alt_struct_reads = 5;
        summary.support_qnames = {"cached_a", "cached_b"};

        Pipeline::ConsensusInputCounts inputs;
        inputs.full_context_input_reads = 1;
        inputs.partial_context_input_reads = 3;
        inputs.left_anchor_input_reads = 2;
        inputs.right_anchor_input_reads = 2;
        inputs.input_event_reads = 1;

        EventReadEvidence evidence;
        evidence.bp_left = 2000;
        evidence.bp_right = 2010;
        evidence.alt_split_reads = 1;
        evidence.alt_indel_reads = 2;
        evidence.alt_struct_reads = 5;
        evidence.ref_span_reads = 4;
        evidence.support_qnames = {"cached_a", "cached_b"};
        summary.event_evidence = evidence;

        const auto validator = pipeline.collect_hypothesis_validator_evidence(
            summary,
            inputs,
            2004);
        assert(validator.feasible_for_expensive_stage);
        assert(validator.event_evidence.bp_left == 2000);
        assert(validator.event_evidence.bp_right == 2010);
        assert(validator.event_evidence.alt_struct_reads == 5);
        assert(validator.event_evidence.ref_span_reads == 4);
        assert(validator.event_evidence.support_qnames.size() == 2);
        assert(validator.event_evidence.support_qnames.front() == "cached_a");

        Pipeline::HypothesisValidatorEvidence primary = validator;
        primary.summary.original_index = 0;
        primary.summary.bp_left = 1000;
        primary.summary.bp_right = 1000;
        primary.summary.support_qnames = {"p1", "p2", "p3"};
        primary.event_evidence.bp_left = 1000;
        primary.event_evidence.bp_right = 1000;
        primary.precise_support = 4;

        Pipeline::HypothesisValidatorEvidence challenger = validator;
        challenger.summary.original_index = 1;
        challenger.summary.bp_left = 2000;
        challenger.summary.bp_right = 2010;
        challenger.summary.support_qnames = {"cached_a", "cached_b"};
        challenger.event_evidence.bp_left = 2000;
        challenger.event_evidence.bp_right = 2010;
        challenger.precise_support = 2;

        const auto shortlist = pipeline.build_expensive_stage_shortlist(
            {primary, challenger});
        assert(shortlist.size() == 2);
        assert(shortlist[0].validator.event_evidence.bp_left == 1000);
        assert(shortlist[1].validator.event_evidence.bp_left == 2000);
        assert(shortlist[1].validator.event_evidence.support_qnames.front() == "cached_a");
    }

    {
        Pipeline::HypothesisValidatorEvidence infeasible;
        infeasible.summary.original_index = 3;
        infeasible.qc_reason = "VALIDATOR_NO_BILATERAL_ANCHOR";

        const auto shortlist = pipeline.build_expensive_stage_shortlist({infeasible});
        assert(shortlist.empty());
    }

    return 0;
}
