#include <cassert>
#include <cstring>
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

}  // namespace

int main() {
    using namespace placer;

    PipelineConfig config;
    config.event_consensus_poa_min_reads = 1;
    Pipeline pipeline(config, nullptr);

    {
        const auto split_left_record = make_record(
            "split_left",
            100,
            60,
            {
                bam_cigar_gen(30, BAM_CSOFT_CLIP),
                bam_cigar_gen(110, BAM_CMATCH),
            },
            std::string(140, 'A'),
            "chr1,240,+,110M30S,60,1;");
        const auto split_right_record = make_record(
            "split_right",
            0,
            60,
            {
                bam_cigar_gen(110, BAM_CMATCH),
                bam_cigar_gen(30, BAM_CSOFT_CLIP),
            },
            std::string(140, 'C'),
            "chr1,100,+,30S110M,60,1;");
        const auto clip_left_noise = make_record(
            "clip_left_noise",
            70,
            60,
            {
                bam_cigar_gen(30, BAM_CSOFT_CLIP),
                bam_cigar_gen(110, BAM_CMATCH),
            },
            std::string(140, 'G'));
        const auto clip_right_noise = make_record(
            "clip_right_noise",
            50,
            60,
            {
                bam_cigar_gen(110, BAM_CMATCH),
                bam_cigar_gen(30, BAM_CSOFT_CLIP),
            },
            std::string(140, 'T'));
        const auto ref_record = make_record(
            "ref_read",
            40,
            60,
            {bam_cigar_gen(200, BAM_CMATCH)},
            std::string(200, 'G'));

        std::vector<const bam1_t*> local_records = {
            clip_left_noise.get(),
            clip_right_noise.get(),
            ref_record.get(),
            split_left_record.get(),
            split_right_record.get(),
        };
        std::vector<ReadReferenceSpan> read_spans = {
            make_span(0, 70, 180),
            make_span(0, 50, 160),
            make_span(0, 40, 240),
            make_span(0, 100, 210),
            make_span(0, 0, 110),
        };

        ComponentCall component;
        component.chrom = "chr1";
        component.tid = 0;
        component.anchor_pos = 200;
        BreakpointCandidate misleading_bp;
        misleading_bp.pos = 200;
        component.breakpoint_candidates.push_back(misleading_bp);

        const EventReadEvidence evidence = pipeline.collect_event_read_evidence(
            component,
            local_records,
            read_spans,
            {});
        assert(evidence.bp_left == 100);
        assert(evidence.bp_right == 110);
        assert(evidence.alt_split_reads == 2);
        assert(evidence.alt_left_clip_reads == 1);
        assert(evidence.alt_right_clip_reads == 1);
        assert(evidence.alt_struct_reads == 2);
        assert(evidence.ref_span_reads == 1);
        assert(evidence.low_mapq_ref_span_reads == 0);
        assert(evidence.support_qnames.size() == 2);
        assert(evidence.support_qnames[0] == "split_left");
        assert(evidence.support_qnames[1] == "split_right");
        assert(evidence.ref_span_qnames.size() == 1);
        assert(evidence.ref_span_qnames.front() == "ref_read");
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
        const auto clip_left_record = make_record(
            "clip_left",
            95,
            60,
            {
                bam_cigar_gen(30, BAM_CSOFT_CLIP),
                bam_cigar_gen(110, BAM_CMATCH),
            },
            std::string(140, 'C'));
        const auto clip_right_record = make_record(
            "clip_right",
            0,
            60,
            {
                bam_cigar_gen(105, BAM_CMATCH),
                bam_cigar_gen(30, BAM_CSOFT_CLIP),
            },
            std::string(135, 'G'));

        std::vector<const bam1_t*> local_records = {
            indel_record.get(),
            clip_left_record.get(),
            clip_right_record.get(),
        };
        std::vector<ReadReferenceSpan> read_spans = {
            make_span(0, 50, 150),
            make_span(0, 95, 205),
            make_span(0, 0, 105),
        };

        ComponentCall component;
        component.chrom = "chr1";
        component.tid = 0;
        component.anchor_pos = 100;

        const EventReadEvidence evidence = pipeline.collect_event_read_evidence(
            component,
            local_records,
            read_spans,
            {});
        assert(evidence.bp_left == 100);
        assert(evidence.bp_right == 100);
        assert(evidence.alt_indel_reads == 1);
        assert(evidence.alt_left_clip_reads == 1);
        assert(evidence.alt_right_clip_reads == 1);
        assert(evidence.alt_struct_reads == 3);
    }

    {
        const auto indel_record = make_record(
            "indel_support",
            50,
            60,
            {
                bam_cigar_gen(50, BAM_CMATCH),
                bam_cigar_gen(60, BAM_CINS),
                bam_cigar_gen(50, BAM_CMATCH),
            },
            std::string(160, 'A'));
        const auto left_clip_record = make_record(
            "left_clip_support",
            95,
            60,
            {
                bam_cigar_gen(30, BAM_CSOFT_CLIP),
                bam_cigar_gen(110, BAM_CMATCH),
            },
            std::string(140, 'C'));

        std::vector<const bam1_t*> local_records = {
            indel_record.get(),
            left_clip_record.get(),
        };
        std::vector<ReadReferenceSpan> read_spans = {
            make_span(0, 50, 150),
            make_span(0, 95, 205),
        };

        ComponentCall component;
        component.chrom = "chr1";
        component.tid = 0;
        component.anchor_pos = 100;

        const EventReadEvidence evidence = pipeline.collect_event_read_evidence(
            component,
            local_records,
            read_spans,
            {});

        assert(evidence.alt_indel_reads == 1);
        assert(evidence.alt_left_clip_reads == 1);
        assert(evidence.alt_right_clip_reads == 0);
        assert(evidence.alt_struct_reads == 2);
        assert(evidence.support_qnames.size() == 2);
        assert(evidence.support_qnames[0] == "indel_support");
        assert(evidence.support_qnames[1] == "left_clip_support");
    }

    {
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

        const EventReadEvidence evidence = pipeline.collect_event_read_evidence(
            component,
            local_records,
            read_spans,
            {});

        assert(evidence.bp_left == 100);
        assert(evidence.bp_right == 100);
        assert(evidence.alt_indel_reads == 1);
        assert(evidence.alt_left_clip_reads == 1);
        assert(evidence.alt_right_clip_reads == 1);
        assert(evidence.alt_struct_reads == 3);
        assert(evidence.ref_span_reads == 1);
    }

    {
        const auto misleading_split_record = make_record(
            "misleading_split",
            10,
            60,
            {
                bam_cigar_gen(40, BAM_CSOFT_CLIP),
                bam_cigar_gen(90, BAM_CMATCH),
                bam_cigar_gen(40, BAM_CSOFT_CLIP),
            },
            std::string(170, 'G'),
            "chr1,260,+,40S90M40S,60,1;");
        const auto indel_record = make_record(
            "local_indel",
            50,
            60,
            {
                bam_cigar_gen(50, BAM_CMATCH),
                bam_cigar_gen(60, BAM_CINS),
                bam_cigar_gen(50, BAM_CMATCH),
            },
            std::string(160, 'A'));

        InsertionFragment local_fragment;
        local_fragment.read_id = "local_indel";
        local_fragment.read_index = 1;
        local_fragment.source = InsertionFragmentSource::kCigarInsertion;
        local_fragment.ref_junc_pos = 100;
        local_fragment.sequence = std::string(60, 'T');

        std::vector<const bam1_t*> local_records = {
            misleading_split_record.get(),
            indel_record.get(),
        };
        std::vector<ReadReferenceSpan> read_spans = {
            make_span(0, 10, 100),
            make_span(0, 50, 150),
        };

        ComponentCall component;
        component.chrom = "chr1";
        component.tid = 0;
        component.anchor_pos = 100;
        BreakpointCandidate seed_bp;
        seed_bp.pos = 100;
        component.breakpoint_candidates.push_back(seed_bp);

        const EventReadEvidence evidence = pipeline.collect_event_read_evidence(
            component,
            local_records,
            read_spans,
            {local_fragment});
        assert(evidence.bp_left == 100);
        assert(evidence.bp_right == 100);
    }

    {
        std::vector<BamPtr> owned_records;
        std::vector<const bam1_t*> local_records;
        std::vector<ReadReferenceSpan> read_spans;
        std::vector<InsertionFragment> fragments;

        auto add_indel_read = [&](const std::string& qname, int32_t start0, int32_t junction0) {
            owned_records.push_back(make_record(
                qname,
                start0,
                60,
                {
                    bam_cigar_gen(150, BAM_CMATCH),
                    bam_cigar_gen(60, BAM_CINS),
                    bam_cigar_gen(150, BAM_CMATCH),
                },
                std::string(360, 'A')));
            local_records.push_back(owned_records.back().get());
            read_spans.push_back(make_span(0, start0, start0 + 300));

            InsertionFragment fragment;
            fragment.read_id = qname;
            fragment.read_index = static_cast<int32_t>(local_records.size() - 1);
            fragment.source = InsertionFragmentSource::kCigarInsertion;
            fragment.ref_junc_pos = junction0;
            fragment.sequence = std::string(60, 'T');
            fragments.push_back(std::move(fragment));
        };

        for (int i = 0; i < 4; ++i) {
            add_indel_read("near_alt_" + std::to_string(i), 850 + i, 1000);
        }
        for (int i = 0; i < 4; ++i) {
            add_indel_read("far_noise_" + std::to_string(i), 750 + i, 900);
        }
        for (int i = 0; i < 6; ++i) {
            owned_records.push_back(make_record(
                "ref_span_" + std::to_string(i),
                700 + i,
                60,
                {bam_cigar_gen(500, BAM_CMATCH)},
                std::string(500, 'C')));
            local_records.push_back(owned_records.back().get());
            read_spans.push_back(make_span(0, 700 + i, 1200 + i));
        }

        ComponentCall component;
        component.chrom = "chr1";
        component.tid = 0;
        component.anchor_pos = 1000;
        BreakpointCandidate far_bp;
        far_bp.pos = 900;
        component.breakpoint_candidates.push_back(far_bp);
        BreakpointCandidate near_bp;
        near_bp.pos = 1000;
        component.breakpoint_candidates.push_back(near_bp);

        const EventReadEvidence evidence = pipeline.collect_event_read_evidence(
            component,
            local_records,
            read_spans,
            fragments);
        assert(evidence.bp_left == 1000);
        assert(evidence.bp_right == 1000);
        assert(evidence.alt_indel_reads == 4);
        assert(evidence.alt_struct_reads == 4);
        assert(evidence.ref_span_reads == 10);
        assert(evidence.support_qnames.size() == 4);
        assert(evidence.support_qnames.front() == "near_alt_0");

        const GenotypeCall genotype = pipeline.genotype_call(component, evidence);
        assert(genotype.genotype == "0/1");
        assert(genotype.gq >= 20);
    }

    {
        std::vector<BamPtr> owned_records;
        std::vector<const bam1_t*> local_records;
        std::vector<ReadReferenceSpan> read_spans;
        std::vector<InsertionFragment> fragments;

        auto add_near_indel = [&](const std::string& qname, int32_t start0) {
            owned_records.push_back(make_record(
                qname,
                start0,
                60,
                {
                    bam_cigar_gen(150, BAM_CMATCH),
                    bam_cigar_gen(60, BAM_CINS),
                    bam_cigar_gen(150, BAM_CMATCH),
                },
                std::string(360, 'A')));
            local_records.push_back(owned_records.back().get());
            read_spans.push_back(make_span(0, start0, start0 + 300));

            InsertionFragment fragment;
            fragment.read_id = qname;
            fragment.read_index = static_cast<int32_t>(local_records.size() - 1);
            fragment.source = InsertionFragmentSource::kCigarInsertion;
            fragment.ref_junc_pos = 1000;
            fragment.sequence = std::string(60, 'T');
            fragments.push_back(std::move(fragment));
        };

        for (int i = 0; i < 4; ++i) {
            add_near_indel("near_indel_" + std::to_string(i), 850 + i);
        }

        owned_records.push_back(make_record(
            "far_split",
            900,
            60,
            {
                bam_cigar_gen(50, BAM_CSOFT_CLIP),
                bam_cigar_gen(200, BAM_CMATCH),
            },
            std::string(250, 'G'),
            "chr1,1200,+,200M50S,60,1;"));
        local_records.push_back(owned_records.back().get());
        read_spans.push_back(make_span(0, 900, 1100));

        InsertionFragment far_split_fragment;
        far_split_fragment.read_id = "far_split";
        far_split_fragment.read_index = static_cast<int32_t>(local_records.size() - 1);
        far_split_fragment.source = InsertionFragmentSource::kSplitSa;
        far_split_fragment.ref_side = ReferenceSide::kRefLeft;
        far_split_fragment.ref_junc_pos = 900;
        far_split_fragment.sequence = std::string(80, 'G');
        fragments.push_back(std::move(far_split_fragment));

        for (int i = 0; i < 6; ++i) {
            owned_records.push_back(make_record(
                "crossing_ref_" + std::to_string(i),
                700 + i,
                60,
                {bam_cigar_gen(500, BAM_CMATCH)},
                std::string(500, 'C')));
            local_records.push_back(owned_records.back().get());
            read_spans.push_back(make_span(0, 700 + i, 1200 + i));
        }

        ComponentCall component;
        component.chrom = "chr1";
        component.tid = 0;
        component.anchor_pos = 1000;
        BreakpointCandidate far_bp;
        far_bp.pos = 900;
        component.breakpoint_candidates.push_back(far_bp);
        BreakpointCandidate near_bp;
        near_bp.pos = 1000;
        component.breakpoint_candidates.push_back(near_bp);

        const EventReadEvidence evidence = pipeline.collect_event_read_evidence(
            component,
            local_records,
            read_spans,
            fragments);

        assert(evidence.bp_left == 1000);
        assert(evidence.bp_right == 1000);
        assert(evidence.alt_split_reads == 0);
        assert(evidence.alt_indel_reads == 4);
        assert(evidence.alt_struct_reads == 4);
        assert(evidence.ref_span_reads == 7);

        const GenotypeCall genotype = pipeline.genotype_call(component, evidence);
        assert(genotype.genotype == "0/1");
        assert(genotype.gq >= 20);
    }

    {
        const std::string alt_seq =
            std::string(60, 'A') + std::string(20, 'T') + std::string(60, 'C');
        const auto alt_record = make_record(
            "consensus_read",
            95,
            60,
            {
                bam_cigar_gen(30, BAM_CSOFT_CLIP),
                bam_cigar_gen(110, BAM_CMATCH),
            },
            alt_seq,
            "chr1,180,+,60M80S,60,1;");

        InsertionFragment fragment;
        fragment.read_id = "consensus_read";
        fragment.read_index = 9;
        fragment.source = InsertionFragmentSource::kSplitSa;
        fragment.start = 60;
        fragment.length = 20;
        fragment.ref_junc_pos = 100;
        fragment.sequence = std::string(20, 'T');

        EventReadEvidence evidence;
        evidence.bp_left = 100;
        evidence.bp_right = 100;
        evidence.support_qnames = {"consensus_read"};

        const EventConsensus consensus = pipeline.build_event_consensus(
            ComponentCall{},
            {alt_record.get()},
            {fragment},
            evidence);

        assert(consensus.input_event_reads == 1);
        assert(consensus.qc_pass);
        assert(consensus.qc_reason == "PASS_EVENT_CONSENSUS");
        assert(consensus.consensus_len >= 140);
        assert(consensus.consensus_seq == alt_seq);
    }

    {
        const std::string local_seq =
            std::string(80, 'A') + std::string(20, 'T') + std::string(80, 'C') + std::string(80, 'G');
        const std::vector<uint32_t> cigar = {
            static_cast<uint32_t>(bam_cigar_gen(static_cast<int>(local_seq.size()), BAM_CMATCH)),
        };
        const auto local_record = make_record(
            "mixed_fragment_read",
            95,
            60,
            cigar,
            local_seq);

        InsertionFragment local_fragment;
        local_fragment.read_id = "mixed_fragment_read";
        local_fragment.read_index = 1;
        local_fragment.source = InsertionFragmentSource::kCigarInsertion;
        local_fragment.start = 80;
        local_fragment.length = 20;
        local_fragment.ref_junc_pos = 100;
        local_fragment.sequence = std::string(20, 'T');

        InsertionFragment diagnostic_fragment;
        diagnostic_fragment.read_id = "mixed_fragment_read";
        diagnostic_fragment.read_index = 1;
        diagnostic_fragment.source = InsertionFragmentSource::kSplitSa;
        diagnostic_fragment.start = 100;
        diagnostic_fragment.length = 80;
        diagnostic_fragment.ref_junc_pos = -1;
        diagnostic_fragment.sequence = std::string(80, 'C');

        EventReadEvidence evidence;
        evidence.bp_left = 100;
        evidence.bp_right = 100;
        evidence.support_qnames = {"mixed_fragment_read"};

        const EventConsensus consensus = pipeline.build_event_consensus(
            ComponentCall{},
            {local_record.get()},
            {local_fragment, diagnostic_fragment},
            evidence);

        assert(consensus.input_event_reads == 1);
        assert(consensus.qc_pass);
        assert(consensus.qc_reason == "PASS_EVENT_CONSENSUS");
        assert(consensus.consensus_seq == (std::string(80, 'A') + std::string(20, 'T') + std::string(80, 'C')));
    }

    {
        const std::string left_flank(80, 'A');
        const std::string insert_seq(20, 'T');
        const std::string right_flank(80, 'C');

        const auto left_clip_record = make_record(
            "left_clip_consensus",
            100,
            60,
            {
                bam_cigar_gen(80, BAM_CMATCH),
                bam_cigar_gen(20, BAM_CSOFT_CLIP),
            },
            left_flank + insert_seq);
        const auto right_clip_record = make_record(
            "right_clip_consensus",
            100,
            60,
            {
                bam_cigar_gen(20, BAM_CSOFT_CLIP),
                bam_cigar_gen(80, BAM_CMATCH),
            },
            insert_seq + right_flank);

        InsertionFragment left_fragment;
        left_fragment.read_id = "left_clip_consensus";
        left_fragment.read_index = 0;
        left_fragment.source = InsertionFragmentSource::kClipRefRight;
        left_fragment.start = 80;
        left_fragment.length = 20;
        left_fragment.ref_junc_pos = 100;
        left_fragment.sequence = insert_seq;

        InsertionFragment right_fragment;
        right_fragment.read_id = "right_clip_consensus";
        right_fragment.read_index = 1;
        right_fragment.source = InsertionFragmentSource::kClipRefLeft;
        right_fragment.start = 0;
        right_fragment.length = 20;
        right_fragment.ref_junc_pos = 100;
        right_fragment.sequence = insert_seq;

        EventReadEvidence evidence;
        evidence.bp_left = 100;
        evidence.bp_right = 100;
        evidence.support_qnames = {
            "left_clip_consensus",
            "right_clip_consensus",
        };

        const EventConsensus consensus = pipeline.build_event_consensus(
            ComponentCall{},
            {left_clip_record.get(), right_clip_record.get()},
            {left_fragment, right_fragment},
            evidence);

        assert(consensus.input_event_reads == 2);
        assert(consensus.qc_pass);
        assert(consensus.qc_reason == "PASS_EVENT_CONSENSUS");
        assert(consensus.consensus_seq == (left_flank + insert_seq + right_flank));
    }

    {
        const std::string left_flank(80, 'A');
        const std::string insert_seq(20, 'T');
        const std::string right_flank(80, 'C');
        const std::string full_event = left_flank + insert_seq + right_flank;

        const auto full_record = make_record(
            "full_context_read",
            95,
            60,
            {
                bam_cigar_gen(180, BAM_CMATCH),
            },
            full_event);
        const auto clip_record = make_record(
            "clip_rescue_read",
            100,
            60,
            {
                bam_cigar_gen(20, BAM_CSOFT_CLIP),
                bam_cigar_gen(80, BAM_CMATCH),
            },
            insert_seq + right_flank);

        InsertionFragment full_fragment;
        full_fragment.read_id = "full_context_read";
        full_fragment.read_index = 0;
        full_fragment.source = InsertionFragmentSource::kCigarInsertion;
        full_fragment.start = 80;
        full_fragment.length = 20;
        full_fragment.ref_junc_pos = 100;
        full_fragment.sequence = insert_seq;

        InsertionFragment clip_fragment;
        clip_fragment.read_id = "clip_rescue_read";
        clip_fragment.read_index = 1;
        clip_fragment.source = InsertionFragmentSource::kClipRefLeft;
        clip_fragment.start = 0;
        clip_fragment.length = 20;
        clip_fragment.ref_junc_pos = 100;
        clip_fragment.sequence = insert_seq;

        EventReadEvidence evidence;
        evidence.bp_left = 100;
        evidence.bp_right = 100;
        evidence.support_qnames = {
            "full_context_read",
            "clip_rescue_read",
        };

        const EventConsensus consensus = pipeline.build_event_consensus(
            ComponentCall{},
            {full_record.get(), clip_record.get()},
            {full_fragment, clip_fragment},
            evidence);

        assert(consensus.input_event_reads == 1);
        assert(consensus.qc_pass);
        assert(consensus.qc_reason == "PASS_EVENT_CONSENSUS");
        assert(consensus.consensus_seq == full_event);
    }

    {
        const std::string left_flank =
            "AAAACCCCGGGGAAAACCCCGGGGAAAACCCCGGGGAAAACCCCGGGGAAAACCCCGGGGAAAACCCCGGGGAAAA";
        const std::string insert_seq =
            "ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA";
        const std::string right_flank =
            "TTTTGGGGCCCCTTTTGGGGCCCCTTTTGGGGCCCCTTTTGGGGCCCCTTTTGGGGCCCCTTTTGGGGCCCCTTTT";
        const std::string reverse_read_seq = left_flank + insert_seq + right_flank;

        const auto reverse_record = make_record(
            "reverse_full_context",
            95,
            60,
            {
                static_cast<uint32_t>(bam_cigar_gen(static_cast<int>(left_flank.size()), BAM_CMATCH)),
                static_cast<uint32_t>(bam_cigar_gen(static_cast<int>(insert_seq.size()), BAM_CINS)),
                static_cast<uint32_t>(bam_cigar_gen(static_cast<int>(right_flank.size()), BAM_CMATCH)),
            },
            reverse_read_seq,
            nullptr,
            BAM_FREVERSE);

        InsertionFragment reverse_fragment;
        reverse_fragment.read_id = "reverse_full_context";
        reverse_fragment.read_index = 0;
        reverse_fragment.source = InsertionFragmentSource::kCigarInsertion;
        reverse_fragment.start = static_cast<int32_t>(left_flank.size());
        reverse_fragment.length = static_cast<int32_t>(insert_seq.size());
        reverse_fragment.ref_junc_pos = 100;
        reverse_fragment.sequence = insert_seq;
        reverse_fragment.is_reverse = true;

        EventReadEvidence evidence;
        evidence.bp_left = 100;
        evidence.bp_right = 100;
        evidence.support_qnames = {"reverse_full_context"};

        const EventConsensus consensus = pipeline.build_event_consensus(
            ComponentCall{},
            {reverse_record.get()},
            {reverse_fragment},
            evidence);
        assert(consensus.input_event_reads == 1);
        assert(consensus.qc_pass);
        assert(consensus.qc_reason == "PASS_EVENT_CONSENSUS");
        assert(consensus.consensus_seq == (left_flank + insert_seq + right_flank));
    }

    return 0;
}
