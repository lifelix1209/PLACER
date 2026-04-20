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
    uint16_t flag,
    int32_t pos0,
    int32_t mapq,
    const std::vector<uint32_t>& cigar,
    const std::string& seq,
    const char* sa_tag = nullptr) {
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
            static_cast<int>(std::char_traits<char>::length(sa_tag) + 1),
            reinterpret_cast<const uint8_t*>(sa_tag));
        assert(aux_ret == 0);
    }

    return record;
}

std::string make_sequence(size_t len, char base) {
    return std::string(len, base);
}

}  // namespace

int main() {
    using namespace placer;

    {
        const auto read = make_record(
            "single_ins",
            0,
            1000,
            60,
            {
                bam_cigar_gen(70, BAM_CMATCH),
                bam_cigar_gen(90, BAM_CINS),
                bam_cigar_gen(70, BAM_CMATCH),
            },
            make_sequence(230, 'A'));

        std::vector<const bam1_t*> records = {read.get()};

        DbscanComponentModule module;
        const auto components = module.build(records, "chr1", 0);

        assert(components.size() == 1);
        const auto& component = components.front();
        assert(component.chrom == "chr1");
        assert(component.tid == 0);
        assert(component.anchor_pos == 1070);
        assert(component.read_indices.size() == 1);
        assert(component.insertion_read_indices.size() == 1);
        assert(component.breakpoint_candidates.size() == 1);
        assert(component.breakpoint_candidates.front().ins_len == 90);
        assert(component.breakpoint_candidates.front().clip_len == 0);
    }

    {
        const auto split_read = make_record(
            "split_like",
            0,
            1000,
            60,
            {
                bam_cigar_gen(90, BAM_CMATCH),
                bam_cigar_gen(60, BAM_CSOFT_CLIP),
            },
            make_sequence(150, 'C'),
            "chr1,1091,+,60S90M,60,1;");

        std::vector<const bam1_t*> records = {split_read.get()};
        DbscanComponentModule module;
        const auto components = module.build(records, "chr1", 0);

        assert(components.size() == 1);
        assert(components.front().split_sa_read_indices.size() == 1);
        assert(components.front().soft_clip_read_indices.empty());
        assert(components.front().breakpoint_candidates.size() == 1);
        assert(
            (components.front().breakpoint_candidates.front().class_mask &
             kCandidateSplitSaSupplementary) != 0);
    }

    {
        const auto clip_read = make_record(
            "clip_hint",
            0,
            2000,
            60,
            {
                bam_cigar_gen(120, BAM_CSOFT_CLIP),
                bam_cigar_gen(180, BAM_CMATCH),
            },
            make_sequence(300, 'G'));

        std::vector<const bam1_t*> records = {clip_read.get()};
        DbscanComponentModule module;
        const auto components = module.build(records, "chr1", 0);

        assert(components.size() == 1);
        assert(components.front().soft_clip_read_indices.size() == 1);
        assert(components.front().split_sa_read_indices.empty());
        assert(components.front().breakpoint_candidates.front().clip_len == 120);
    }

    {
        const auto read_a = make_record(
            "cluster_a_1",
            0,
            1000,
            60,
            {
                bam_cigar_gen(80, BAM_CMATCH),
                bam_cigar_gen(90, BAM_CINS),
                bam_cigar_gen(80, BAM_CMATCH),
            },
            make_sequence(250, 'A'));
        const auto read_b = make_record(
            "cluster_a_2",
            0,
            1010,
            60,
            {
                bam_cigar_gen(75, BAM_CMATCH),
                bam_cigar_gen(92, BAM_CINS),
                bam_cigar_gen(75, BAM_CMATCH),
            },
            make_sequence(242, 'C'));
        const auto read_c = make_record(
            "cluster_b_1",
            0,
            1005,
            60,
            {
                bam_cigar_gen(80, BAM_CMATCH),
                bam_cigar_gen(400, BAM_CINS),
                bam_cigar_gen(80, BAM_CMATCH),
            },
            make_sequence(560, 'G'));
        const auto read_d = make_record(
            "cluster_b_2",
            0,
            1015,
            60,
            {
                bam_cigar_gen(75, BAM_CMATCH),
                bam_cigar_gen(405, BAM_CINS),
                bam_cigar_gen(75, BAM_CMATCH),
            },
            make_sequence(555, 'T'));

        std::vector<const bam1_t*> records = {
            read_a.get(),
            read_b.get(),
            read_c.get(),
            read_d.get(),
        };
        DbscanComponentModule module;
        const auto components = module.build(records, "chr1", 0);

        assert(components.size() == 2);
        assert(components[0].read_indices.size() == 2);
        assert(components[1].read_indices.size() == 2);
        assert(components[0].anchor_pos != components[1].anchor_pos ||
               components[0].breakpoint_candidates.front().ins_len !=
                   components[1].breakpoint_candidates.front().ins_len);
        const int32_t len_a = components[0].breakpoint_candidates.front().ins_len;
        const int32_t len_b = components[1].breakpoint_candidates.front().ins_len;
        assert(std::abs(len_a - len_b) >= 200);
    }

    return 0;
}
