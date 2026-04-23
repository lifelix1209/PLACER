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
        const auto read_with_remote_softclip = make_record(
            "read_with_remote_softclip",
            0,
            1000,
            60,
            {
                bam_cigar_gen(60, BAM_CSOFT_CLIP),
                bam_cigar_gen(500, BAM_CMATCH),
                bam_cigar_gen(200, BAM_CINS),
                bam_cigar_gen(100, BAM_CMATCH),
            },
            make_sequence(860, 'A'));
        const auto local_indel_read = make_record(
            "local_indel_read",
            0,
            1450,
            60,
            {
                bam_cigar_gen(50, BAM_CMATCH),
                bam_cigar_gen(200, BAM_CINS),
                bam_cigar_gen(100, BAM_CMATCH),
            },
            make_sequence(350, 'C'));

        const std::vector<const bam1_t*> records = {
            read_with_remote_softclip.get(),
            local_indel_read.get(),
        };

        LinearBinComponentModule module;
        const auto components = module.build(records, "chr1", 0, 900, 1800);
        assert(components.size() == 1);
        const auto& component = components.front();
        assert(component.soft_clip_read_indices.empty());
        assert(component.insertion_read_indices.size() == 2);
        assert(component.breakpoint_candidates.size() == 2);
        for (const auto& bp : component.breakpoint_candidates) {
            assert(bp.pos >= component.bin_start);
            assert(bp.pos <= component.bin_end);
            assert(bp.ins_len >= 200);
            assert(bp.clip_len == 0);
        }
        assert(component.anchor_pos >= 1500);
    }

    {
        const auto supplementary_only = make_record(
            "supplementary_only",
            BAM_FSUPPLEMENTARY,
            2000,
            60,
            {
                bam_cigar_gen(60, BAM_CSOFT_CLIP),
                bam_cigar_gen(220, BAM_CMATCH),
            },
            make_sequence(280, 'G'),
            "chr1,2400,+,220M60S,60,1;");

        const std::vector<const bam1_t*> records = {
            supplementary_only.get(),
        };

        LinearBinComponentModule module;
        const auto components = module.build(records, "chr1", 0, 1800, 2600);
        assert(components.empty());
    }

    {
        const auto ambiguous_softclip_vs_indel = make_record(
            "ambiguous_softclip_vs_indel",
            0,
            1000,
            60,
            {
                bam_cigar_gen(500, BAM_CSOFT_CLIP),
                bam_cigar_gen(500, BAM_CMATCH),
                bam_cigar_gen(400, BAM_CINS),
                bam_cigar_gen(100, BAM_CMATCH),
            },
            make_sequence(1500, 'T'));
        const auto softclip_support = make_record(
            "softclip_support",
            0,
            1000,
            60,
            {
                bam_cigar_gen(500, BAM_CSOFT_CLIP),
                bam_cigar_gen(300, BAM_CMATCH),
            },
            make_sequence(800, 'A'));
        const auto indel_support = make_record(
            "indel_support",
            0,
            1450,
            60,
            {
                bam_cigar_gen(50, BAM_CMATCH),
                bam_cigar_gen(400, BAM_CINS),
                bam_cigar_gen(100, BAM_CMATCH),
            },
            make_sequence(550, 'C'));

        const std::vector<const bam1_t*> records = {
            ambiguous_softclip_vs_indel.get(),
            softclip_support.get(),
            indel_support.get(),
        };

        LinearBinComponentModule module;
        const auto components = module.build(records, "chr1", 0, 900, 1800);

        assert(components.size() == 2);

        const ComponentCall* softclip_component = nullptr;
        const ComponentCall* indel_component = nullptr;
        for (const auto& component : components) {
            if (component.anchor_pos < 1200) {
                softclip_component = &component;
            } else {
                indel_component = &component;
            }
        }

        assert(softclip_component != nullptr);
        assert(indel_component != nullptr);

        assert(softclip_component->read_indices.size() == 1);
        assert(softclip_component->soft_clip_read_indices.size() == 1);
        assert(softclip_component->insertion_read_indices.empty());

        assert(indel_component->anchor_pos >= 1450);
        assert(indel_component->anchor_pos <= 1550);
        assert(indel_component->read_indices.size() == 2);
        assert(indel_component->soft_clip_read_indices.empty());
        assert(indel_component->insertion_read_indices.size() == 2);
        assert(indel_component->breakpoint_candidates.size() == 2);
        for (const auto& bp : indel_component->breakpoint_candidates) {
            assert(bp.pos >= 1500);
            assert(bp.pos <= 1500);
            assert(bp.ins_len >= 400);
            assert(bp.clip_len == 0);
        }
    }

    return 0;
}
