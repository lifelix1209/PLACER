#include <cassert>
#include <memory>
#include <string>
#include <vector>

#include <htslib/sam.h>

#include "gate1_module.h"

namespace {

using BamPtr = std::unique_ptr<bam1_t, decltype(&bam_destroy1)>;

BamPtr make_record(
    const std::string& qname,
    uint16_t flag,
    int32_t pos0,
    int32_t mapq,
    const std::vector<uint32_t>& cigar,
    const std::string& seq,
    int32_t nm = -1) {
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

    if (nm >= 0) {
        const int aux_ret = bam_aux_append(
            record.get(),
            "NM",
            'i',
            sizeof(nm),
            reinterpret_cast<const uint8_t*>(&nm));
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

    SignalFirstGate1Module module;

    {
        const auto low_mapq_long_insertion = make_record(
            "low_mapq_long_insertion",
            0,
            1000,
            5,
            {
                bam_cigar_gen(250, BAM_CMATCH),
                bam_cigar_gen(400, BAM_CINS),
                bam_cigar_gen(250, BAM_CMATCH),
            },
            make_sequence(900, 'A'));
        ReadView read(low_mapq_long_insertion.get());
        assert(module.pass_preliminary(read));
    }

    {
        const auto low_mapq_background = make_record(
            "low_mapq_background",
            0,
            1000,
            5,
            {
                bam_cigar_gen(600, BAM_CMATCH),
            },
            make_sequence(600, 'C'));
        ReadView read(low_mapq_background.get());
        assert(!module.pass_preliminary(read));
    }

    {
        const auto long_insertion_with_bad_clip_and_nm = make_record(
            "long_insertion_with_bad_clip_and_nm",
            0,
            1000,
            1,
            {
                bam_cigar_gen(400, BAM_CSOFT_CLIP),
                bam_cigar_gen(50, BAM_CMATCH),
                bam_cigar_gen(400, BAM_CINS),
                bam_cigar_gen(300, BAM_CMATCH),
                bam_cigar_gen(40, BAM_CSOFT_CLIP),
            },
            make_sequence(1190, 'G'),
            220);
        ReadView read(long_insertion_with_bad_clip_and_nm.get());
        assert(module.pass_preliminary(read));
    }

    return 0;
}
