#ifdef NDEBUG
#undef NDEBUG
#endif
#include <algorithm>
#include <cassert>
#include <memory>
#include <string>
#include <vector>

#include <htslib/sam.h>

#include "pipeline.h"

namespace {

using BamPtr = std::unique_ptr<bam1_t, decltype(&bam_destroy1)>;

BamPtr make_record(
    const std::string& qname,
    int32_t pos0,
    int32_t mapq,
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

bool has_insertion_fragment_at(
    const std::vector<placer::InsertionFragment>& fragments,
    int32_t ref_junc_pos,
    int32_t length) {
    return std::any_of(
        fragments.begin(),
        fragments.end(),
        [&](const placer::InsertionFragment& fragment) {
            return fragment.source == placer::InsertionFragmentSource::kCigarInsertion &&
                   fragment.ref_junc_pos == ref_junc_pos &&
                   fragment.length == length;
        });
}

}  // namespace

int main() {
    using namespace placer;

    PipelineConfig config;
    config.min_long_ins_for_seq_extract = 50;
    CigarInsertionFragmentModule module(config);

    ComponentCall component;
    component.chrom = "4";
    component.tid = 0;
    component.anchor_pos = 21752342;
    component.read_indices = {0};
    component.insertion_read_indices = {0};

    BreakpointCandidate truth_candidate;
    truth_candidate.chrom = "4";
    truth_candidate.pos = 21752342;
    truth_candidate.ins_len = 387;
    truth_candidate.read_id = "multi_insert_read";
    component.breakpoint_candidates.push_back(truth_candidate);

    const int32_t pos0 = 21734874;
    const std::vector<uint32_t> cigar = {
        bam_cigar_gen(4837, BAM_CMATCH),
        bam_cigar_gen(861, BAM_CINS),
        bam_cigar_gen(10910, BAM_CMATCH),
        bam_cigar_gen(1689, BAM_CINS),
        bam_cigar_gen(1721, BAM_CMATCH),
        bam_cigar_gen(387, BAM_CINS),
        bam_cigar_gen(120, BAM_CMATCH),
    };
    const std::string seq(4837 + 861 + 10910 + 1689 + 1721 + 387 + 120, 'A');
    const auto read = make_record("multi_insert_read", pos0, 60, cigar, seq);
    const std::vector<const bam1_t*> records = {read.get()};

    const auto fragments = module.extract(component, records);

    assert(has_insertion_fragment_at(fragments, 21750621, 1689));
    assert(has_insertion_fragment_at(fragments, 21739711, 861));
    assert(has_insertion_fragment_at(fragments, 21752342, 387));

    {
        const auto low_mapq_read = make_record("low_mapq_insert_read", pos0, 59, cigar, seq);
        const std::vector<const bam1_t*> low_mapq_records = {low_mapq_read.get()};
        const auto low_mapq_fragments = module.extract(component, low_mapq_records);

        assert(low_mapq_fragments.empty());
    }

    return 0;
}
