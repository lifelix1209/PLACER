#include <algorithm>
#include <cassert>
#include <string>
#include <vector>

#define private public
#include "pipeline.h"
#undef private

namespace {

placer::InsertionFragment make_split_fragment(
    const std::string& read_id,
    int32_t ref_junc_pos,
    placer::ReferenceSide ref_side) {
    placer::InsertionFragment frag;
    frag.chrom = "chr9";
    frag.anchor_pos = 28720475;
    frag.read_id = read_id;
    frag.source = placer::InsertionFragmentSource::kSplitSa;
    frag.ref_side = ref_side;
    frag.ref_junc_pos = ref_junc_pos;
    frag.start = 40;
    frag.length = 120;
    frag.read_len = 400;
    frag.sequence = std::string(120, 'A');
    frag.split_sa_reliable = true;
    return frag;
}

bool has_hypothesis(
    const std::vector<placer::Pipeline::BreakpointHypothesis>& hypotheses,
    int32_t left,
    int32_t right) {
    return std::any_of(
        hypotheses.begin(),
        hypotheses.end(),
        [&](const auto& h) {
            return h.valid && h.left == left && h.right == right;
        });
}

bool has_hypothesis_in_window(
    const std::vector<placer::Pipeline::BreakpointHypothesis>& hypotheses,
    int32_t left_lo,
    int32_t left_hi,
    int32_t right_lo,
    int32_t right_hi) {
    return std::any_of(
        hypotheses.begin(),
        hypotheses.end(),
        [&](const auto& h) {
            return h.valid &&
                   h.left >= left_lo &&
                   h.left <= left_hi &&
                   h.right >= right_lo &&
                   h.right <= right_hi;
        });
}

}  // namespace

int main() {
    using namespace placer;

    PipelineConfig config;
    Pipeline pipeline(config, nullptr);

    {
        Pipeline::BreakpointHypothesis clip_raw;
        clip_raw.valid = true;
        clip_raw.left = 1000;
        clip_raw.right = 1020;
        clip_raw.center = 1010;
        clip_raw.support = 39;
        clip_raw.priority = 5;

        Pipeline::BreakpointHypothesis clip_fragment = clip_raw;
        clip_fragment.left = 1003;
        clip_fragment.right = 1029;
        clip_fragment.center = 1016;
        clip_fragment.support = 32;
        clip_fragment.priority = 3;

        Pipeline::BreakpointHypothesis distant_clip = clip_raw;
        distant_clip.left = 1300;
        distant_clip.right = 1400;
        distant_clip.center = 1350;
        distant_clip.support = 19;
        distant_clip.priority = 3;

        Pipeline::BreakpointHypothesis precise = clip_raw;
        precise.left = 1100;
        precise.right = 1100;
        precise.center = 1100;
        precise.support = 2;
        precise.priority = 1;

        const auto selected = pipeline.select_diverse_breakpoint_hypotheses(
            {clip_raw, clip_fragment, distant_clip, precise},
            3);
        assert(selected.size() == 3);
        assert(has_hypothesis_in_window(selected, 1000, 1003, 1020, 1029));
        assert(has_hypothesis(selected, 1300, 1400));
        assert(has_hypothesis(selected, 1100, 1100));
    }

    ComponentCall component;
    component.chrom = "chr9";
    component.tid = 0;
    component.anchor_pos = 28720475;

    std::vector<const bam1_t*> local_records;
    std::vector<InsertionFragment> fragments;

    fragments.push_back(make_split_fragment("truth_l1", 28719854, ReferenceSide::kRefLeft));
    fragments.push_back(make_split_fragment("truth_l2", 28719854, ReferenceSide::kRefLeft));
    fragments.push_back(make_split_fragment("truth_l3", 28719913, ReferenceSide::kRefLeft));
    fragments.push_back(make_split_fragment("truth_l4", 28719913, ReferenceSide::kRefLeft));
    fragments.push_back(make_split_fragment("truth_r1", 28719918, ReferenceSide::kRefRight));

    fragments.push_back(make_split_fragment("fail_l1", 28720684, ReferenceSide::kRefLeft));
    fragments.push_back(make_split_fragment("fail_l2", 28720684, ReferenceSide::kRefLeft));
    fragments.push_back(make_split_fragment("fail_r1", 28720494, ReferenceSide::kRefRight));
    fragments.push_back(make_split_fragment("fail_r2", 28720494, ReferenceSide::kRefRight));

    const auto hypotheses = pipeline.enumerate_breakpoint_hypotheses(
        component,
        local_records,
        fragments,
        28718727,
        28722368,
        6);

    assert(!hypotheses.empty());
    assert(has_hypothesis_in_window(hypotheses, 28719854, 28719913, 28719918, 28719918));
    assert(has_hypothesis(hypotheses, 28720494, 28720684));
    return 0;
}
