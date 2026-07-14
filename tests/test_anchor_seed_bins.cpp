#include <cassert>
#include <string>
#include <vector>

#define private public
#include "pipeline.h"
#undef private

namespace {

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

}  // namespace

int main() {
    using namespace placer;

    PipelineConfig config;
    Pipeline pipeline(config, nullptr);

    std::string reference = build_reference(700);
    const std::string repeated_left = reference.substr(120, 70);
    reference.replace(250, 70, repeated_left);

    const auto bins = pipeline.collect_anchor_seed_bins(
        repeated_left,
        320,
        100,
        reference.substr(100, 260),
        true);

    assert(!bins.empty());
    assert(bins.size() <= 8);
    assert(bins.front().ref_bin_start == 248);
    assert(bins.front().support > 0);

    return 0;
}
