#ifndef PLACER_GATE1_MODULE_H
#define PLACER_GATE1_MODULE_H

#include "pipeline.h"

#include <cstdint>

namespace placer {

struct Gate1SignalConfig {
    // Hard filters
    int32_t min_seq_len = 50;

    // Signal criteria
    int32_t long_soft_clip_min = 100;

    // Fallback background retention when no signal is present.
    int32_t background_mapq_min = 20;

    // Fuse 1: require at least one solid reference anchor.
    int32_t min_anchor_match_bases = 200;

    // Fuse 2: for long soft-clip signals, require clip-adjacent flank anchor.
    int32_t min_clip_flank_match_bases = 120;

    // Fuse 3: if NM exists, reject very poor alignments.
    double max_nm_rate = 0.20;
};

class SignalFirstGate1Module final : public IGate1Module {
public:
    explicit SignalFirstGate1Module(Gate1SignalConfig config = Gate1SignalConfig{});

    bool pass_preliminary(const ReadView& read) const override;

private:
    Gate1SignalConfig config_;
};

}  // namespace placer

#endif  // PLACER_GATE1_MODULE_H
