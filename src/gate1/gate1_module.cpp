#include "gate1_module.h"

#include <algorithm>
#include <cstdint>

#include <htslib/sam.h>

namespace placer {
namespace {

struct CigarSummary {
    int32_t total_match_bases = 0;
    int32_t max_match_block = 0;

    int32_t max_soft_clip = 0;
    int32_t leading_soft_clip = 0;
    int32_t trailing_soft_clip = 0;

    int32_t right_anchor_after_leading = 0;
    int32_t left_anchor_before_trailing = 0;
};

bool is_match_like(int op) {
    return op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF;
}

int find_first_non_hard_clip(const uint32_t* cigar, int32_t n_cigar) {
    for (int32_t i = 0; i < n_cigar; ++i) {
        if (bam_cigar_op(cigar[i]) != BAM_CHARD_CLIP) {
            return i;
        }
    }
    return -1;
}

int find_last_non_hard_clip(const uint32_t* cigar, int32_t n_cigar) {
    for (int32_t i = n_cigar - 1; i >= 0; --i) {
        if (bam_cigar_op(cigar[i]) != BAM_CHARD_CLIP) {
            return i;
        }
    }
    return -1;
}

CigarSummary summarize_cigar(const ReadView& read) {
    CigarSummary s;

    const uint32_t* cigar = read.cigar();
    const int32_t n_cigar = read.n_cigar();
    if (!cigar || n_cigar <= 0) {
        return s;
    }

    int32_t current_match_block = 0;

    for (int32_t i = 0; i < n_cigar; ++i) {
        const int op = bam_cigar_op(cigar[i]);
        const int32_t len = static_cast<int32_t>(bam_cigar_oplen(cigar[i]));

        if (is_match_like(op)) {
            s.total_match_bases += len;
            current_match_block += len;
            s.max_match_block = std::max(s.max_match_block, current_match_block);
        } else {
            current_match_block = 0;
        }

        if (op == BAM_CSOFT_CLIP) {
            s.max_soft_clip = std::max(s.max_soft_clip, len);
        }
    }

    const int first = find_first_non_hard_clip(cigar, n_cigar);
    const int last = find_last_non_hard_clip(cigar, n_cigar);
    if (first < 0 || last < 0 || first > last) {
        return s;
    }

    if (bam_cigar_op(cigar[first]) == BAM_CSOFT_CLIP) {
        s.leading_soft_clip = static_cast<int32_t>(bam_cigar_oplen(cigar[first]));
        int32_t flank = 0;
        for (int i = first + 1; i <= last; ++i) {
            const int op = bam_cigar_op(cigar[i]);
            if (!is_match_like(op)) {
                break;
            }
            flank += static_cast<int32_t>(bam_cigar_oplen(cigar[i]));
        }
        s.right_anchor_after_leading = flank;
    }

    if (bam_cigar_op(cigar[last]) == BAM_CSOFT_CLIP) {
        s.trailing_soft_clip = static_cast<int32_t>(bam_cigar_oplen(cigar[last]));
        int32_t flank = 0;
        for (int i = last - 1; i >= first; --i) {
            const int op = bam_cigar_op(cigar[i]);
            if (!is_match_like(op)) {
                break;
            }
            flank += static_cast<int32_t>(bam_cigar_oplen(cigar[i]));
        }
        s.left_anchor_before_trailing = flank;
    }

    return s;
}

}  // namespace

SignalFirstGate1Module::SignalFirstGate1Module(Gate1SignalConfig config)
    : config_(config) {}

bool SignalFirstGate1Module::pass_preliminary(const ReadView& read) const {
    const uint16_t flag = read.flag();

    // Hard drop list requested by design.
    if ((flag & BAM_FUNMAP) != 0 || (flag & BAM_FSECONDARY) != 0) {
        return false;
    }

    if (read.seq_len() < config_.min_seq_len) {
        return false;
    }

    const CigarSummary cigar = summarize_cigar(read);

    const bool has_supplementary = (flag & BAM_FSUPPLEMENTARY) != 0;
    const bool has_sa = read.has_sa_tag();
    const bool has_long_soft_clip = cigar.max_soft_clip >= config_.long_soft_clip_min;

    const bool has_signal = has_supplementary || has_sa || has_long_soft_clip;

    if (!has_signal) {
        return read.mapq() > config_.background_mapq_min;
    }

    // Fuse 1: at least one solid anchor on reference.
    if (cigar.max_match_block < config_.min_anchor_match_bases) {
        return false;
    }

    // Fuse 2: for long soft-clips, require acceptable flank anchors.
    if (has_long_soft_clip) {
        const bool leading_long = cigar.leading_soft_clip >= config_.long_soft_clip_min;
        const bool trailing_long = cigar.trailing_soft_clip >= config_.long_soft_clip_min;

        if (leading_long && cigar.right_anchor_after_leading < config_.min_clip_flank_match_bases) {
            return false;
        }
        if (trailing_long && cigar.left_anchor_before_trailing < config_.min_clip_flank_match_bases) {
            return false;
        }
    }

    // Fuse 3: if NM exists, reject overly noisy alignments.
    int64_t nm = -1;
    if (read.get_int_tag("NM", nm) && cigar.total_match_bases > 0) {
        const double nm_rate = static_cast<double>(nm) /
            static_cast<double>(cigar.total_match_bases);
        if (nm_rate > config_.max_nm_rate) {
            return false;
        }
    }

    return true;
}

}  // namespace placer
