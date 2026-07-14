#include "pipeline.h"
#include "contig_alias.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <string>
#include <utility>

#include <htslib/faidx.h>

namespace placer {
namespace {

std::string to_upper_acgt(std::string s) {
    for (char& c : s) {
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    }
    return s;
}

bool has_only_acgt(const std::string& s) {
    for (char c : s) {
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
            return false;
        }
    }
    return !s.empty();
}

bool sequence_is_n_rich_reference_context(const std::string& seq) {
    constexpr int32_t kMinContinuousNRun = 20;
    constexpr double kMaxAllowedNFraction = 0.50;
    if (seq.empty()) {
        return false;
    }

    int32_t n_bases = 0;
    int32_t current_run = 0;
    int32_t max_run = 0;
    for (char c : seq) {
        if (c == 'N') {
            ++n_bases;
            ++current_run;
            max_run = std::max(max_run, current_run);
        } else {
            current_run = 0;
        }
    }

    const double n_fraction =
        static_cast<double>(n_bases) / static_cast<double>(seq.size());
    return max_run >= kMinContinuousNRun || n_fraction >= kMaxAllowedNFraction;
}

double background_occurrence_fraction(
    const std::string& region,
    const std::string& motif) {
    if (region.empty() || motif.empty() || region.size() < motif.size()) {
        return 1.0;
    }

    int32_t total = 0;
    int32_t hit = 0;
    for (size_t i = 0; i + motif.size() <= region.size(); ++i) {
        ++total;
        if (region.compare(i, motif.size(), motif) == 0) {
            ++hit;
        }
    }
    if (total <= 0) {
        return 1.0;
    }
    return static_cast<double>(hit) / static_cast<double>(total);
}

}  // namespace

struct TSDDetector::Impl {
    explicit Impl(const std::string& fasta_path) {
        if (fasta_path.empty()) {
            return;
        }
        fai = fai_load(fasta_path.c_str());
        if (fai == nullptr) {
            if (fai_build(fasta_path.c_str()) == 0) {
                fai = fai_load(fasta_path.c_str());
            }
        }
    }

    ~Impl() {
        if (fai != nullptr) {
            fai_destroy(fai);
            fai = nullptr;
        }
    }

    bool valid() const {
        return fai != nullptr;
    }

    std::string resolve_chrom(const std::string& chrom) const {
        if (!valid()) {
            return {};
        }
        for (const std::string& alias : contig_name_aliases(chrom)) {
            if (faidx_has_seq(fai, alias.c_str()) > 0) {
                return alias;
            }
        }
        return {};
    }

    std::string fetch(const std::string& chrom, int32_t start, int32_t end) const {
        if (!valid()) {
            return {};
        }
        if (chrom.empty()) {
            return {};
        }
        if (end <= start) {
            return {};
        }

        const int32_t qstart = std::max(0, start);
        const int32_t qend = std::max(qstart, end - 1);
        const std::string resolved_chrom = resolve_chrom(chrom);
        if (resolved_chrom.empty()) {
            return {};
        }
        int len = 0;
        char* seq = faidx_fetch_seq(fai, resolved_chrom.c_str(), qstart, qend, &len);
        if (seq == nullptr || len <= 0) {
            if (seq != nullptr) {
                std::free(seq);
            }
            return {};
        }

        std::string out(seq, static_cast<size_t>(len));
        std::free(seq);
        return to_upper_acgt(std::move(out));
    }

    faidx_t* fai = nullptr;
};

TSDDetector::TSDDetector(PipelineConfig config)
    : config_(std::move(config)) {
    if (config_.reference_fasta_path.empty()) {
        return;
    }
    // Each Pipeline instance owns its own faidx-backed detector. Parallel
    // workers therefore construct worker-local Pipeline objects instead of
    // sharing a detector across threads.
    auto impl = std::make_shared<Impl>(config_.reference_fasta_path);
    if (impl->valid()) {
        impl_ = std::move(impl);
    }
}

bool TSDDetector::is_enabled() const {
    return can_fetch_reference() && config_.tsd_enable;
}

bool TSDDetector::can_fetch_reference() const {
    return static_cast<bool>(impl_);
}

bool TSDDetector::reference_position_is_poly_n(
    const std::string& chrom,
    int32_t pos) const {
    if (!impl_ || chrom.empty() || pos < 0) {
        return false;
    }
    const std::string base = impl_->fetch(chrom, pos, pos + 1);
    return base.size() == 1 && base[0] == 'N';
}

bool TSDDetector::reference_interval_is_n_rich(
    const std::string& chrom,
    int32_t start,
    int32_t end) const {
    if (!impl_ || chrom.empty()) {
        return false;
    }
    if (end <= start) {
        return false;
    }
    const std::string seq = impl_->fetch(chrom, start, end);
    return sequence_is_n_rich_reference_context(seq);
}

std::string TSDDetector::fetch_window(
    const std::string& chrom,
    int32_t start,
    int32_t end) const {
    if (!impl_) {
        return {};
    }
    return impl_->fetch(chrom, start, end);
}

TsdDetection TSDDetector::detect(
    const std::string& chrom,
    int32_t left_bp,
    int32_t right_bp) const {
    TsdDetection out;
    if (!is_enabled() || chrom.empty()) {
        return out;
    }

    if (left_bp > right_bp) {
        std::swap(left_bp, right_bp);
    }
    left_bp = std::max(0, left_bp);
    right_bp = std::max(0, right_bp);

    const int32_t min_len = std::max(1, config_.tsd_min_len);
    const int32_t max_len = std::max(min_len, config_.tsd_max_len);
    const int32_t flank = std::max(10, config_.tsd_flank_window);
    const double bg_p_max = std::clamp(config_.tsd_bg_p_max, 0.0, 1.0);

    for (int32_t len = max_len; len >= min_len; --len) {
        if (left_bp - len < 0) {
            continue;
        }
        const std::string left = impl_->fetch(chrom, left_bp - len, left_bp);
        const std::string right = impl_->fetch(chrom, right_bp, right_bp + len);
        if (left.size() != static_cast<size_t>(len) || right.size() != static_cast<size_t>(len)) {
            continue;
        }
        if (!has_only_acgt(left) || !has_only_acgt(right)) {
            continue;
        }
        if (left != right) {
            continue;
        }

        const std::string bg_region = impl_->fetch(
            chrom,
            std::max(0, left_bp - flank),
            right_bp + flank);
        const double p = background_occurrence_fraction(bg_region, left);

        out.type = (p <= bg_p_max) ? "DUP" : "UNCERTAIN";
        out.length = len;
        out.sequence = left;
        out.bg_p = p;
        out.significant = (p <= bg_p_max);
        return out;
    }

    const int32_t delta = right_bp - left_bp;
    if (delta >= min_len && delta <= max_len) {
        const std::string del_seq = impl_->fetch(chrom, left_bp, right_bp);
        if (static_cast<int32_t>(del_seq.size()) == delta && has_only_acgt(del_seq)) {
            const std::string bg_region = impl_->fetch(
                chrom,
                std::max(0, left_bp - flank),
                right_bp + flank);
            const double p = background_occurrence_fraction(bg_region, del_seq);
            out.type = (p <= bg_p_max) ? "DEL" : "UNCERTAIN";
            out.length = delta;
            out.sequence = del_seq;
            out.bg_p = p;
            out.significant = (p <= bg_p_max);
            return out;
        }
    }

    return out;
}

}  // namespace placer
