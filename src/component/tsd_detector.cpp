#include "pipeline.h"

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
        int len = 0;
        char* seq = faidx_fetch_seq(fai, chrom.c_str(), qstart, qend, &len);
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
    if (!config_.tsd_enable || config_.reference_fasta_path.empty()) {
        return;
    }
    auto impl = std::make_shared<Impl>(config_.reference_fasta_path);
    if (impl->valid()) {
        impl_ = std::move(impl);
    }
}

bool TSDDetector::is_enabled() const {
    return static_cast<bool>(impl_) && config_.tsd_enable;
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
