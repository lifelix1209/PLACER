#if PLACER_HAS_ABPOA
extern "C" {
#include <abpoa.h>
}
#endif

#include "abpoa_wrapper.h"
#include <algorithm>
#include <array>
#include <stdexcept>
#include <utility>

namespace placer {

namespace {

inline uint8_t base_to_code(char base) {
    switch (base) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4;
    }
}

inline char code_to_base(uint8_t code) {
    static constexpr char kDecode[] = {'A', 'C', 'G', 'T', 'N'};
    return kDecode[code > 4 ? 4 : code];
}

#if !PLACER_HAS_ABPOA
std::string majority_vote_consensus(const std::vector<std::string>& sequences) {
    if (sequences.empty()) return {};

    size_t max_len = 0;
    for (const auto& seq : sequences) {
        max_len = std::max(max_len, seq.size());
    }
    if (max_len == 0) return {};

    std::string consensus;
    consensus.reserve(max_len);

    for (size_t pos = 0; pos < max_len; ++pos) {
        std::array<int, 5> counts{0, 0, 0, 0, 0};
        for (const auto& seq : sequences) {
            if (pos >= seq.size()) continue;
            counts[base_to_code(seq[pos])]++;
        }

        int best_idx = 4;
        int best_count = counts[4];
        for (int idx = 0; idx < 4; ++idx) {
            if (counts[idx] > best_count) {
                best_count = counts[idx];
                best_idx = idx;
            }
        }
        if (best_count > 0) {
            consensus.push_back(code_to_base(static_cast<uint8_t>(best_idx)));
        }
    }

    return consensus;
}
#endif

}  // namespace

// ============================================================================
// AbPOAWrapper Implementation
// ============================================================================

AbPOAWrapper::AbPOAWrapper() : AbPOAWrapper(Config()) {}

AbPOAWrapper::AbPOAWrapper(Config config) : config_(std::move(config)) {
#if PLACER_HAS_ABPOA
    para_ = abpoa_init_para();
    if (!para_) {
        throw std::runtime_error("Failed to initialize abPOA parameters");
    }

    para_->align_mode = config_.align_mode;
    para_->gap_mode = config_.gap_mode;
    para_->match = config_.match;
    para_->mismatch = config_.mismatch;
    para_->gap_open1 = config_.gap_open;
    para_->gap_open2 = config_.gap_open;
    para_->gap_ext1 = config_.gap_extend;
    para_->gap_ext2 = config_.gap_extend;
    para_->wb = config_.extra_bw;
    para_->out_msa = config_.output_msa ? 1 : 0;
    para_->out_cons = config_.output_consensus ? 1 : 0;
    para_->out_gfa = config_.output_gfa ? 1 : 0;
    para_->max_n_cons = config_.max_n_cons;
    para_->min_freq = config_.min_freq;

    abpoa_post_set_para(para_);

    ab_ = abpoa_init();
    if (!ab_) {
        abpoa_free_para(para_);
        para_ = nullptr;
        throw std::runtime_error("Failed to initialize abPOA object");
    }
#endif
}

AbPOAWrapper::~AbPOAWrapper() {
#if PLACER_HAS_ABPOA
    if (ab_) abpoa_free(ab_);
    if (para_) abpoa_free_para(para_);
#endif
}

std::vector<uint8_t> AbPOAWrapper::encode_sequence(const std::string& seq) const {
    std::vector<uint8_t> encoded(seq.size());
    for (size_t i = 0; i < seq.size(); ++i) {
        encoded[i] = base_to_code(seq[i]);
    }
    return encoded;
}

std::string AbPOAWrapper::decode_consensus(const uint8_t* cons, int len) const {
    if (!cons || len <= 0) return {};
    std::string result;
    result.reserve(static_cast<size_t>(len));
    for (int i = 0; i < len; ++i) {
        result.push_back(code_to_base(cons[i]));
    }
    return result;
}

AbPOAWrapper::Result AbPOAWrapper::align(const std::vector<std::string>& sequences) {
    Result result;
    result.sequences = sequences;

    if (sequences.empty()) {
        best_score_ = 0;
        last_consensus_.clear();
        return result;
    }

    reset();

    if (sequences.size() == 1) {
        result.consensus = sequences[0];
        best_score_ = static_cast<int>(result.consensus.size());
        result.best_score = best_score_;
        last_consensus_ = result.consensus;
        return result;
    }

#if PLACER_HAS_ABPOA
    int n_seqs = static_cast<int>(sequences.size());
    std::vector<int> seq_lens(n_seqs);
    std::vector<uint8_t*> seqs(n_seqs);
    std::vector<std::vector<uint8_t>> encoded_data;

    encoded_data.reserve(static_cast<size_t>(n_seqs));
    for (int i = 0; i < n_seqs; ++i) {
        seq_lens[i] = static_cast<int>(sequences[i].size());
        encoded_data.push_back(encode_sequence(sequences[i]));
        seqs[i] = encoded_data[static_cast<size_t>(i)].data();
    }

    best_score_ = abpoa_msa(ab_, para_, n_seqs, nullptr, seq_lens.data(),
                            seqs.data(), nullptr, nullptr);

    if (para_->out_cons && ab_->abc && ab_->abc->n_cons > 0) {
        result.consensus = decode_consensus(ab_->abc->cons_base[0], ab_->abc->cons_len[0]);
    } else {
        result.consensus = sequences[0];
    }
#else
    result.consensus = majority_vote_consensus(sequences);
    if (result.consensus.empty()) {
        result.consensus = sequences[0];
    }
    best_score_ = static_cast<int>(result.consensus.size());
#endif

    result.best_score = best_score_;
    last_consensus_ = result.consensus;
    return result;
}

void AbPOAWrapper::reset() {
#if PLACER_HAS_ABPOA
    if (ab_ && para_) {
        abpoa_reset(ab_, para_, 0);
    }
#endif
    best_score_ = 0;
    last_consensus_.clear();
}

int AbPOAWrapper::add_sequences(const std::vector<std::string>& sequences) {
    if (sequences.empty()) return 0;

#if PLACER_HAS_ABPOA
    int n_seqs = static_cast<int>(sequences.size());
    std::vector<int> seq_lens(n_seqs);
    std::vector<uint8_t*> seqs(n_seqs);
    std::vector<std::vector<uint8_t>> encoded_data;

    encoded_data.reserve(static_cast<size_t>(n_seqs));
    for (int i = 0; i < n_seqs; ++i) {
        seq_lens[i] = static_cast<int>(sequences[i].size());
        encoded_data.push_back(encode_sequence(sequences[i]));
        seqs[i] = encoded_data[static_cast<size_t>(i)].data();
    }

    best_score_ = abpoa_msa(ab_, para_, n_seqs, nullptr, seq_lens.data(),
                            seqs.data(), nullptr, nullptr);
    last_consensus_ = get_consensus();
    return best_score_;
#else
    auto result = align(sequences);
    return result.best_score;
#endif
}

std::string AbPOAWrapper::get_consensus() const {
#if PLACER_HAS_ABPOA
    if (!ab_ || !ab_->abc || ab_->abc->n_cons == 0) {
        return {};
    }
    return decode_consensus(ab_->abc->cons_base[0], ab_->abc->cons_len[0]);
#else
    return last_consensus_;
#endif
}

// ============================================================================
// BatchPOA Implementation
// ============================================================================

BatchPOA::BatchPOA(int) : wrapper_(std::make_unique<AbPOAWrapper>()) {}

BatchPOA::BatchResult BatchPOA::process(const std::vector<std::string>& sequences) {
    BatchResult result;
    auto align_result = wrapper_->align(sequences);
    result.consensus = std::move(align_result.consensus);
    result.score = wrapper_->get_best_score();
    return result;
}

void BatchPOA::clear() {
    wrapper_->reset();
}

}  // namespace placer
