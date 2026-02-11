// Include abPOA C headers FIRST (before wrapper header)
extern "C" {
#include <abpoa.h>
}

#include "abpoa_wrapper.h"
#include <algorithm>
#include <cstring>
#include <iostream>
#include <cstdlib>

namespace placer {

// ============================================================================
// AbPOAWrapper Implementation
// ============================================================================

AbPOAWrapper::AbPOAWrapper() : AbPOAWrapper(Config()) {}

AbPOAWrapper::AbPOAWrapper(Config config) : config_(std::move(config)) {
    // Initialize parameters
    para_ = abpoa_init_para();
    if (!para_) {
        throw std::runtime_error("Failed to initialize abPOA parameters");
    }

    // Configure alignment mode
    para_->align_mode = config_.align_mode;
    para_->gap_mode = config_.gap_mode;

    // Configure scoring
    para_->match = config_.match;
    para_->mismatch = config_.mismatch;
    para_->gap_open1 = config_.gap_open;
    para_->gap_open2 = config_.gap_open;
    para_->gap_ext1 = config_.gap_extend;
    para_->gap_ext2 = config_.gap_extend;

    // Configure banding
    para_->wb = config_.extra_bw;

    // Configure output
    para_->out_msa = config_.output_msa ? 1 : 0;
    para_->out_cons = config_.output_consensus ? 1 : 0;
    para_->out_gfa = config_.output_gfa ? 1 : 0;

    // Configure consensus
    para_->max_n_cons = config_.max_n_cons;
    para_->min_freq = config_.min_freq;

    // Post-process parameters
    abpoa_post_set_para(para_);

    // Initialize abPOA object
    ab_ = abpoa_init();
    if (!ab_) {
        abpoa_free_para(para_);
        throw std::runtime_error("Failed to initialize abPOA object");
    }
}

AbPOAWrapper::~AbPOAWrapper() {
    if (ab_) abpoa_free(ab_);
    if (para_) abpoa_free_para(para_);
}

std::vector<uint8_t> AbPOAWrapper::encode_sequence(const std::string& seq) const {
    std::vector<uint8_t> encoded(seq.size());
    for (size_t i = 0; i < seq.size(); ++i) {
        switch (seq[i]) {
            case 'A': case 'a': encoded[i] = 0; break;
            case 'C': case 'c': encoded[i] = 1; break;
            case 'G': case 'g': encoded[i] = 2; break;
            case 'T': case 't': encoded[i] = 3; break;
            default:  encoded[i] = 4; break;  // N or other
        }
    }
    return encoded;
}

std::string AbPOAWrapper::decode_consensus(const uint8_t* cons, int len) const {
    if (!cons || len <= 0) return "";
    std::string result;
    result.reserve(len);
    static const char decode[] = {'A', 'C', 'G', 'T', 'N'};
    for (int i = 0; i < len; ++i) {
        uint8_t code = cons[i];
        if (code > 4) code = 4;
        result.push_back(decode[code]);
    }
    return result;
}

AbPOAWrapper::Result AbPOAWrapper::align(const std::vector<std::string>& sequences) {
    Result result;

    if (sequences.empty()) {
        return result;
    }

    // Reset the graph
    reset();

    if (sequences.size() == 1) {
        // Single sequence: consensus is the sequence itself
        result.consensus = sequences[0];
        result.sequences = sequences;
        return result;
    }

    // Prepare sequences for abPOA
    int n_seqs = static_cast<int>(sequences.size());

    // Allocate name array (not used, set to NULL)
    char** seq_names = nullptr;

    // Encode sequences
    std::vector<int> seq_lens(n_seqs);
    std::vector<uint8_t*> seqs(n_seqs);
    std::vector<std::vector<uint8_t>> encoded_data;

    encoded_data.reserve(n_seqs);
    for (int i = 0; i < n_seqs; ++i) {
        seq_lens[i] = static_cast<int>(sequences[i].size());
        encoded_data.push_back(encode_sequence(sequences[i]));
        seqs[i] = encoded_data[i].data();
    }

    // No quality weights - use NULL
    int** qual_weights = nullptr;

    // Run abPOA MSA
    best_score_ = abpoa_msa(ab_, para_, n_seqs, seq_names, seq_lens.data(),
                            seqs.data(), qual_weights, nullptr);

    // Extract consensus
    if (para_->out_cons && ab_->abc && ab_->abc->n_cons > 0) {
        // Get first consensus
        int cons_len = ab_->abc->cons_len[0];
        uint8_t* cons_seq = ab_->abc->cons_base[0];
        result.consensus = decode_consensus(cons_seq, cons_len);
    }

    // Copy sequences
    result.sequences = sequences;

    return result;
}

void AbPOAWrapper::reset() {
    if (ab_ && para_) {
        abpoa_reset(ab_, para_, 0);
    }
    best_score_ = 0;
}

int AbPOAWrapper::add_sequences(const std::vector<std::string>& sequences) {
    if (sequences.empty()) return 0;

    // Encode sequences
    int n_seqs = static_cast<int>(sequences.size());
    std::vector<int> seq_lens(n_seqs);
    std::vector<uint8_t*> seqs(n_seqs);
    std::vector<std::vector<uint8_t>> encoded_data;

    for (int i = 0; i < n_seqs; ++i) {
        seq_lens[i] = static_cast<int>(sequences[i].size());
        encoded_data.push_back(encode_sequence(sequences[i]));
        seqs[i] = encoded_data[i].data();
    }

    // Align sequences progressively
    int** qual_weights = nullptr;
    int score = abpoa_msa(ab_, para_, n_seqs, nullptr, seq_lens.data(),
                          seqs.data(), qual_weights, nullptr);

    return score;
}

std::string AbPOAWrapper::get_consensus() const {
    if (!ab_->abc || ab_->abc->n_cons == 0) {
        return "";
    }
    return decode_consensus(ab_->abc->cons_base[0], ab_->abc->cons_len[0]);
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
