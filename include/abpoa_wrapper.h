#ifndef PLACER_ABPOA_WRAPPER_H
#define PLACER_ABPOA_WRAPPER_H

#include <string>
#include <vector>
#include <memory>
#include <cstdint>

// Include abPOA types (must come before any usage)
#include <abpoa.h>

namespace placer {

/**
 * abPOA Wrapper for PLACER
 *
 * Uses the abPOA library (yangao07/abPOA) for:
 * - Multiple sequence alignment (MSA)
 * - Consensus generation using heaviest bundling
 * - Adaptive banded dynamic programming
 * - SIMD optimization (SSE2/SSE4.1/AVX2)
 */
class AbPOAWrapper {
public:
    struct Config {
        // Alignment mode: 0=global, 1=local, 2=extend
        int align_mode = 1;  // local mode for TE detection
        // Gap mode: 0=linear, 1=affine
        int gap_mode = 1;

        // Scoring (must be positive for abPOA)
        int match = 2;
        int mismatch = 3;    // positive: penalty is mismatch * -1
        int gap_open = 5;    // positive: penalty is gap_open * -1
        int gap_extend = 1;  // positive: penalty is gap_extend * -1

        // Banding
        int extra_bw = 10;  // extra band width

        // Output options
        bool output_msa = false;
        bool output_consensus = true;
        bool output_gfa = false;

        // Consensus
        int max_n_cons = 1;       // max consensus sequences
        double min_freq = 0.0;    // min frequency for consensus
    };

    struct Result {
        std::string consensus;              // consensus sequence
        std::vector<std::string> sequences; // aligned sequences
        int best_score = 0;
    };

    AbPOAWrapper();
    explicit AbPOAWrapper(Config config);
    ~AbPOAWrapper();

    // Non-copyable
    AbPOAWrapper(const AbPOAWrapper&) = delete;
    AbPOAWrapper& operator=(const AbPOAWrapper&) = delete;

    /**
     * Run POA on a set of sequences
     * @param sequences Vector of sequences (ACGT only)
     * @return Result containing consensus and aligned MSA
     */
    Result align(const std::vector<std::string>& sequences);

    /**
     * Reset the POA graph for a new alignment
     */
    void reset();

    /**
     * Add sequences incrementally to existing graph
     */
    int add_sequences(const std::vector<std::string>& sequences);

    /**
     * Get the current consensus from the graph
     */
    std::string get_consensus() const;

    /**
     * Get alignment score
     */
    int get_best_score() const { return best_score_; }

private:
    Config config_;
    abpoa_para_t* para_ = nullptr;
    abpoa_t* ab_ = nullptr;
    int best_score_ = 0;

    // Convert DNA string to abPOA format (uint8_t array with 0-3 encoding)
    std::vector<uint8_t> encode_sequence(const std::string& seq) const;

    // Convert abPOA consensus back to string
    std::string decode_consensus(const uint8_t* cons, int len) const;
};

/**
 * Simple POA wrapper for quick consensus generation
 */
inline std::string generate_consensus(const std::vector<std::string>& sequences) {
    AbPOAWrapper::Config config;
    config.output_consensus = true;
    config.output_msa = false;

    AbPOAWrapper poa(config);
    auto result = poa.align(sequences);
    return result.consensus;
}

/**
 * Batch consensus generator for components
 */

/**
 * Batch consensus generator for components
 */
class BatchPOA {
public:
    struct BatchResult {
        std::string consensus;
        int score;
    };

    explicit BatchPOA(int max_batch_size = 100);

    /**
     * Process a batch of sequences and return consensus
     */
    BatchResult process(const std::vector<std::string>& sequences);

    /**
     * Clear the batch for reuse
     */
    void clear();

private:
    std::unique_ptr<AbPOAWrapper> wrapper_;
};

}  // namespace placer

#endif  // PLACER_ABPOA_WRAPPER_H
