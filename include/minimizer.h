#ifndef PLACER_MINIMIZER_H
#define PLACER_MINIMIZER_H

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <cstring>
#include <numeric>
#include <cmath>
#include <random>
#include <random>

namespace placer {

/**
 * Minimizer-based k-mer set for similarity estimation
 *
 * Uses the "minimizer" scheme: for each k-mer window, take the minimum hash
 * This provides a uniform sampling of k-mers regardless of position
 */
class MinimizerSet {
public:
    struct Config {
        int k = 15;           // k-mer size
        int w = 10;           // window size (number of k-mers in window)
        int min_hash = 64;     // number of minimizers to keep (sorted by hash)
    };

    explicit MinimizerSet(Config config) : config_(config) {}

    /**
     * Extract minimizer set from a sequence
     * @param seq DNA sequence
     * @return unordered_set of minimizer hashes
     */
    std::unordered_set<uint64_t> extract(const std::string& seq) const {
        std::unordered_set<uint64_t> minimizers;
        if (seq.size() < static_cast<size_t>(config_.k)) {
            return minimizers;
        }

        int k = config_.k;
        int w = config_.w;
        int n_kmers = static_cast<int>(seq.size()) - k + 1;

        if (n_kmers <= 0) return minimizers;

        // Sliding window over k-mers
        // For each window of w k-mers, find the minimizer (minimum hash)
        // Keep top min_hash minimizers across entire sequence

        struct HashPos {
            uint64_t hash;
            int position;
        };

        std::vector<HashPos> window_minimizers;
        window_minimizers.reserve(n_kmers / w + 10);

        for (int i = 0; i < n_kmers; ++i) {
            uint64_t hash = hash_kmer(seq.c_str() + i, k);
            // Keep a small sliding window of minimizer candidates
            window_minimizers.push_back({hash, i});

            // When we have w k-mers, find the minimizer and slide
            if (static_cast<int>(window_minimizers.size()) == w) {
                // Find minimizer (minimum hash) in window
                auto min_it = std::min_element(
                    window_minimizers.begin(), window_minimizers.end(),
                    [](const HashPos& a, const HashPos& b) { return a.hash < b.hash; }
                );

                // Add to global minimizers
                minimizers.insert(min_it->hash);

                // Slide: remove k-mers that are now out of window
                int cutoff_pos = min_it->position - w + 1;
                window_minimizers.erase(
                    std::remove_if(window_minimizers.begin(), window_minimizers.end(),
                        [cutoff_pos](const HashPos& hp) { return hp.position < cutoff_pos; }),
                    window_minimizers.end()
                );
            }
        }

        // Add remaining minimizers from last window
        for (const auto& hp : window_minimizers) {
            minimizers.insert(hp.hash);
        }

        return minimizers;
    }

    /**
     * Compute Jaccard similarity between two minimizer sets
     */
    static double jaccard(const std::unordered_set<uint64_t>& a,
                          const std::unordered_set<uint64_t>& b) {
        if (a.empty() && b.empty()) return 1.0;
        if (a.empty() || b.empty()) return 0.0;

        // Intersection
        size_t inter = 0;
        for (const auto& x : a) {
            if (b.count(x)) ++inter;
        }

        // Union
        size_t union_size = a.size() + b.size() - inter;

        return static_cast<double>(inter) / static_cast<double>(union_size);
    }

    /**
     * Compute average Jaccard of one set against a collection
     */
    static double avg_jaccard(const std::unordered_set<uint64_t>& target,
                             const std::vector<std::unordered_set<uint64_t>>& all) {
        if (all.empty()) return 0.0;
        double sum = 0.0;
        for (const auto& s : all) {
            sum += jaccard(target, s);
        }
        return sum / static_cast<double>(all.size());
    }

private:
    Config config_;

    // 2-bit encoding: A=0, C=1, G=2, T=3
    static inline uint8_t char_to_bit(char c) {
        switch (c) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': return 3;
            default: return 0;
        }
    }

    // Simple rolling hash for k-mer (polynomial hash)
    static uint64_t hash_kmer(const char* seq, int k) {
        uint64_t h = 0;
        for (int i = 0; i < k; ++i) {
            h = (h << 2) | char_to_bit(seq[i]);
        }
        return h;
    }
};

/**
 * Jaccard statistics for component similarity analysis
 */
struct JaccardStats {
    double median = 0.0;
    double p10 = 0.0;      // 10th percentile
    double p25 = 0.0;
    double p75 = 0.0;
    double p90 = 0.0;
    double mean = 0.0;
    double std = 0.0;
    size_t n_pairs = 0;

    // Statistics of individual reads
    double avg_individual_jaccard = 0.0;  // average J(i, others)

    // Bucket statistics if split detected
    bool split_triggered = false;
    double bucket_separation = 0.0;  // J(r*, s*) - measure of separation

    // Check if component should be split (thresholds need calibration)
    bool should_split(double median_threshold = 0.75, double p10_threshold = 0.60) const {
        return median <= median_threshold && p10 <= p10_threshold;
    }
};

/**
 * Component similarity analyzer using minimizer Jaccard
 */
class ComponentSimilarityAnalyzer {
public:
    struct Config {
        int k = 15;
        int w = 10;
        int sample_pairs = 200;     // Number of random pairs to sample
        double median_threshold = 0.75;   // Trigger if median <= this
        double p10_threshold = 0.60;     // Trigger if P10 <= this
        int min_bucket_size = 3;           // Minimum reads per bucket
        double min_bucket_frac = 0.2;      // Minimum fraction per bucket
    };

    struct BucketResult {
        std::vector<size_t> bucket_a;  // Indices of reads in bucket A
        std::vector<size_t> bucket_b;  // Indices of reads in bucket B
        double separation = 0.0;        // J(r*, s*) - measure of separation
        bool valid = false;             // If split is valid (meets constraints)
    };

    explicit ComponentSimilarityAnalyzer(Config config) : config_(config) {}

    /**
     * Analyze a component's reads for similarity structure
     */
    JaccardStats analyze(const std::vector<std::string>& sequences) {
        JaccardStats stats;

        if (sequences.size() < 4) {
            return stats;  // Need at least 4 reads for meaningful analysis
        }

        // Extract minimizer sets
        MinimizerSet minimizer({config_.k, config_.w, 64});
        std::vector<std::unordered_set<uint64_t>> minimizer_sets;
        minimizer_sets.reserve(sequences.size());

        for (const auto& seq : sequences) {
            minimizer_sets.push_back(minimizer.extract(seq));
        }

        // Compute pairwise Jaccard for sampled pairs
        std::vector<double> jaccards;
        size_t n = sequences.size();

        if (config_.sample_pairs > 0) {
            size_t max_pairs = n * (n - 1) / 2;
            size_t sample_n = std::min<size_t>(config_.sample_pairs, max_pairs);

            // Random sampling of pairs
            jaccards.reserve(sample_n);

            std::mt19937 rng(42);  // Fixed seed for reproducibility
            std::uniform_int_distribution<size_t> dist(0, n - 1);

            std::unordered_set<uint64_t> seen_pairs;
            seen_pairs.reserve(sample_n * 2);

            size_t attempts = 0;
            while (jaccards.size() < sample_n && attempts < sample_n * 3) {
                ++attempts;
                size_t i = dist(rng);
                size_t j = dist(rng);
                if (i == j) continue;

                // Create unique pair key (smaller first)
                uint64_t key = (i < j) ? (i << 32) | j : (j << 32) | i;
                if (seen_pairs.count(key)) continue;
                seen_pairs.insert(key);

                double jacc = MinimizerSet::jaccard(minimizer_sets[i], minimizer_sets[j]);
                jaccards.push_back(jacc);
            }
        }

        // Compute statistics
        stats.n_pairs = jaccards.size();
        if (!jaccards.empty()) {
            std::sort(jaccards.begin(), jaccards.end());
            stats.mean = std::accumulate(jaccards.begin(), jaccards.end(), 0.0) / jaccards.size();

            // Percentiles
            auto get_percentile = [&](double p) -> double {
                double idx = p * (jaccards.size() - 1);
                size_t low = static_cast<size_t>(idx);
                size_t high = std::min(low + 1, jaccards.size() - 1);
                return jaccards[low] * (1 - (idx - low)) + jaccards[high] * (idx - low);
            };

            stats.p10 = get_percentile(0.10);
            stats.p25 = get_percentile(0.25);
            stats.median = get_percentile(0.50);
            stats.p75 = get_percentile(0.75);
            stats.p90 = get_percentile(0.90);

            // Std deviation
            double sq_sum = 0.0;
            for (double j : jaccards) {
                sq_sum += (j - stats.mean) * (j - stats.mean);
            }
            stats.std = std::sqrt(sq_sum / jaccards.size());
        }

        // Compute average individual Jaccard (how "central" each read is)
        std::vector<double> avg_jaccards(sequences.size(), 0.0);
        for (size_t i = 0; i < sequences.size(); ++i) {
            double sum = 0.0;
            for (size_t j = 0; j < sequences.size(); ++j) {
                if (i != j) {
                    sum += MinimizerSet::jaccard(minimizer_sets[i], minimizer_sets[j]);
                }
            }
            avg_jaccards[i] = sum / (sequences.size() - 1);
        }
        stats.avg_individual_jaccard =
            std::accumulate(avg_jaccards.begin(), avg_jaccards.end(), 0.0) / avg_jaccards.size();

        // Check split trigger
        stats.split_triggered = stats.should_split(config_.median_threshold, config_.p10_threshold);

        return stats;
    }

    /**
     * Split reads into two buckets if similarity structure suggests it
     */
    BucketResult split(const std::vector<std::string>& sequences,
                       const JaccardStats& stats) {
        BucketResult result;

        if (sequences.size() < 6) {
            return result;  // Need enough reads for meaningful split
        }

        if (!stats.split_triggered) {
            return result;
        }

        // Extract minimizer sets
        MinimizerSet minimizer({config_.k, config_.w, 64});
        std::vector<std::unordered_set<uint64_t>> minimizer_sets;
        minimizer_sets.reserve(sequences.size());

        for (const auto& seq : sequences) {
            minimizer_sets.push_back(minimizer.extract(seq));
        }

        size_t n = sequences.size();

        // Find r*: read with maximum average Jaccard (most "central")
        size_t r_star = 0;
        double max_avg_j = -1.0;
        for (size_t i = 0; i < n; ++i) {
            double sum = 0.0;
            for (size_t j = 0; j < n; ++j) {
                if (i != j) {
                    sum += MinimizerSet::jaccard(minimizer_sets[i], minimizer_sets[j]);
                }
            }
            double avg = sum / (n - 1);
            if (avg > max_avg_j) {
                max_avg_j = avg;
                r_star = i;
            }
        }

        // Find s*: read with minimum Jaccard to r* (most "distant")
        size_t s_star = 0;
        double min_j = 2.0;
        for (size_t i = 0; i < n; ++i) {
            if (i == r_star) continue;
            double j = MinimizerSet::jaccard(minimizer_sets[r_star], minimizer_sets[i]);
            if (j < min_j) {
                min_j = j;
                s_star = i;
            }
        }

        result.separation = min_j;

        // Bucket assignment: J(i, r*) >= J(i, s*) → A, else → B
        double threshold = min_j;
        for (size_t i = 0; i < n; ++i) {
            if (i == r_star || i == s_star) continue;  // Handle these specially

            double j_r = MinimizerSet::jaccard(minimizer_sets[i], minimizer_sets[r_star]);
            double j_s = MinimizerSet::jaccard(minimizer_sets[i], minimizer_sets[s_star]);

            if (j_r >= j_s) {
                result.bucket_a.push_back(i);
            } else {
                result.bucket_b.push_back(i);
            }
        }

        // Handle r* and s* based on their Jaccard to each bucket's reads
        // Simplified: put r* in A, s* in B
        result.bucket_a.push_back(r_star);
        result.bucket_b.push_back(s_star);

        // Check constraints
        size_t min_bucket = std::min(result.bucket_a.size(), result.bucket_b.size());
        double min_frac = static_cast<double>(min_bucket) / n;

        result.valid = (min_bucket >= static_cast<size_t>(config_.min_bucket_size)) &&
                       (min_frac >= config_.min_bucket_frac);

        return result;
    }

private:
    Config config_;
};

}  // namespace placer

#endif  // PLACER_MINIMIZER_H
