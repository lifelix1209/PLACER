#ifndef PLACER_GATE1_H
#define PLACER_GATE1_H

#include "bam_reader.h"
#include <string_view>
#include <vector>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <cstdint>

namespace placer {

struct TEHit {
    uint32_t hit_count = 0;
    double hit_density = 0.0;
    std::unordered_map<int, uint32_t> family_hits;
};

struct ProbeFragment {
    std::string_view sequence;
    int source_type;
    int read_offset;
};

// Forward declare Config before use in default arguments
struct HashTEIndexConfig {
    int kmer_size = 15;
    int min_hit_count = 3;
    double min_hit_density = 0.05;
};

struct Gate1Config {
    int probe_len = 200;
    int min_clip_len = 100;
    int min_ins_len = 50;
    int ins_neighborhood = 100;
    int min_total_hits = 5;
    double min_density = 0.1;
    // TE-proxy specific thresholds
    int min_hit_count = 3;
    double min_hit_density = 0.05;
};

class HashTEIndex {
public:
    explicit HashTEIndex(HashTEIndexConfig c);

    bool build_from_sequences(const std::vector<std::pair<int, std::string>>& te_sequences);
    TEHit query(std::string_view seq) const;
    bool passes_gate1(const TEHit& hit) const {
        return hit.hit_count >= (uint32_t)config_.min_hit_count && hit.hit_density >= config_.min_hit_density;
    }

    static std::unique_ptr<HashTEIndex> build_from_fasta(const std::string& fasta_path, const HashTEIndexConfig& config);

    size_t size() const { return kmer_map.size(); }

    // Public for implementation access
    const HashTEIndexConfig& config() const { return config_; }
    const std::unordered_map<uint64_t, int>& kmers() const { return kmer_map; }
    const std::unordered_map<std::string, int>& family_names() const { return family_name_to_id; }
    const std::vector<std::string>& family_ids() const { return id_to_family_name; }

    // Register a family name and return its ID
    int get_or_create_family_id(const std::string& family_name) {
        auto it = family_name_to_id.find(family_name);
        if (it != family_name_to_id.end()) {
            return it->second;
        }
        int new_id = (int)id_to_family_name.size();
        family_name_to_id[family_name] = new_id;
        id_to_family_name.push_back(family_name);
        return new_id;
    }

protected:
    HashTEIndexConfig config_;
    // Key: 2-bit encoded k-mer (k <= 31), Value: Family ID
    std::unordered_map<uint64_t, int> kmer_map;
    // Family name -> ID mapping (for robust FASTA parsing)
    std::unordered_map<std::string, int> family_name_to_id;
    std::vector<std::string> id_to_family_name;

    // Static utility methods for bit-packed k-mer
    static inline uint8_t char_to_2bit(char c) {
        switch(c) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': return 3;
            default: return 4; // N or other invalid
        }
    }

    static inline uint64_t make_kmer_key(const char* seq, int kmer_size) {
        uint64_t kmer = 0;
        for (int i = 0; i < kmer_size; ++i) {
            kmer = (kmer << 2) | char_to_2bit(seq[i]);
        }
        return kmer;
    }

    // Rolling update for sliding window (assumes valid chars only)
    static inline uint64_t roll_kmer(uint64_t prev, uint8_t outgoing, uint8_t incoming, uint64_t mask) {
        return ((prev << 2) | incoming) & mask;
    }
};

class Gate1 {
public:
    struct Result {
        bool passed = false;
        uint32_t total_hits = 0;
        std::string dominant_family;
        std::vector<ProbeFragment> probes;
    };

    Gate1(std::shared_ptr<HashTEIndex> index, Gate1Config c) : index_(std::move(index)), config_(c) {}

    Result evaluate(const ReadSketch& read) const;
    bool passes(const ReadSketch& read) const;

    const Gate1Config& get_config() const { return config_; }

    // Probe extraction - public for testing
    void extract_probes(const ReadSketch& read, std::vector<ProbeFragment>& out_probes) const;

private:
    std::shared_ptr<HashTEIndex> index_;
    Gate1Config config_;
};

}  // namespace placer

#endif  // PLACER_GATE1_H
