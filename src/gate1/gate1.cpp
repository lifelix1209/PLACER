#include "gate1.h"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <cstring>

namespace placer {

// ============= HashTEIndex 实现 =============

HashTEIndex::HashTEIndex(HashTEIndexConfig c) : config_(c) {}

bool HashTEIndex::build_from_sequences(const std::vector<std::pair<int, std::string>>& te_sequences) {
    kmer_map.clear();
    // NOTE: Do NOT clear family_name_to_id or id_to_family_name here
    // These are set up by build_from_fasta before calling this method

    for (const auto& [family_id, seq] : te_sequences) {
        const char* seq_ptr = seq.c_str();
        int seq_len = (int)seq.size();
        int k = config_.kmer_size;

        if (seq_len < k) continue;

        uint64_t kmer = 0;
        bool has_invalid = false;
        int valid_len = 0;

        // First k-mer
        for (int i = 0; i < k; ++i) {
            uint8_t code = char_to_2bit(seq_ptr[i]);
            if (code > 3) {
                has_invalid = true;
                break;
            }
            kmer = (kmer << 2) | code;
        }

        if (!has_invalid) {
            kmer_map[kmer] = family_id;
            valid_len = k;
        }

        // Sliding window
        for (int i = k; i < seq_len; ++i) {
            uint8_t code = char_to_2bit(seq_ptr[i]);
            if (code > 3) {
                // Reset on N or invalid char
                has_invalid = true;
                valid_len = 0;
                continue;
            }

            if (has_invalid) {
                // Need to rebuild kmer from scratch
                kmer = 0;
                valid_len = 0;
                for (int j = 0; j < k && (i + j) < seq_len; ++j) {
                    uint8_t c2 = char_to_2bit(seq_ptr[i + j]);
                    if (c2 > 3) {
                        valid_len = 0;
                        break;
                    }
                    kmer = (kmer << 2) | c2;
                    valid_len = j + 1;
                }
                has_invalid = (valid_len < k);
            } else {
                // Rolling update
                uint64_t mask = (1ULL << (2 * k)) - 1;
                uint8_t outgoing = char_to_2bit(seq_ptr[i - k]);
                kmer = roll_kmer(kmer, outgoing, code, mask);
                valid_len = k;
            }

            if (valid_len >= k && !has_invalid) {
                kmer_map[kmer] = family_id;
            }
        }
    }

    return !kmer_map.empty();
}

TEHit HashTEIndex::query(std::string_view seq) const {
    TEHit result;

    if (seq.empty() || kmer_map.empty()) return result;

    const char* seq_ptr = seq.data();
    int seq_len = (int)seq.size();
    int k = config_.kmer_size;

    if (seq_len < k) return result;

    uint64_t mask = (1ULL << (2 * k)) - 1;
    size_t total_kmers = 0;

    uint64_t kmer = 0;
    bool has_invalid = false;
    int valid_len = 0;

    // First k-mer
    for (int i = 0; i < k; ++i) {
        uint8_t code = char_to_2bit(seq_ptr[i]);
        if (code > 3) {
            has_invalid = true;
            break;
        }
        kmer = (kmer << 2) | code;
    }

    if (!has_invalid) {
        ++total_kmers;
        auto it = kmer_map.find(kmer);
        if (it != kmer_map.end()) {
            ++result.hit_count;
            ++result.family_hits[it->second];
        }
        valid_len = k;
    }

    // Sliding window - O(L) single pass
    for (int i = k; i < seq_len; ++i) {
        uint8_t code = char_to_2bit(seq_ptr[i]);
        if (code > 3) {
            // Reset on N or invalid char
            has_invalid = true;
            valid_len = 0;
            continue;
        }

        if (has_invalid) {
            // Rebuild kmer from scratch
            kmer = 0;
            valid_len = 0;
            for (int j = 0; j < k && (i + j) < seq_len; ++j) {
                uint8_t c2 = char_to_2bit(seq_ptr[i + j]);
                if (c2 > 3) {
                    valid_len = 0;
                    break;
                }
                kmer = (kmer << 2) | c2;
                valid_len = j + 1;
            }
            has_invalid = (valid_len < k);
        } else {
            // Rolling update
            uint8_t outgoing = char_to_2bit(seq_ptr[i - k]);
            kmer = roll_kmer(kmer, outgoing, code, mask);
            valid_len = k;
        }

        if (valid_len >= k && !has_invalid) {
            ++total_kmers;
            auto it = kmer_map.find(kmer);
            if (it != kmer_map.end()) {
                ++result.hit_count;
                ++result.family_hits[it->second];
            }
        }
    }

    if (total_kmers > 0) {
        result.hit_density = (double)result.hit_count / total_kmers;
    }

    return result;
}

std::unique_ptr<HashTEIndex> HashTEIndex::build_from_fasta(const std::string& fasta_path, const HashTEIndexConfig& config) {
    std::ifstream infile(fasta_path);
    if (!infile.is_open()) return nullptr;

    std::vector<std::pair<int, std::string>> sequences;
    std::vector<std::string> headers;
    std::string line;
    std::string current_seq;

    while (std::getline(infile, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            // Save previous entry
            if (!current_seq.empty()) {
                sequences.emplace_back(-1, current_seq);
            }
            headers.push_back(line.substr(1));
            current_seq.clear();
        } else {
            current_seq += line;
        }
    }

    // Save last entry
    if (!current_seq.empty()) {
        sequences.emplace_back(-1, current_seq);
    }

    infile.close();

    auto index = std::make_unique<HashTEIndex>(config);

    // Register family names and update sequence family IDs
    for (size_t i = 0; i < sequences.size() && i < headers.size(); ++i) {
        int family_id = index->get_or_create_family_id(headers[i]);
        sequences[i].first = family_id;
    }

    if (!index->build_from_sequences(sequences)) {
        return nullptr;
    }

    return index;
}

// ============= Gate1 实现 =============

void Gate1::extract_probes(const ReadSketch& read, std::vector<ProbeFragment>& out_probes) const {
    if (read.sequence.empty()) return;

    std::string_view full_seq = read.sequence;
    int seq_len = (int)full_seq.length();

    // 端部探针
    int end5_len = std::min(seq_len, config_.probe_len);
    out_probes.push_back({full_seq.substr(0, end5_len), 0, 0});

    if (seq_len > config_.probe_len) {
        out_probes.push_back({full_seq.substr(seq_len - config_.probe_len), 0, seq_len - config_.probe_len});
    }

    // CIGAR 驱动的内部探针
    int q_offset = 0;
    for (const auto& [op, len] : read.cigar_ops) {
        bool consume_query = (op == 'M' || op == 'I' || op == 'S' || op == '=' || op == 'X');
        if (consume_query) {
            if (op == 'S' && len >= config_.min_clip_len) {
                if (q_offset + len <= seq_len) {
                    out_probes.push_back({full_seq.substr(q_offset, len), 1, q_offset});
                }
            } else if (op == 'I' && len >= config_.min_ins_len) {
                int start = std::max(0, q_offset - config_.ins_neighborhood / 2);
                int end = std::min(seq_len, q_offset + len + config_.ins_neighborhood / 2);
                out_probes.push_back({full_seq.substr(start, end - start), 2, start});
            }
            q_offset += len;
        }
    }
}

Gate1::Result Gate1::evaluate(const ReadSketch& read) const {
    Result result;
    std::vector<ProbeFragment> probes;
    probes.reserve(16);

    extract_probes(read, probes);
    result.probes = std::move(probes);

    std::unordered_map<int, uint32_t> family_votes;

    for (const auto& probe : result.probes) {
        auto hit = index_->query(probe.sequence);

        if (index_->passes_gate1(hit)) {
            result.passed = true;
        }

        result.total_hits += hit.hit_count;

        for (const auto& [fam_id, count] : hit.family_hits) {
            family_votes[fam_id] += count;
        }
    }

    if (!result.passed && result.total_hits >= (uint32_t)config_.min_total_hits) {
        result.passed = true;
    }

    if (result.passed && !family_votes.empty()) {
        using Pair = decltype(family_votes)::value_type;
        auto max_it = std::max_element(family_votes.begin(), family_votes.end(),
            [](const Pair& a, const Pair& b) { return a.second < b.second; });
        // Convert family ID back to name
        const auto& family_id_to_name = index_->family_ids();
        int fam_id = max_it->first;
        if (fam_id >= 0 && fam_id < (int)family_id_to_name.size()) {
            result.dominant_family = family_id_to_name[fam_id];
        } else {
            result.dominant_family = std::to_string(fam_id);
        }
    }

    return result;
}

bool Gate1::passes(const ReadSketch& read) const {
    std::vector<ProbeFragment> probes;
    probes.reserve(16);
    extract_probes(read, probes);

    for (const auto& probe : probes) {
        auto hit = index_->query(probe.sequence);
        if (index_->passes_gate1(hit)) return true;
    }

    return false;
}

}  // namespace placer
