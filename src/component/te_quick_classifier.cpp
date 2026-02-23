#include "pipeline.h"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace placer {
namespace {

uint8_t char_to_2bit(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return 4;
    }
}

bool build_kmer(const std::string& s, int32_t start, int32_t k, uint64_t& out) {
    out = 0;
    for (int32_t i = 0; i < k; ++i) {
        const uint8_t code = char_to_2bit(s[static_cast<size_t>(start + i)]);
        if (code > 3) {
            return false;
        }
        out = (out << 2) | code;
    }
    return true;
}

std::string take_header_token(const std::string& header) {
    size_t i = 0;
    while (i < header.size() && std::isspace(static_cast<unsigned char>(header[i]))) {
        ++i;
    }
    size_t j = i;
    while (j < header.size() && !std::isspace(static_cast<unsigned char>(header[j]))) {
        ++j;
    }
    if (j > i) {
        return header.substr(i, j - i);
    }
    return header;
}

std::string upper_acgt(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        out.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
    }
    return out;
}

char complement_base(char c) {
    switch (c) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return 'N';
    }
}

std::string reverse_complement(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (auto it = s.rbegin(); it != s.rend(); ++it) {
        out.push_back(complement_base(*it));
    }
    return out;
}

bool is_softclip_source(InsertionFragmentSource source) {
    return source == InsertionFragmentSource::kClipRefLeft ||
           source == InsertionFragmentSource::kClipRefRight;
}

int32_t max_homopolymer_run(const std::string& seq) {
    if (seq.empty()) {
        return 0;
    }

    int32_t best = 1;
    int32_t run = 1;
    for (size_t i = 1; i < seq.size(); ++i) {
        if (seq[i] == seq[i - 1]) {
            ++run;
            best = std::max(best, run);
        } else {
            run = 1;
        }
    }
    return best;
}

double at_fraction(const std::string& seq) {
    int32_t at = 0;
    int32_t total = 0;
    for (char c : seq) {
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
            continue;
        }
        ++total;
        if (c == 'A' || c == 'T') {
            ++at;
        }
    }
    if (total <= 0) {
        return 0.0;
    }
    return static_cast<double>(at) / static_cast<double>(total);
}

bool is_low_complexity_softclip(
    const InsertionFragment& fragment,
    const std::string& seq,
    double at_fraction_min,
    int32_t homopolymer_run_min) {
    if (!is_softclip_source(fragment.source) || seq.empty()) {
        return false;
    }

    return at_fraction(seq) >= at_fraction_min ||
           max_homopolymer_run(seq) >= homopolymer_run_min;
}

std::mutex g_fragment_hits_tsv_mutex;

}  // namespace

struct TEKmerQuickClassifierModule::Index {
    explicit Index(int32_t kmer_size) : k(kmer_size) {}

    int32_t k = 15;
    std::vector<std::string> te_names;

    // Value >=0: TE id
    // Value == -1: ambiguous (k-mer occurs in >=2 TE entries)
    std::unordered_map<uint64_t, int32_t> kmer_to_id;

    bool build_from_fasta(const std::string& fasta_path) {
        std::ifstream in(fasta_path);
        if (!in.is_open()) {
            return false;
        }

        std::string line;
        std::string header;
        std::string seq;

        auto flush = [&]() {
            if (header.empty() || seq.empty()) {
                header.clear();
                seq.clear();
                return;
            }
            const int32_t id = static_cast<int32_t>(te_names.size());
            te_names.push_back(take_header_token(header));
            const std::string normalized = upper_acgt(seq);
            // Index both strands so read-orientation does not drop TE hits.
            add_sequence(id, normalized);
            add_sequence(id, reverse_complement(normalized));
            header.clear();
            seq.clear();
        };

        while (std::getline(in, line)) {
            if (line.empty()) {
                continue;
            }
            if (line[0] == '>') {
                flush();
                header = line.substr(1);
                continue;
            }
            seq += line;
        }
        flush();

        return !te_names.empty() && !kmer_to_id.empty();
    }

    const std::string& te_name(int32_t id) const {
        static const std::string kEmpty;
        if (id < 0 || id >= static_cast<int32_t>(te_names.size())) {
            return kEmpty;
        }
        return te_names[static_cast<size_t>(id)];
    }

    int32_t lookup(uint64_t key) const {
        const auto it = kmer_to_id.find(key);
        if (it == kmer_to_id.end()) {
            return -2;
        }
        return it->second;
    }

    void add_sequence(int32_t te_id, const std::string& seq) {
        if (static_cast<int32_t>(seq.size()) < k) {
            return;
        }

        for (int32_t i = 0; i + k <= static_cast<int32_t>(seq.size()); ++i) {
            uint64_t key;
            if (!build_kmer(seq, i, k, key)) {
                continue;
            }
            auto it = kmer_to_id.find(key);
            if (it == kmer_to_id.end()) {
                kmer_to_id.emplace(key, te_id);
                continue;
            }
            if (it->second != te_id) {
                it->second = -1;
            }
        }
    }
};

TEKmerQuickClassifierModule::TEKmerQuickClassifierModule(PipelineConfig config)
    : config_(std::move(config)) {
    // Prepare TSV output.
    if (!config_.ins_fragment_hits_tsv_path.empty()) {
        std::lock_guard<std::mutex> lock(g_fragment_hits_tsv_mutex);
        std::ofstream out(config_.ins_fragment_hits_tsv_path, std::ios::out | std::ios::trunc);
        if (out.is_open()) {
            out << "fragment_id\tte\tfrag_len\thit_kmers\ttotal_kmers\tcoverage\taligned_len_est\tkmer_support\n";
        }
    }

    if (config_.te_fasta_path.empty()) {
        return;
    }

    auto idx = std::make_shared<Index>(config_.te_kmer_size);
    if (!idx->build_from_fasta(config_.te_fasta_path)) {
        std::cerr << "[TEQuick] failed to build k-mer index from: " << config_.te_fasta_path << "\n";
        return;
    }
    index_ = std::move(idx);
}

bool TEKmerQuickClassifierModule::is_enabled() const {
    return static_cast<bool>(index_);
}

std::vector<FragmentTEHit> TEKmerQuickClassifierModule::classify(
    const std::vector<InsertionFragment>& fragments) const {
    std::vector<FragmentTEHit> hits;
    hits.reserve(fragments.size());

    if (!is_enabled() || !index_ || index_->k <= 0) {
        return hits;
    }

    const bool write_tsv = !config_.ins_fragment_hits_tsv_path.empty();
    const double low_complexity_at_fraction_min =
        std::clamp(config_.te_softclip_low_complexity_at_frac_min, 0.0, 1.0);
    const int32_t low_complexity_homopolymer_min =
        std::max(1, config_.te_softclip_low_complexity_homopolymer_min);
    std::ostringstream tsv_buffer;

    for (const auto& frag : fragments) {
        FragmentTEHit hit;
        hit.fragment_id = frag.fragment_id;
        hit.fragment_len = static_cast<int32_t>(frag.sequence.size());

        const std::string seq = upper_acgt(frag.sequence);
        if (is_low_complexity_softclip(
                frag,
                seq,
                low_complexity_at_fraction_min,
                low_complexity_homopolymer_min)) {
            if (write_tsv) {
                tsv_buffer << hit.fragment_id << "\t" << hit.te_name << "\t" << hit.fragment_len << "\t"
                           << hit.hit_kmers << "\t" << hit.total_kmers << "\t"
                           << hit.coverage << "\t" << hit.aligned_len_est << "\t" << hit.kmer_support << "\n";
            }
            hits.push_back(std::move(hit));
            continue;
        }

        const int32_t k = index_->k;
        if (static_cast<int32_t>(seq.size()) < k) {
            hits.push_back(hit);
            continue;
        }

        std::unordered_map<int32_t, int32_t> counts;
        int32_t total = 0;

        for (int32_t i = 0; i + k <= static_cast<int32_t>(seq.size()); ++i) {
            uint64_t key;
            if (!build_kmer(seq, i, k, key)) {
                continue;
            }
            ++total;
            const int32_t id = index_->lookup(key);
            if (id >= 0) {
                counts[id] += 1;
            }
        }

        hit.total_kmers = total;

        int32_t best_id = -1;
        int32_t best_hits = 0;
        for (const auto& kv : counts) {
            if (kv.second > best_hits) {
                best_hits = kv.second;
                best_id = kv.first;
            }
        }

        hit.hit_kmers = best_hits;
        hit.coverage = (total > 0) ? static_cast<double>(best_hits) / static_cast<double>(total) : 0.0;
        hit.kmer_support = hit.coverage;
        hit.te_name = index_->te_name(best_id);

        // Second pass: longest consecutive run for best_id
        int32_t run = 0;
        int32_t max_run = 0;
        if (best_id >= 0) {
            for (int32_t i = 0; i + k <= static_cast<int32_t>(seq.size()); ++i) {
                uint64_t key;
                if (!build_kmer(seq, i, k, key)) {
                    run = 0;
                    continue;
                }
                const int32_t id = index_->lookup(key);
                if (id == best_id) {
                    ++run;
                    max_run = std::max(max_run, run);
                } else {
                    run = 0;
                }
            }
        }
        hit.aligned_len_est = (max_run > 0) ? (k + max_run - 1) : 0;

        if (write_tsv) {
            tsv_buffer << hit.fragment_id << "\t" << hit.te_name << "\t" << hit.fragment_len << "\t"
                       << hit.hit_kmers << "\t" << hit.total_kmers << "\t"
                       << hit.coverage << "\t" << hit.aligned_len_est << "\t" << hit.kmer_support << "\n";
        }

        hits.push_back(std::move(hit));
    }

    if (write_tsv) {
        const std::string rows = tsv_buffer.str();
        if (!rows.empty()) {
            std::lock_guard<std::mutex> lock(g_fragment_hits_tsv_mutex);
            std::ofstream out(config_.ins_fragment_hits_tsv_path, std::ios::out | std::ios::app);
            if (out.is_open()) {
                out << rows;
            }
        }
    }

    return hits;
}

ClusterTECall TEKmerQuickClassifierModule::vote_cluster(
    const std::vector<FragmentTEHit>& hits) const {
    ClusterTECall call;

    std::vector<FragmentTEHit> filtered;
    filtered.reserve(hits.size());

    for (const auto& h : hits) {
        if (!h.te_name.empty()) {
            filtered.push_back(h);
        }
    }

    call.fragment_count = static_cast<int32_t>(filtered.size());
    if (call.fragment_count < config_.te_min_fragments_for_vote) {
        return call;
    }

    std::unordered_map<std::string, int32_t> votes;
    for (const auto& h : filtered) {
        votes[h.te_name] += 1;
    }

    std::string best_te;
    int32_t best_votes = 0;
    for (const auto& kv : votes) {
        if (kv.second > best_votes) {
            best_te = kv.first;
            best_votes = kv.second;
        }
    }

    call.te_name = best_te;
    call.vote_fraction = static_cast<double>(best_votes) / static_cast<double>(filtered.size());

    std::vector<double> supports;
    for (const auto& h : filtered) {
        if (h.te_name == best_te) {
            supports.push_back(h.kmer_support);
        }
    }
    std::sort(supports.begin(), supports.end());
    if (!supports.empty()) {
        const size_t mid = supports.size() / 2;
        call.median_identity = (supports.size() % 2 == 1)
            ? supports[mid]
            : 0.5 * (supports[mid - 1] + supports[mid]);
    }

    call.passed =
        call.vote_fraction >= config_.te_vote_fraction_min &&
        call.median_identity >= config_.te_median_identity_min;

    return call;
}

}  // namespace placer
