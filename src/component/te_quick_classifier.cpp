#include "pipeline.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <mutex>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace placer {
namespace {

struct TeEntry {
    std::string name;
    std::string sequence;
};

uint8_t char_to_2bit(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return 4;
    }
}

template <typename Fn>
void for_each_valid_kmer(const std::string& seq, int32_t k, Fn&& fn) {
    if (k <= 0 || static_cast<int32_t>(seq.size()) < k) {
        return;
    }

    const uint64_t mask = (k >= 32)
        ? std::numeric_limits<uint64_t>::max()
        : ((uint64_t{1} << (2 * k)) - 1);
    uint64_t key = 0;
    int32_t valid_bases = 0;
    for (int32_t i = 0; i < static_cast<int32_t>(seq.size()); ++i) {
        const uint8_t code = char_to_2bit(seq[static_cast<size_t>(i)]);
        if (code > 3) {
            key = 0;
            valid_bases = 0;
            continue;
        }
        key = ((key << 2) | code) & mask;
        if (valid_bases < k) {
            ++valid_bases;
        }
        if (valid_bases >= k) {
            fn(i - k + 1, key);
        }
    }
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

double shannon_entropy_acgt(const std::string& seq) {
    int32_t counts[4] = {0, 0, 0, 0};
    int32_t total = 0;
    for (char c : seq) {
        switch (c) {
            case 'A':
                counts[0] += 1;
                total += 1;
                break;
            case 'C':
                counts[1] += 1;
                total += 1;
                break;
            case 'G':
                counts[2] += 1;
                total += 1;
                break;
            case 'T':
                counts[3] += 1;
                total += 1;
                break;
            default:
                break;
        }
    }
    if (total <= 0) {
        return 0.0;
    }

    double entropy = 0.0;
    for (int i = 0; i < 4; ++i) {
        if (counts[i] <= 0) {
            continue;
        }
        const double p = static_cast<double>(counts[i]) / static_cast<double>(total);
        entropy -= p * std::log2(p);
    }
    return entropy;
}

double kmer_uniqueness_ratio(const std::string& seq, int32_t k) {
    if (k <= 0 || static_cast<int32_t>(seq.size()) < k) {
        return 0.0;
    }
    int32_t total = 0;
    std::unordered_set<uint64_t> uniq;
    uniq.reserve(seq.size());
    for_each_valid_kmer(seq, k, [&](int32_t /*start*/, uint64_t key) {
        uniq.insert(key);
        total += 1;
    });
    if (total <= 0) {
        return 0.0;
    }
    return static_cast<double>(uniq.size()) / static_cast<double>(total);
}

bool is_low_complexity_softclip(
    const InsertionFragment& fragment,
    const std::string& seq,
    double at_fraction_min,
    int32_t homopolymer_run_min,
    double entropy_min,
    double kmer_uniqueness_min) {
    if (!is_softclip_source(fragment.source) || seq.empty()) {
        return false;
    }

    return at_fraction(seq) >= at_fraction_min ||
           max_homopolymer_run(seq) >= homopolymer_run_min ||
           shannon_entropy_acgt(seq) < std::max(0.0, entropy_min) ||
           kmer_uniqueness_ratio(seq, 5) < std::clamp(kmer_uniqueness_min, 0.0, 1.0);
}

double semiglobal_edit_identity(const std::string& query, const std::string& target) {
    const int32_t n = static_cast<int32_t>(query.size());
    const int32_t m = static_cast<int32_t>(target.size());
    if (n <= 0 || m <= 0) {
        return 0.0;
    }

    std::vector<int32_t> prev(static_cast<size_t>(m + 1), 0);
    std::vector<int32_t> curr(static_cast<size_t>(m + 1), 0);

    for (int32_t i = 1; i <= n; ++i) {
        curr[0] = i;
        for (int32_t j = 1; j <= m; ++j) {
            const int32_t sub_cost =
                prev[static_cast<size_t>(j - 1)] +
                ((query[static_cast<size_t>(i - 1)] == target[static_cast<size_t>(j - 1)]) ? 0 : 1);
            const int32_t del_cost = prev[static_cast<size_t>(j)] + 1;
            const int32_t ins_cost = curr[static_cast<size_t>(j - 1)] + 1;
            curr[static_cast<size_t>(j)] = std::min({sub_cost, del_cost, ins_cost});
        }
        std::swap(prev, curr);
    }

    int32_t best = prev[0];
    for (int32_t j = 1; j <= m; ++j) {
        best = std::min(best, prev[static_cast<size_t>(j)]);
    }
    const double identity =
        1.0 - (static_cast<double>(best) / static_cast<double>(std::max(1, n)));
    return std::clamp(identity, 0.0, 1.0);
}

std::vector<int32_t> parse_kmer_sizes_csv(const std::string& csv, int32_t fallback_k) {
    std::vector<int32_t> out;
    std::stringstream ss(csv);
    std::string token;
    while (std::getline(ss, token, ',')) {
        size_t b = 0;
        while (b < token.size() && std::isspace(static_cast<unsigned char>(token[b]))) {
            ++b;
        }
        size_t e = token.size();
        while (e > b && std::isspace(static_cast<unsigned char>(token[e - 1]))) {
            --e;
        }
        if (e <= b) {
            continue;
        }

        const std::string trimmed = token.substr(b, e - b);
        char* end = nullptr;
        const long parsed = std::strtol(trimmed.c_str(), &end, 10);
        if (end == trimmed.c_str() || (end && *end != '\0')) {
            continue;
        }
        if (parsed < 7 || parsed > 31) {
            continue;
        }
        out.push_back(static_cast<int32_t>(parsed));
    }

    if (fallback_k >= 7 && fallback_k <= 31) {
        out.push_back(fallback_k);
    }
    if (out.empty()) {
        out.push_back(13);
    }

    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

bool load_te_entries_from_fasta(
    const std::string& fasta_path,
    std::vector<TeEntry>& out_entries) {
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
        TeEntry e;
        e.name = take_header_token(header);
        e.sequence = upper_acgt(seq);
        out_entries.push_back(std::move(e));
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
    return !out_entries.empty();
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

    bool build_from_entries(const std::vector<TeEntry>& entries) {
        te_names.clear();
        kmer_to_id.clear();
        te_names.reserve(entries.size());

        for (size_t i = 0; i < entries.size(); ++i) {
            const int32_t te_id = static_cast<int32_t>(i);
            te_names.push_back(entries[i].name);
            add_sequence(te_id, entries[i].sequence);
            add_sequence(te_id, reverse_complement(entries[i].sequence));
        }

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

        for_each_valid_kmer(seq, k, [&](int32_t /*start*/, uint64_t key) {
            auto it = kmer_to_id.find(key);
            if (it == kmer_to_id.end()) {
                kmer_to_id.emplace(key, te_id);
                return;
            }
            if (it->second != te_id) {
                it->second = -1;
            }
        });
    }
};

TEKmerQuickClassifierModule::TEKmerQuickClassifierModule(PipelineConfig config)
    : config_(std::move(config)) {
    if (!config_.ins_fragment_hits_tsv_path.empty()) {
        std::lock_guard<std::mutex> lock(g_fragment_hits_tsv_mutex);
        std::ofstream out(config_.ins_fragment_hits_tsv_path, std::ios::out | std::ios::trunc);
        if (out.is_open()) {
            out << "fragment_id\tte\tfrag_len\thit_kmers\ttotal_kmers\tcoverage\taligned_len_est\tkmer_support\tmultik_support\trescue_used\n";
        }
    }

    if (config_.te_fasta_path.empty()) {
        return;
    }

    std::vector<TeEntry> entries;
    if (!load_te_entries_from_fasta(config_.te_fasta_path, entries)) {
        std::cerr << "[TEQuick] failed to parse TE FASTA: " << config_.te_fasta_path << "\n";
        return;
    }

    te_names_.reserve(entries.size());
    te_sequences_.reserve(entries.size());
    for (const auto& e : entries) {
        te_names_.push_back(e.name);
        te_sequences_.push_back(e.sequence);
    }

    const auto ks = parse_kmer_sizes_csv(config_.te_kmer_sizes_csv, config_.te_kmer_size);
    for (int32_t k : ks) {
        auto idx = std::make_shared<Index>(k);
        if (!idx->build_from_entries(entries)) {
            continue;
        }
        indices_.push_back(std::move(idx));
    }

    if (indices_.empty()) {
        std::cerr << "[TEQuick] failed to build any k-mer index from: " << config_.te_fasta_path << "\n";
        return;
    }

    auto best_primary = indices_.front();
    for (const auto& idx : indices_) {
        if (std::abs(idx->k - config_.te_kmer_size) < std::abs(best_primary->k - config_.te_kmer_size)) {
            best_primary = idx;
        }
    }
    primary_index_ = std::move(best_primary);
}

bool TEKmerQuickClassifierModule::is_enabled() const {
    return static_cast<bool>(primary_index_) && !indices_.empty();
}

std::vector<FragmentTEHit> TEKmerQuickClassifierModule::classify(
    const std::vector<InsertionFragment>& fragments) const {
    std::vector<FragmentTEHit> hits;
    hits.reserve(fragments.size());

    if (!is_enabled()) {
        return hits;
    }

    const bool write_tsv = !config_.ins_fragment_hits_tsv_path.empty();
    const double low_complexity_at_fraction_min =
        std::clamp(config_.te_softclip_low_complexity_at_frac_min, 0.0, 1.0);
    const int32_t low_complexity_homopolymer_min =
        std::max(1, config_.te_softclip_low_complexity_homopolymer_min);
    const double low_complexity_entropy_min =
        std::max(0.0, config_.te_softclip_entropy_min);
    const double low_complexity_kmer_uniqueness_min =
        std::clamp(config_.te_softclip_kmer_uniqueness_min, 0.0, 1.0);
    const int32_t rescue_topn = std::max(1, config_.te_low_kmer_rescue_topn);
    const int32_t rescue_min_frag_len = std::max(1, config_.te_low_kmer_rescue_min_frag_len);
    const double rescue_identity_min =
        std::clamp(config_.te_low_kmer_rescue_identity_min, 0.0, 1.0);
    const double rescue_margin_max =
        std::clamp(config_.te_low_kmer_rescue_margin_max, 0.0, 1.0);
    const double support_gate =
        std::clamp(config_.te_median_identity_min, 0.0, 1.0);

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
                low_complexity_homopolymer_min,
                low_complexity_entropy_min,
                low_complexity_kmer_uniqueness_min)) {
            if (write_tsv) {
                tsv_buffer << hit.fragment_id << "\t" << hit.te_name << "\t" << hit.fragment_len << "\t"
                           << hit.hit_kmers << "\t" << hit.total_kmers << "\t"
                           << hit.coverage << "\t" << hit.aligned_len_est << "\t" << hit.kmer_support << "\t"
                           << hit.multik_support << "\t" << (hit.rescue_used ? 1 : 0) << "\n";
            }
            hits.push_back(std::move(hit));
            continue;
        }

        std::unordered_map<int32_t, double> weighted_support_sum;
        double total_weight = 0.0;

        int32_t primary_total = 0;
        std::unordered_map<int32_t, int32_t> primary_counts;
        int32_t fallback_total = 0;
        std::unordered_map<int32_t, int32_t> fallback_counts;

        for (const auto& idx_ptr : indices_) {
            const Index& idx = *idx_ptr;
            const int32_t k = idx.k;
            if (static_cast<int32_t>(seq.size()) < k) {
                continue;
            }

            std::unordered_map<int32_t, int32_t> counts;
            int32_t total = 0;
            for_each_valid_kmer(seq, k, [&](int32_t /*start*/, uint64_t key) {
                ++total;
                const int32_t id = idx.lookup(key);
                if (id >= 0) {
                    counts[id] += 1;
                }
            });
            if (total <= 0) {
                continue;
            }

            const double w = static_cast<double>(k);
            total_weight += w;
            for (const auto& kv : counts) {
                const double support_k =
                    static_cast<double>(kv.second) / static_cast<double>(total);
                weighted_support_sum[kv.first] += (w * support_k);
            }

            if (idx_ptr.get() == primary_index_.get()) {
                primary_total = total;
                primary_counts = counts;
            }
            if (fallback_total <= 0) {
                fallback_total = total;
                fallback_counts = counts;
            }
        }

        std::vector<std::pair<int32_t, double>> ranked;
        ranked.reserve(weighted_support_sum.size());
        if (total_weight > 0.0) {
            for (const auto& kv : weighted_support_sum) {
                ranked.push_back({kv.first, kv.second / total_weight});
            }
            std::sort(ranked.begin(), ranked.end(), [](const auto& a, const auto& b) {
                if (a.second != b.second) {
                    return a.second > b.second;
                }
                return a.first < b.first;
            });
        }

        int32_t best_id = -1;
        double best_score = 0.0;
        double second_score = 0.0;
        if (!ranked.empty()) {
            best_id = ranked.front().first;
            best_score = ranked.front().second;
            if (ranked.size() >= 2) {
                second_score = ranked[1].second;
            }
        }

        hit.multik_support = best_score;
        hit.coverage = best_score;
        hit.kmer_support = best_score;
        hit.te_name = (best_id >= 0 && best_id < static_cast<int32_t>(te_names_.size()))
            ? te_names_[static_cast<size_t>(best_id)]
            : std::string();

        const auto& count_source = (primary_total > 0) ? primary_counts : fallback_counts;
        hit.total_kmers = (primary_total > 0) ? primary_total : fallback_total;
        if (best_id >= 0) {
            const auto it = count_source.find(best_id);
            if (it != count_source.end()) {
                hit.hit_kmers = it->second;
            }
        }

        const Index* run_index = nullptr;
        for (const auto& idx_ptr : indices_) {
            if (static_cast<int32_t>(seq.size()) < idx_ptr->k) {
                continue;
            }
            if (!run_index || idx_ptr->k > run_index->k) {
                run_index = idx_ptr.get();
            }
        }
        if (run_index && best_id >= 0) {
            int32_t run = 0;
            int32_t max_run = 0;
            int32_t prev_start = -2;
            for_each_valid_kmer(seq, run_index->k, [&](int32_t start, uint64_t key) {
                if (start != (prev_start + 1)) {
                    run = 0;
                }
                const int32_t id = run_index->lookup(key);
                if (id == best_id) {
                    ++run;
                    max_run = std::max(max_run, run);
                } else {
                    run = 0;
                }
                prev_start = start;
            });
            hit.aligned_len_est = (max_run > 0) ? (run_index->k + max_run - 1) : 0;
        }

        if (config_.te_low_kmer_rescue_enable &&
            static_cast<int32_t>(seq.size()) >= rescue_min_frag_len &&
            !ranked.empty()) {
            const bool support_trigger = best_score < support_gate;
            const bool margin_trigger = (ranked.size() >= 2) &&
                ((best_score - second_score) < rescue_margin_max);
            if (support_trigger || margin_trigger) {
                int32_t rescue_best_id = -1;
                double rescue_best_identity = 0.0;
                const int32_t n = std::min(
                    rescue_topn,
                    static_cast<int32_t>(ranked.size()));
                for (int32_t i = 0; i < n; ++i) {
                    const int32_t te_id = ranked[static_cast<size_t>(i)].first;
                    if (te_id < 0 || te_id >= static_cast<int32_t>(te_sequences_.size())) {
                        continue;
                    }
                    const double iden = semiglobal_edit_identity(
                        seq,
                        te_sequences_[static_cast<size_t>(te_id)]);
                    if (iden > rescue_best_identity ||
                        (iden == rescue_best_identity && te_id < rescue_best_id)) {
                        rescue_best_identity = iden;
                        rescue_best_id = te_id;
                    }
                }
                if (rescue_best_id >= 0 && rescue_best_identity >= rescue_identity_min) {
                    hit.rescue_used = true;
                    hit.te_name = te_names_[static_cast<size_t>(rescue_best_id)];
                    hit.kmer_support = std::max(hit.kmer_support, rescue_best_identity);
                    hit.coverage = hit.kmer_support;
                    if (hit.total_kmers > 0) {
                        const auto it = count_source.find(rescue_best_id);
                        if (it != count_source.end()) {
                            hit.hit_kmers = it->second;
                        }
                    }
                }
            }
        }

        if (write_tsv) {
            tsv_buffer << hit.fragment_id << "\t" << hit.te_name << "\t" << hit.fragment_len << "\t"
                       << hit.hit_kmers << "\t" << hit.total_kmers << "\t"
                       << hit.coverage << "\t" << hit.aligned_len_est << "\t" << hit.kmer_support << "\t"
                       << hit.multik_support << "\t" << (hit.rescue_used ? 1 : 0) << "\n";
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

    std::vector<const FragmentTEHit*> filtered;
    filtered.reserve(hits.size());
    for (const auto& h : hits) {
        if (!h.te_name.empty()) {
            filtered.push_back(&h);
        }
    }

    call.fragment_count = static_cast<int32_t>(filtered.size());
    if (call.fragment_count <= 0) {
        return call;
    }

    std::unordered_map<std::string, double> weighted_votes;
    std::unordered_map<std::string, int32_t> raw_votes;
    double total_weight = 0.0;

    for (const auto* h : filtered) {
        const double w = std::max(1e-6, std::clamp(h->kmer_support, 0.0, 1.0));
        weighted_votes[h->te_name] += w;
        raw_votes[h->te_name] += 1;
        total_weight += w;
    }

    std::string best_te;
    std::string second_te;
    double best_w = -1.0;
    double second_w = -1.0;
    int32_t best_raw = -1;
    int32_t second_raw = -1;

    for (const auto& kv : weighted_votes) {
        const auto raw_it = raw_votes.find(kv.first);
        const int32_t raw = (raw_it == raw_votes.end()) ? 0 : raw_it->second;
        const bool better_than_best =
            kv.second > best_w ||
            (kv.second == best_w && raw > best_raw) ||
            (kv.second == best_w && raw == best_raw && (best_te.empty() || kv.first < best_te));
        if (better_than_best) {
            second_te = best_te;
            second_w = best_w;
            second_raw = best_raw;
            best_te = kv.first;
            best_w = kv.second;
            best_raw = raw;
            continue;
        }

        const bool better_than_second =
            kv.second > second_w ||
            (kv.second == second_w && raw > second_raw) ||
            (kv.second == second_w && raw == second_raw && (second_te.empty() || kv.first < second_te));
        if (better_than_second) {
            second_te = kv.first;
            second_w = kv.second;
            second_raw = raw;
        }
    }

    call.top1_te_name = best_te;
    call.top2_te_name = second_te;
    if (total_weight > 0.0) {
        call.posterior_top1 = std::clamp(best_w / total_weight, 0.0, 1.0);
        call.posterior_top2 = std::clamp(std::max(0.0, second_w) / total_weight, 0.0, 1.0);
    }
    call.posterior_margin = std::max(0.0, call.posterior_top1 - call.posterior_top2);

    if (call.fragment_count < config_.te_min_fragments_for_vote) {
        return call;
    }

    call.te_name = best_te;
    call.vote_fraction = call.posterior_top1;

    std::vector<double> supports;
    std::vector<double> multik_supports;
    int32_t rescue_count = 0;
    for (const auto* h : filtered) {
        if (h->te_name == best_te) {
            supports.push_back(h->kmer_support);
            multik_supports.push_back(h->multik_support);
            if (h->rescue_used) {
                rescue_count += 1;
            }
        }
    }
    std::sort(supports.begin(), supports.end());
    std::sort(multik_supports.begin(), multik_supports.end());
    if (!supports.empty()) {
        const size_t mid = supports.size() / 2;
        call.median_identity = (supports.size() % 2 == 1)
            ? supports[mid]
            : 0.5 * (supports[mid - 1] + supports[mid]);
    }
    if (!multik_supports.empty()) {
        const size_t mid = multik_supports.size() / 2;
        call.multik_support = (multik_supports.size() % 2 == 1)
            ? multik_supports[mid]
            : 0.5 * (multik_supports[mid - 1] + multik_supports[mid]);
        call.rescue_frac =
            static_cast<double>(rescue_count) /
            static_cast<double>(multik_supports.size());
    }

    call.passed =
        call.vote_fraction >= config_.te_vote_fraction_min &&
        call.median_identity >= config_.te_median_identity_min;

    return call;
}

}  // namespace placer
