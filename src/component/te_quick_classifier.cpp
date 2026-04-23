#include "pipeline.h"
#include "te_family_alignment.h"

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

constexpr int32_t kInsertSequenceShortlistTopN = 16;
constexpr int32_t kInsertSequenceMinLen = 80;
constexpr double kInsertSequenceResolvedMinIdentity = 0.78;
constexpr double kInsertSequenceResolvedMinQueryCoverage = 0.72;
constexpr double kInsertSequenceFamilyOnlyMinIdentity = 0.68;
constexpr double kInsertSequenceFamilyOnlyMinQueryCoverage = 0.55;
constexpr double kInsertSequenceUnknownMinIdentity = 0.55;
constexpr double kInsertSequenceUnknownMinQueryCoverage = 0.45;
constexpr double kInsertSequenceUnknownCoarsePrefilterMin = 0.35;
constexpr double kInsertSequenceUnknownCoarseChainCoverageMin = 0.55;
constexpr double kInsertSequenceUnknownHighOccupancyMinQueryCoverage = 0.50;
constexpr double kInsertSequenceFamilyWinnerScoreMin = 0.52;
constexpr double kInsertSequenceExactTemplateMinIdentity = 0.99;
constexpr double kInsertSequenceExactTemplateMinQueryCoverage = 0.99;
constexpr double kInsertSequenceExactTemplateRunnerUpMaxScore = 0.97;

struct TeEntry {
    std::string name;
    std::string sequence;
    std::string reverse_complement_sequence;
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

void populate_kmer_positions(
    const std::string& seq,
    int32_t k,
    KmerPositionMap& out) {
    out.clear();
    out.reserve(seq.size());
    for_each_valid_kmer(seq, k, [&](int32_t start, uint64_t key) {
        out[key].push_back(start);
    });
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

std::string upper_ascii(std::string s) {
    for (char& c : s) {
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    }
    return s;
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

struct TeNameParts {
    std::string exact_name = "NA";
    std::string family = "NA";
    std::string family_key = "NA";
    std::string subfamily = "NA";
};

TeNameParts parse_te_name_parts(const std::string& te_name) {
    TeNameParts parts;
    const std::string token = take_header_token(te_name);
    if (token.empty()) {
        return parts;
    }
    parts.exact_name = token;

    std::string family;
    const size_t hash = token.find('#');
    if (hash != std::string::npos) {
        const std::string left = token.substr(0, hash);
        const std::string right = token.substr(hash + 1);
        parts.subfamily = left.empty() ? token : left;
        family = right.empty() ? parts.subfamily : right;

        const size_t slash = family.rfind('/');
        if (slash != std::string::npos && (slash + 1) < family.size()) {
            family = family.substr(slash + 1);
        }
    } else {
        const size_t colon = token.find(':');
        if (colon != std::string::npos) {
            family = token.substr(0, colon);
            if ((colon + 1) < token.size()) {
                parts.subfamily = token.substr(colon + 1);
            } else {
                parts.subfamily = token;
            }
        } else {
            parts.subfamily = token;
            family = token;
        }
    }

    if (parts.subfamily.empty()) {
        parts.subfamily = token;
    }
    parts.exact_name = parts.subfamily;
    if (family.empty()) {
        family = parts.subfamily;
    }

    const std::string family_key = upper_ascii(family);
    if (family_key.rfind("ALU", 0) == 0) {
        parts.exact_name = parts.subfamily;
        parts.family = "ALU";
        parts.family_key = "ALU";
        return parts;
    }
    if (family_key.rfind("L1", 0) == 0) {
        parts.exact_name = parts.subfamily;
        parts.family = "L1";
        parts.family_key = "L1";
        return parts;
    }
    if (family_key.rfind("SVA", 0) == 0) {
        parts.exact_name = parts.subfamily;
        parts.family = "SVA";
        parts.family_key = "SVA";
        return parts;
    }
    if (family_key.rfind("HERV", 0) == 0) {
        parts.exact_name = parts.subfamily;
        parts.family = "HERV";
        parts.family_key = "HERV";
        return parts;
    }

    parts.family = family;
    parts.family_key = family_key.empty() ? "NA" : family_key;
    return parts;
}

struct LocalAlignmentCell {
    int32_t score = 0;
    int32_t matches = 0;
    int32_t alignment_len = 0;
    int32_t query_bases = 0;
};

struct LocalAlignmentSummary {
    double identity = 0.0;
    double query_coverage = 0.0;
    double score = 0.0;
};

bool better_local_alignment_cell(
    const LocalAlignmentCell& lhs,
    const LocalAlignmentCell& rhs,
    int32_t query_len) {
    if (lhs.score != rhs.score) {
        return lhs.score > rhs.score;
    }

    const double lhs_cov = query_len > 0
        ? static_cast<double>(lhs.query_bases) / static_cast<double>(query_len)
        : 0.0;
    const double rhs_cov = query_len > 0
        ? static_cast<double>(rhs.query_bases) / static_cast<double>(query_len)
        : 0.0;
    if (lhs_cov != rhs_cov) {
        return lhs_cov > rhs_cov;
    }

    const double lhs_identity = lhs.alignment_len > 0
        ? static_cast<double>(lhs.matches) / static_cast<double>(lhs.alignment_len)
        : 0.0;
    const double rhs_identity = rhs.alignment_len > 0
        ? static_cast<double>(rhs.matches) / static_cast<double>(rhs.alignment_len)
        : 0.0;
    if (lhs_identity != rhs_identity) {
        return lhs_identity > rhs_identity;
    }

    if (lhs.query_bases != rhs.query_bases) {
        return lhs.query_bases > rhs.query_bases;
    }
    return lhs.alignment_len > rhs.alignment_len;
}

LocalAlignmentSummary local_alignment_summary(
    const std::string& query,
    const std::string& target) {
    LocalAlignmentSummary summary;
    const int32_t n = static_cast<int32_t>(query.size());
    const int32_t m = static_cast<int32_t>(target.size());
    if (n <= 0 || m <= 0) {
        return summary;
    }

    constexpr int32_t kMatchScore = 2;
    constexpr int32_t kMismatchPenalty = 2;
    constexpr int32_t kGapPenalty = 2;

    std::vector<LocalAlignmentCell> prev(static_cast<size_t>(m + 1));
    std::vector<LocalAlignmentCell> curr(static_cast<size_t>(m + 1));
    LocalAlignmentCell best;

    for (int32_t i = 1; i <= n; ++i) {
        curr[0] = LocalAlignmentCell{};
        for (int32_t j = 1; j <= m; ++j) {
            LocalAlignmentCell cell;

            const LocalAlignmentCell& diag_prev = prev[static_cast<size_t>(j - 1)];
            LocalAlignmentCell diag = diag_prev;
            diag.score +=
                (query[static_cast<size_t>(i - 1)] == target[static_cast<size_t>(j - 1)])
                ? kMatchScore
                : -kMismatchPenalty;
            diag.alignment_len += 1;
            diag.query_bases += 1;
            if (query[static_cast<size_t>(i - 1)] == target[static_cast<size_t>(j - 1)]) {
                diag.matches += 1;
            }
            if (diag.score > 0) {
                cell = diag;
            }

            const LocalAlignmentCell& up_prev = prev[static_cast<size_t>(j)];
            LocalAlignmentCell up = up_prev;
            up.score -= kGapPenalty;
            up.alignment_len += 1;
            up.query_bases += 1;
            if (up.score > 0 && better_local_alignment_cell(up, cell, n)) {
                cell = up;
            }

            const LocalAlignmentCell& left_prev = curr[static_cast<size_t>(j - 1)];
            LocalAlignmentCell left = left_prev;
            left.score -= kGapPenalty;
            left.alignment_len += 1;
            if (left.score > 0 && better_local_alignment_cell(left, cell, n)) {
                cell = left;
            }

            curr[static_cast<size_t>(j)] = cell;
            if (better_local_alignment_cell(cell, best, n)) {
                best = cell;
            }
        }
        std::swap(prev, curr);
        std::fill(curr.begin(), curr.end(), LocalAlignmentCell{});
    }

    if (best.alignment_len <= 0 || best.query_bases <= 0) {
        return summary;
    }
    summary.identity = static_cast<double>(best.matches) /
        static_cast<double>(best.alignment_len);
    summary.query_coverage = static_cast<double>(best.query_bases) /
        static_cast<double>(n);
    summary.score = summary.identity * summary.query_coverage;
    return summary;
}

bool better_local_alignment_summary(
    const LocalAlignmentSummary& lhs,
    const LocalAlignmentSummary& rhs) {
    if (lhs.score != rhs.score) {
        return lhs.score > rhs.score;
    }
    if (lhs.identity != rhs.identity) {
        return lhs.identity > rhs.identity;
    }
    return lhs.query_coverage > rhs.query_coverage;
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
        e.reverse_complement_sequence = reverse_complement(e.sequence);
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

void fnv1a_append_byte(uint64_t& hash, uint8_t byte) {
    constexpr uint64_t kFnvPrime = 1099511628211ull;
    hash ^= static_cast<uint64_t>(byte);
    hash *= kFnvPrime;
}

void fnv1a_append_string(uint64_t& hash, const std::string& value) {
    for (unsigned char c : value) {
        fnv1a_append_byte(hash, c);
    }
    fnv1a_append_byte(hash, 0xff);
}

void fnv1a_append_int32(uint64_t& hash, int32_t value) {
    for (int shift = 0; shift < 32; shift += 8) {
        fnv1a_append_byte(hash, static_cast<uint8_t>((value >> shift) & 0xff));
    }
}

std::string build_te_library_cache_key(
    const std::vector<TeEntry>& entries,
    const std::vector<int32_t>& ks,
    int32_t requested_kmer_size,
    int32_t family_representatives) {
    uint64_t hash = 14695981039346656037ull;
    fnv1a_append_int32(hash, static_cast<int32_t>(entries.size()));
    for (const auto& entry : entries) {
        fnv1a_append_string(hash, entry.name);
        fnv1a_append_string(hash, entry.sequence);
    }
    fnv1a_append_int32(hash, requested_kmer_size);
    fnv1a_append_int32(hash, static_cast<int32_t>(ks.size()));
    for (int32_t k : ks) {
        fnv1a_append_int32(hash, k);
    }
    fnv1a_append_int32(hash, family_representatives);

    std::ostringstream out;
    out << std::hex << hash;
    return out.str();
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
            add_sequence(te_id, entries[i].reverse_complement_sequence);
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

struct TEKmerQuickClassifierModule::AlignmentShortlistDb {
    explicit AlignmentShortlistDb(int32_t kmer_size) : k(kmer_size) {}

    int32_t k = 11;
    std::unordered_map<uint64_t, std::vector<int32_t>> kmer_to_te_ids;

    bool build_from_entries(const std::vector<TeEntry>& entries) {
        kmer_to_te_ids.clear();

        for (size_t i = 0; i < entries.size(); ++i) {
            std::unordered_set<uint64_t> te_kmers;
            te_kmers.reserve(entries[i].sequence.size());

            for_each_valid_kmer(entries[i].sequence, k, [&](int32_t /*start*/, uint64_t key) {
                te_kmers.insert(key);
            });

            for_each_valid_kmer(entries[i].reverse_complement_sequence, k, [&](int32_t /*start*/, uint64_t key) {
                te_kmers.insert(key);
            });

            for (uint64_t key : te_kmers) {
                kmer_to_te_ids[key].push_back(static_cast<int32_t>(i));
            }
        }

        return !kmer_to_te_ids.empty();
    }
};

struct TEKmerQuickClassifierModule::TemplateSeedDb {
    explicit TemplateSeedDb(int32_t kmer_size) : k(kmer_size) {}

    int32_t k = 11;
    std::vector<KmerPositionMap> forward_positions;
    std::vector<KmerPositionMap> reverse_positions;

    bool build_from_entries(const std::vector<TeEntry>& entries) {
        forward_positions.clear();
        reverse_positions.clear();
        forward_positions.resize(entries.size());
        reverse_positions.resize(entries.size());

        for (size_t i = 0; i < entries.size(); ++i) {
            if (static_cast<int32_t>(entries[i].sequence.size()) >= k) {
                forward_positions[i].reserve(entries[i].sequence.size());
                populate_kmer_positions(entries[i].sequence, k, forward_positions[i]);
            }
            if (static_cast<int32_t>(entries[i].reverse_complement_sequence.size()) >= k) {
                reverse_positions[i].reserve(entries[i].reverse_complement_sequence.size());
                populate_kmer_positions(
                    entries[i].reverse_complement_sequence,
                    k,
                    reverse_positions[i]);
            }
        }

        return !forward_positions.empty() && !reverse_positions.empty();
    }
};

TEKmerQuickClassifierModule::TEKmerQuickClassifierModule(PipelineConfig config)
    : config_(std::move(config)) {
    struct CacheState {
        std::shared_ptr<const std::vector<std::string>> te_names;
        std::shared_ptr<const std::vector<std::string>> te_sequences;
        std::shared_ptr<const std::vector<std::string>> te_reverse_complement_sequences;
        std::vector<std::shared_ptr<const Index>> indices;
        std::shared_ptr<const Index> primary_index;
        std::shared_ptr<const AlignmentShortlistDb> alignment_shortlist_db;
        std::shared_ptr<const TemplateSeedDb> template_seed_db;
        std::shared_ptr<const TeFamilyAlignmentIndex> family_alignment_index;
        std::shared_ptr<const TeFamilyGroupCache> family_rep_groups;
    };
    static std::mutex cache_mutex;
    static std::unordered_map<std::string, std::shared_ptr<const CacheState>> cache_by_key;

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

    const auto ks = parse_kmer_sizes_csv(config_.te_kmer_sizes_csv, config_.te_kmer_size);
    const std::string cache_key = build_te_library_cache_key(
        entries,
        ks,
        std::max(7, std::min(31, config_.te_kmer_size)),
        std::max(1, config_.te_family_representatives));
    {
        std::lock_guard<std::mutex> lock(cache_mutex);
        const auto it = cache_by_key.find(cache_key);
        if (it != cache_by_key.end()) {
            const auto& cached = it->second;
            te_names_ = cached->te_names;
            te_sequences_ = cached->te_sequences;
            te_reverse_complement_sequences_ = cached->te_reverse_complement_sequences;
            indices_ = cached->indices;
            primary_index_ = cached->primary_index;
            alignment_shortlist_db_ = cached->alignment_shortlist_db;
            template_seed_db_ = cached->template_seed_db;
            family_alignment_index_ = cached->family_alignment_index;
            family_rep_groups_ = cached->family_rep_groups;
            return;
        }
    }

    auto te_names = std::make_shared<std::vector<std::string>>();
    auto te_sequences = std::make_shared<std::vector<std::string>>();
    auto te_reverse_complement_sequences = std::make_shared<std::vector<std::string>>();
    te_names->reserve(entries.size());
    te_sequences->reserve(entries.size());
    te_reverse_complement_sequences->reserve(entries.size());
    for (const auto& e : entries) {
        te_names->push_back(e.name);
        te_sequences->push_back(e.sequence);
        te_reverse_complement_sequences->push_back(e.reverse_complement_sequence);
    }

    std::vector<std::shared_ptr<const Index>> indices;
    for (int32_t k : ks) {
        auto idx = std::make_shared<Index>(k);
        if (!idx->build_from_entries(entries)) {
            continue;
        }
        indices.push_back(std::move(idx));
    }

    if (indices.empty()) {
        std::cerr << "[TEQuick] failed to build any k-mer index from: " << config_.te_fasta_path << "\n";
        return;
    }

    auto shortlist_db = std::make_shared<AlignmentShortlistDb>(
        std::max(7, std::min(31, config_.te_kmer_size)));
    std::shared_ptr<const AlignmentShortlistDb> alignment_shortlist_db;
    if (shortlist_db->build_from_entries(entries)) {
        alignment_shortlist_db = std::move(shortlist_db);
    }
    auto template_seed_db = std::make_shared<TemplateSeedDb>(
        std::max(7, std::min(31, config_.te_kmer_size)));
    std::shared_ptr<const TemplateSeedDb> cached_template_seed_db;
    if (template_seed_db->build_from_entries(entries)) {
        cached_template_seed_db = std::move(template_seed_db);
    }
    const TeFamilyCacheBundle family_cache = build_te_family_cache(
        *te_names,
        *te_sequences,
        *te_reverse_complement_sequences,
        config_.te_family_representatives,
        config_.te_kmer_size);

    auto best_primary = indices.front();
    for (const auto& idx : indices) {
        if (std::abs(idx->k - config_.te_kmer_size) < std::abs(best_primary->k - config_.te_kmer_size)) {
            best_primary = idx;
        }
    }

    auto cached = std::make_shared<CacheState>();
    cached->te_names = std::move(te_names);
    cached->te_sequences = std::move(te_sequences);
    cached->te_reverse_complement_sequences = std::move(te_reverse_complement_sequences);
    cached->indices = indices;
    cached->primary_index = std::move(best_primary);
    cached->alignment_shortlist_db = std::move(alignment_shortlist_db);
    cached->template_seed_db = std::move(cached_template_seed_db);
    cached->family_alignment_index = family_cache.alignment_index;
    cached->family_rep_groups = family_cache.rep_groups;
    {
        std::lock_guard<std::mutex> lock(cache_mutex);
        cache_by_key[cache_key] = cached;
    }

    te_names_ = cached->te_names;
    te_sequences_ = cached->te_sequences;
    te_reverse_complement_sequences_ = cached->te_reverse_complement_sequences;
    indices_ = cached->indices;
    primary_index_ = cached->primary_index;
    alignment_shortlist_db_ = cached->alignment_shortlist_db;
    template_seed_db_ = cached->template_seed_db;
    family_alignment_index_ = cached->family_alignment_index;
    family_rep_groups_ = cached->family_rep_groups;
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
    const auto& te_names = *te_names_;
    const auto& te_sequences = *te_sequences_;

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
        std::clamp(config_.te_low_kmer_support_trigger, 0.0, 1.0);

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
        if (best_id >= 0 && best_id < static_cast<int32_t>(te_names.size())) {
            hit.te_name = parse_te_name_parts(
                te_names[static_cast<size_t>(best_id)]).exact_name;
        } else {
            hit.te_name.clear();
        }

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
                    if (te_id < 0 || te_id >= static_cast<int32_t>(te_sequences.size())) {
                        continue;
                    }
                    const double iden = semiglobal_edit_identity(
                        seq,
                        te_sequences[static_cast<size_t>(te_id)]);
                    if (iden > rescue_best_identity ||
                        (iden == rescue_best_identity && te_id < rescue_best_id)) {
                        rescue_best_identity = iden;
                        rescue_best_id = te_id;
                    }
                }
                if (rescue_best_id >= 0 && rescue_best_identity >= rescue_identity_min) {
                    hit.rescue_used = true;
                    hit.te_name = parse_te_name_parts(
                        te_names[static_cast<size_t>(rescue_best_id)]).exact_name;
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

TEAlignmentEvidence TEKmerQuickClassifierModule::align_insert_sequence(
    const std::string& raw_insert_seq) const {
    TEAlignmentEvidence evidence;
    last_exact_alignments_ = 0;
    if (!family_alignment_index_ || !te_names_ || !te_sequences_ || !te_reverse_complement_sequences_) {
        evidence.qc_reason = "TE_LIBRARY_UNAVAILABLE";
        return evidence;
    }
    const auto& te_names = *te_names_;
    const auto& te_sequences = *te_sequences_;
    const auto& te_reverse_complement_sequences = *te_reverse_complement_sequences_;
    const auto& family_index = *family_alignment_index_;

    const std::string insert_seq = upper_acgt(raw_insert_seq);
    if (insert_seq.empty()) {
        evidence.qc_reason = "EMPTY_INSERT_SEQUENCE";
        return evidence;
    }
    if (static_cast<int32_t>(insert_seq.size()) < kInsertSequenceMinLen) {
        evidence.qc_reason = "INSERT_SEQ_TOO_SHORT";
        return evidence;
    }

    const auto family_scores = score_family_prefilter(
        family_index,
        insert_seq,
        config_.te_kmer_size,
        config_.te_family_topn);
    if (family_scores.empty()) {
        evidence.qc_reason = "NO_TE_ALIGNMENT_SHORTLIST";
        return evidence;
    }

    struct FamilyCandidate {
        int32_t group_index = -1;
        int32_t representative_index = -1;
        bool reverse = false;
        TeNameParts name_parts;
        LocalAlignmentSummary alignment;
        double prefilter_score = 0.0;
        double total_score = 0.0;
    };
    std::vector<FamilyCandidate> family_candidates;
    family_candidates.reserve(family_scores.size());

    for (const auto& family_score : family_scores) {
        if (family_score.group_index < 0 ||
            family_score.group_index >= static_cast<int32_t>(family_index.groups.size())) {
            continue;
        }
        const auto& group = family_index.groups[static_cast<size_t>(family_score.group_index)];
        if (family_score.representative_index < 0 ||
            family_score.representative_index >=
                static_cast<int32_t>(group.representatives.size())) {
            continue;
        }
        const auto& rep = group.representatives[
            static_cast<size_t>(family_score.representative_index)];

        FamilyCandidate candidate;
        candidate.group_index = family_score.group_index;
        candidate.representative_index = family_score.representative_index;
        candidate.reverse = family_score.reverse;
        candidate.prefilter_score = family_score.family_prefilter_score;
        candidate.name_parts.family = group.family;
        candidate.name_parts.family_key = group.family_key;
        candidate.name_parts.exact_name = rep.exact_name;
        candidate.name_parts.subfamily = rep.exact_name;
        candidate.alignment.identity = family_score.chained_query_coverage;
        candidate.alignment.query_coverage = family_score.chained_query_coverage;
        candidate.alignment.score = family_score.best_rep_chain_norm;
        candidate.total_score = family_score.family_prefilter_score;
        family_candidates.push_back(std::move(candidate));
    }

    if (family_candidates.empty()) {
        evidence.qc_reason = "NO_TE_ALIGNMENT_MATCH";
        return evidence;
    }

    std::sort(family_candidates.begin(), family_candidates.end(), [](const FamilyCandidate& lhs, const FamilyCandidate& rhs) {
        if (lhs.total_score != rhs.total_score) {
            return lhs.total_score > rhs.total_score;
        }
        if (lhs.alignment.score != rhs.alignment.score) {
            return lhs.alignment.score > rhs.alignment.score;
        }
        return lhs.name_parts.family < rhs.name_parts.family;
    });

    const int32_t band_kmer_size = std::max(7, std::min(31, config_.te_kmer_size));
    const bool use_cached_template_seeds =
        template_seed_db_ &&
        template_seed_db_->k == band_kmer_size &&
        static_cast<int32_t>(insert_seq.size()) >= band_kmer_size;
    KmerPositionMap insert_kmer_positions;
    if (use_cached_template_seeds) {
        insert_kmer_positions.reserve(insert_seq.size());
        populate_kmer_positions(insert_seq, band_kmer_size, insert_kmer_positions);
    }

    struct TemplateCandidate {
        int32_t te_id = -1;
        TeNameParts name_parts;
        LocalAlignmentSummary alignment;
        bool reverse = false;
        int32_t band_center = 0;
        int32_t band_width = 32;
        bool exact_used = false;
        double approx_score = 0.0;
        double refine_score = 0.0;
    };

    const auto better_template_band = [](const TemplateBandEstimate& lhs, const TemplateBandEstimate& rhs) {
        if (lhs.seed_score != rhs.seed_score) {
            return lhs.seed_score > rhs.seed_score;
        }
        if (lhs.query_coverage != rhs.query_coverage) {
            return lhs.query_coverage > rhs.query_coverage;
        }
        return lhs.band_width < rhs.band_width;
    };

    const auto sort_template_candidates = [](std::vector<TemplateCandidate>& candidates) {
        std::sort(candidates.begin(), candidates.end(), [](const TemplateCandidate& lhs, const TemplateCandidate& rhs) {
            if (lhs.refine_score != rhs.refine_score) {
                return lhs.refine_score > rhs.refine_score;
            }
            if (lhs.alignment.identity != rhs.alignment.identity) {
                return lhs.alignment.identity > rhs.alignment.identity;
            }
            return lhs.name_parts.exact_name < rhs.name_parts.exact_name;
        });
    };

    struct RefinedFamilyCandidate {
        int32_t group_index = -1;
        double shortlist_score = 0.0;
        std::vector<TemplateCandidate> template_candidates;
    };

    const int32_t refine_family_topn = std::max(1, config_.te_template_refine_topn);
    std::vector<RefinedFamilyCandidate> refined_families;
    refined_families.reserve(std::min(
        static_cast<int32_t>(family_candidates.size()),
        refine_family_topn));
    for (int32_t family_rank = 0;
         family_rank < static_cast<int32_t>(family_candidates.size()) &&
         family_rank < refine_family_topn;
         ++family_rank) {
        const auto& family_candidate = family_candidates[static_cast<size_t>(family_rank)];
        const auto& group = family_index.groups[
            static_cast<size_t>(family_candidate.group_index)];

        RefinedFamilyCandidate refined_family;
        refined_family.group_index = family_candidate.group_index;
        refined_family.shortlist_score = family_candidate.total_score;
        refined_family.template_candidates.reserve(group.member_te_ids.size());
        for (int32_t te_id : group.member_te_ids) {
            if (te_id < 0 || te_id >= static_cast<int32_t>(te_sequences.size()) ||
                te_id >= static_cast<int32_t>(te_reverse_complement_sequences.size()) ||
                te_id >= static_cast<int32_t>(te_names.size())) {
                continue;
            }
            TemplateBandEstimate forward;
            TemplateBandEstimate reverse;
            if (use_cached_template_seeds &&
                te_id < static_cast<int32_t>(template_seed_db_->forward_positions.size()) &&
                te_id < static_cast<int32_t>(template_seed_db_->reverse_positions.size())) {
                forward = estimate_template_band_from_positions(
                    insert_kmer_positions,
                    static_cast<int32_t>(insert_seq.size()),
                    template_seed_db_->forward_positions[static_cast<size_t>(te_id)],
                    static_cast<int32_t>(te_sequences[static_cast<size_t>(te_id)].size()),
                    band_kmer_size);
                reverse = estimate_template_band_from_positions(
                    insert_kmer_positions,
                    static_cast<int32_t>(insert_seq.size()),
                    template_seed_db_->reverse_positions[static_cast<size_t>(te_id)],
                    static_cast<int32_t>(
                        te_reverse_complement_sequences[static_cast<size_t>(te_id)].size()),
                    band_kmer_size);
            } else {
                forward = estimate_template_band(
                    insert_seq,
                    te_sequences[static_cast<size_t>(te_id)],
                    band_kmer_size);
                reverse = estimate_template_band(
                    insert_seq,
                    te_reverse_complement_sequences[static_cast<size_t>(te_id)],
                    band_kmer_size);
            }
            const bool use_reverse = better_template_band(reverse, forward);
            const TemplateBandEstimate band = use_reverse ? reverse : forward;
            if (band.seed_score <= 0.0 && band.query_coverage <= 0.0) {
                continue;
            }

            TemplateCandidate candidate;
            candidate.te_id = te_id;
            candidate.name_parts = parse_te_name_parts(te_names[static_cast<size_t>(te_id)]);
            candidate.alignment.identity = band.query_coverage;
            candidate.alignment.query_coverage = band.query_coverage;
            candidate.alignment.score = band.seed_score;
            candidate.reverse = use_reverse;
            candidate.band_center = band.band_center;
            candidate.band_width = band.band_width;
            candidate.approx_score = band.seed_score;
            candidate.refine_score = band.seed_score;
            refined_family.template_candidates.push_back(std::move(candidate));
        }
        if (refined_family.template_candidates.empty()) {
            continue;
        }
        sort_template_candidates(refined_family.template_candidates);
        refined_families.push_back(std::move(refined_family));
    }

    if (refined_families.empty()) {
        evidence.qc_reason = "NO_TE_ALIGNMENT_MATCH";
        return evidence;
    }

    struct ExactWorkItem {
        size_t family_index = 0;
        size_t template_index = 0;
        double approx_score = 0.0;
        double shortlist_score = 0.0;
    };

    std::vector<ExactWorkItem> exact_work;
    for (size_t family_index = 0; family_index < refined_families.size(); ++family_index) {
        const auto& refined_family = refined_families[family_index];
        for (size_t template_index = 0;
             template_index < refined_family.template_candidates.size();
             ++template_index) {
            const auto& candidate = refined_family.template_candidates[template_index];
            exact_work.push_back({
                family_index,
                template_index,
                candidate.approx_score,
                refined_family.shortlist_score,
            });
        }
    }

    std::sort(exact_work.begin(), exact_work.end(), [](const ExactWorkItem& lhs, const ExactWorkItem& rhs) {
        if (lhs.approx_score != rhs.approx_score) {
            return lhs.approx_score > rhs.approx_score;
        }
        return lhs.shortlist_score > rhs.shortlist_score;
    });

    const int32_t exact_limit = std::max(1, config_.te_exact_align_topn);
    for (int32_t i = 0;
         i < static_cast<int32_t>(exact_work.size()) && i < exact_limit;
         ++i) {
        const ExactWorkItem& work = exact_work[static_cast<size_t>(i)];
        auto& candidate =
            refined_families[work.family_index].template_candidates[work.template_index];
        const std::string& target_seq = candidate.reverse
            ? te_reverse_complement_sequences[static_cast<size_t>(candidate.te_id)]
            : te_sequences[static_cast<size_t>(candidate.te_id)];
        const ExactAlignmentSummary exact = banded_semiglobal_affine_align(
            insert_seq,
            target_seq,
            candidate.band_center,
            candidate.band_width);
        candidate.exact_used = true;
        ++last_exact_alignments_;
        if (exact.score > 0.0) {
            candidate.alignment.identity = exact.identity;
            candidate.alignment.query_coverage = exact.query_coverage;
            candidate.alignment.score = exact.score;
            candidate.refine_score = exact.score;
        }
    }

    for (auto& refined_family : refined_families) {
        sort_template_candidates(refined_family.template_candidates);
    }

    std::sort(refined_families.begin(), refined_families.end(), [](const RefinedFamilyCandidate& lhs, const RefinedFamilyCandidate& rhs) {
        const TemplateCandidate& lhs_best = lhs.template_candidates.front();
        const TemplateCandidate& rhs_best = rhs.template_candidates.front();
        if (lhs_best.refine_score != rhs_best.refine_score) {
            return lhs_best.refine_score > rhs_best.refine_score;
        }
        if (lhs.shortlist_score != rhs.shortlist_score) {
            return lhs.shortlist_score > rhs.shortlist_score;
        }
        return lhs_best.name_parts.family < rhs_best.name_parts.family;
    });

    std::vector<TemplateCandidate>& template_candidates =
        refined_families.front().template_candidates;
    const TemplateCandidate& best_template = template_candidates.front();
    evidence.best_family = best_template.name_parts.family;
    evidence.second_family = "NA";
    evidence.second_score = 0.0;
    if (refined_families.size() >= 2) {
        evidence.second_family =
            refined_families[1].template_candidates.front().name_parts.family;
        evidence.second_score =
            refined_families[1].template_candidates.front().refine_score;
    } else if (family_candidates.size() >= 2) {
        evidence.second_family = family_candidates[1].name_parts.family;
        evidence.second_score = family_candidates[1].alignment.score;
    }

    evidence.best_subfamily = best_template.name_parts.subfamily;
    evidence.best_identity = best_template.alignment.identity;
    evidence.best_query_coverage = best_template.alignment.query_coverage;
    evidence.best_score = best_template.refine_score;
    const double coarse_prefilter_score = family_candidates.front().prefilter_score;
    const double coarse_chain_coverage = family_candidates.front().alignment.query_coverage;
    evidence.coarse_prefilter_score = coarse_prefilter_score;
    evidence.coarse_chain_coverage = coarse_chain_coverage;
    evidence.cross_family_margin = std::max(
        0.0,
        evidence.best_score - evidence.second_score);

    const bool supports_unknown_call =
        evidence.best_identity + 1e-9 >= kInsertSequenceUnknownMinIdentity &&
        evidence.best_query_coverage + 1e-9 >= kInsertSequenceUnknownMinQueryCoverage;
    const bool supports_coarse_unknown_call =
        coarse_prefilter_score + 1e-9 >= kInsertSequenceUnknownCoarsePrefilterMin &&
        coarse_chain_coverage + 1e-9 >= kInsertSequenceUnknownCoarseChainCoverageMin &&
        evidence.best_query_coverage + 1e-9 >= kInsertSequenceUnknownMinQueryCoverage;
    const bool supports_high_occupancy_unknown_call =
        evidence.best_identity + 1e-9 >= 0.50 &&
        evidence.best_query_coverage + 1e-9 >= kInsertSequenceUnknownHighOccupancyMinQueryCoverage &&
        evidence.cross_family_margin >
            (std::max(config_.te_family_margin_min, 0.05) + 0.02);
    if (!supports_unknown_call) {
        if (supports_coarse_unknown_call || supports_high_occupancy_unknown_call) {
            evidence.best_family = "UNKNOWN";
            evidence.best_subfamily = "UNKNOWN";
            evidence.pass = true;
            evidence.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";
            return evidence;
        }
        if (evidence.best_identity + 1e-9 < kInsertSequenceUnknownMinIdentity) {
            evidence.qc_reason = "TE_ALIGNMENT_LOW_IDENTITY";
        } else {
            evidence.qc_reason = "TE_ALIGNMENT_LOW_QUERY_COVERAGE";
        }
        return evidence;
    }

    const bool family_has_winner =
        evidence.best_score + 1e-9 >= kInsertSequenceFamilyWinnerScoreMin &&
        evidence.cross_family_margin > (config_.te_family_margin_min + 1e-4);
    if (!family_has_winner) {
        evidence.best_family = "UNKNOWN";
        evidence.best_subfamily = "UNKNOWN";
        evidence.pass = true;
        evidence.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";
        return evidence;
    }

    const bool supports_family_only =
        evidence.best_identity + 1e-9 >= kInsertSequenceFamilyOnlyMinIdentity &&
        evidence.best_query_coverage + 1e-9 >= kInsertSequenceFamilyOnlyMinQueryCoverage;
    if (!supports_family_only) {
        evidence.best_family = "UNKNOWN";
        evidence.best_subfamily = "UNKNOWN";
        evidence.pass = true;
        evidence.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";
        return evidence;
    }

    double subfamily_margin = 1.0;
    double second_subfamily_score = 0.0;
    if (template_candidates.size() >= 2) {
        second_subfamily_score = template_candidates[1].refine_score;
        subfamily_margin = std::max(
            0.0,
            template_candidates[0].refine_score - template_candidates[1].refine_score);
    }
    const bool exact_template_rescue =
        evidence.best_identity + 1e-9 >= kInsertSequenceExactTemplateMinIdentity &&
        evidence.best_query_coverage + 1e-9 >= kInsertSequenceExactTemplateMinQueryCoverage &&
        (template_candidates.size() < 2 ||
         second_subfamily_score + 1e-9 <= kInsertSequenceExactTemplateRunnerUpMaxScore);
    const bool supports_subfamily =
        ((evidence.best_identity + 1e-9 >= kInsertSequenceResolvedMinIdentity &&
          evidence.best_query_coverage + 1e-9 >= kInsertSequenceResolvedMinQueryCoverage &&
          subfamily_margin + 1e-9 >= config_.te_subfamily_margin_min) ||
         exact_template_rescue);
    if (!supports_subfamily) {
        evidence.best_subfamily.clear();
        evidence.pass = true;
        evidence.qc_reason = "PASS_INSERT_TE_ALIGNMENT_FAMILY_ONLY";
        return evidence;
    }

    evidence.pass = true;
    evidence.qc_reason = "PASS_INSERT_TE_ALIGNMENT";
    return evidence;
}

}  // namespace placer
