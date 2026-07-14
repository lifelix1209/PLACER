#include "pipeline.h"
#include "te_family_alignment.h"
#include "te_sequence_explainer.h"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <unistd.h>

namespace placer {
namespace {

constexpr int32_t kInsertSequenceShortlistTopN = 16;

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

bool is_acgt_base(char c) {
    return char_to_2bit(c) <= 3;
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

uint64_t reverse_complement_kmer_key(uint64_t key, int32_t k) {
    uint64_t rc = 0;
    for (int32_t i = 0; i < k; ++i) {
        const uint64_t code = key & uint64_t{3};
        rc = (rc << 2) | (uint64_t{3} - code);
        key >>= 2;
    }
    return rc;
}

uint64_t canonical_kmer_key(uint64_t key, int32_t k) {
    const uint64_t rc = reverse_complement_kmer_key(key, k);
    return std::min(key, rc);
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
    std::string class_label = "NA";
    std::string order_label = "NA";
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

        const size_t slash = right.find('/');
        if (slash != std::string::npos) {
            parts.class_label = slash > 0 ? right.substr(0, slash) : "NA";
            parts.order_label =
                (slash + 1) < right.size() ? right.substr(slash + 1) : "NA";
            const size_t last_slash = family.rfind('/');
            if (last_slash != std::string::npos && (last_slash + 1) < family.size()) {
                family = family.substr(last_slash + 1);
            }
        } else if (!right.empty()) {
            parts.class_label = right;
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
        parts.family = "Alu";
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

void mark_unresolved_te_annotation(TEAlignmentEvidence& evidence) {
    evidence.best_family = "UNKNOWN";
    evidence.best_subfamily = "UNKNOWN";
    evidence.annotation_confidence = "NA";
    evidence.annotation_class = "NA";
    evidence.annotation_order = "NA";
}

std::string confidence_from_qc_reason(const std::string& qc_reason) {
    if (qc_reason == "PASS_INSERT_TE_ALIGNMENT") {
        return "HIGH";
    }
    if (qc_reason == "PASS_INSERT_TE_ALIGNMENT_FAMILY_ONLY") {
        return "MEDIUM";
    }
    if (qc_reason == "PASS_INSERT_TE_ALIGNMENT_UNKNOWN") {
        return "LOW";
    }
    return "NA";
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

std::string shell_single_quote(const std::string& value) {
    std::string out = "'";
    for (char c : value) {
        if (c == '\'') {
            out += "'\\''";
        } else {
            out.push_back(c);
        }
    }
    out += "'";
    return out;
}

int run_shell_command_checked(const std::string& command) {
    const int status = std::system(command.c_str());
    if (status == -1) {
        return -1;
    }
    return status;
}

bool file_exists_nonempty(const std::filesystem::path& path) {
    std::error_code ec;
    return std::filesystem::exists(path, ec) &&
           !ec &&
           std::filesystem::is_regular_file(path, ec) &&
           !ec &&
           std::filesystem::file_size(path, ec) > 0 &&
           !ec;
}

bool blast_db_files_exist(const std::filesystem::path& db_prefix) {
    return file_exists_nonempty(db_prefix.string() + ".nhr") &&
           file_exists_nonempty(db_prefix.string() + ".nin") &&
           file_exists_nonempty(db_prefix.string() + ".nsq");
}

std::filesystem::path blast_work_dir() {
    std::filesystem::path dir = std::filesystem::temp_directory_path() / "placer_te_blast";
    std::error_code ec;
    std::filesystem::create_directories(dir, ec);
    if (ec) {
        throw std::runtime_error(
            "failed to create PLACER BLAST work directory: " + dir.string() +
            ": " + ec.message());
    }
    return dir;
}

std::string ensure_te_blast_db(
    const std::string& te_fasta_path,
    const std::string& makeblastdb_path,
    const std::string& cache_key) {
    if (te_fasta_path.empty()) {
        return {};
    }
    if (makeblastdb_path.empty()) {
        throw std::runtime_error("BLAST+ makeblastdb path is empty");
    }

    const std::filesystem::path db_prefix =
        blast_work_dir() / ("te_library_" + cache_key);
    if (blast_db_files_exist(db_prefix)) {
        return db_prefix.string();
    }

    static std::mutex makeblastdb_mutex;
    std::lock_guard<std::mutex> lock(makeblastdb_mutex);
    if (blast_db_files_exist(db_prefix)) {
        return db_prefix.string();
    }

    const std::string command =
        shell_single_quote(makeblastdb_path) +
        " -in " + shell_single_quote(te_fasta_path) +
        " -dbtype nucl -out " +
        shell_single_quote(db_prefix.string()) +
        " > /dev/null";
    const int status = run_shell_command_checked(command);
    if (status != 0 || !blast_db_files_exist(db_prefix)) {
        throw std::runtime_error(
            "failed to build BLAST database for TE FASTA '" + te_fasta_path +
            "' with makeblastdb '" + makeblastdb_path + "'");
    }
    return db_prefix.string();
}

std::string write_blast_query_fasta(const std::string& insert_seq) {
    static std::atomic<uint64_t> counter{0};
    const uint64_t id = counter.fetch_add(1, std::memory_order_relaxed);
    const std::filesystem::path path =
        blast_work_dir() /
        ("insert_query_" + std::to_string(static_cast<long long>(::getpid())) +
         "_" + std::to_string(id) + ".fa");
    std::ofstream out(path);
    if (!out.is_open()) {
        throw std::runtime_error("failed to write BLAST query FASTA: " + path.string());
    }
    out << ">insert\n";
    constexpr size_t kLineWidth = 80;
    for (size_t offset = 0; offset < insert_seq.size(); offset += kLineWidth) {
        out << insert_seq.substr(offset, kLineWidth) << "\n";
    }
    return path.string();
}

std::string write_blast_batch_query_fasta(
    const std::vector<std::pair<std::string, std::string>>& queries) {
    static std::atomic<uint64_t> counter{0};
    const uint64_t id = counter.fetch_add(1, std::memory_order_relaxed);
    const std::filesystem::path path =
        blast_work_dir() /
        ("insert_batch_query_" + std::to_string(static_cast<long long>(::getpid())) +
         "_" + std::to_string(id) + ".fa");
    std::ofstream out(path);
    if (!out.is_open()) {
        throw std::runtime_error("failed to write BLAST batch query FASTA: " + path.string());
    }
    constexpr size_t kLineWidth = 80;
    for (const auto& query : queries) {
        out << ">" << query.first << "\n";
        for (size_t offset = 0; offset < query.second.size(); offset += kLineWidth) {
            out << query.second.substr(offset, kLineWidth) << "\n";
        }
    }
    return path.string();
}

struct BlastHsp {
    std::string query_id;
    std::string subject_id;
    double identity = 0.0;
    int32_t alignment_length = 0;
    int32_t query_length = 0;
    int32_t query_start = -1;
    int32_t query_end = -1;
    int32_t target_start = -1;
    int32_t target_end = -1;
    double bitscore = 0.0;
    double evalue = 1.0;
};

struct BlastSubjectHit {
    std::string subject_id;
    TeNameParts name_parts;
    double identity = 0.0;
    double query_coverage = 0.0;
    double score = 0.0;
    double bitscore = 0.0;
    double best_evalue = 1.0;
    int32_t query_start = -1;
    int32_t query_end = -1;
    int32_t target_start = -1;
    int32_t target_end = -1;
};

bool parse_blast_hsp_line(const std::string& line, BlastHsp& out) {
    std::string pident;
    std::string length;
    std::string qlen;
    std::string qstart;
    std::string qend;
    std::string sstart;
    std::string send;
    std::string bitscore;
    std::string evalue;
    std::string query_id;
    std::stringstream fields(line);
    if (!(fields >> query_id >> out.subject_id >> pident >> length >> qlen >>
          qstart >> qend >> sstart >> send >> bitscore >> evalue)) {
        return false;
    }
    out.query_id = query_id;
    try {
        out.identity = std::clamp(std::stod(pident) / 100.0, 0.0, 1.0);
        out.alignment_length = std::max(0, std::stoi(length));
        out.query_length = std::max(0, std::stoi(qlen));
        const int32_t qa = std::stoi(qstart);
        const int32_t qb = std::stoi(qend);
        const int32_t ta = std::stoi(sstart);
        const int32_t tb = std::stoi(send);
        out.query_start = std::min(qa, qb) - 1;
        out.query_end = std::max(qa, qb);
        out.target_start = std::min(ta, tb) - 1;
        out.target_end = std::max(ta, tb);
        out.bitscore = std::stod(bitscore);
        out.evalue = std::stod(evalue);
    } catch (const std::exception&) {
        return false;
    }
    return out.query_length > 0 &&
           out.alignment_length > 0 &&
           out.query_end > out.query_start;
}

double covered_fraction_from_intervals(
    int32_t length,
    std::vector<std::pair<int32_t, int32_t>> intervals) {
    if (length <= 0 || intervals.empty()) {
        return 0.0;
    }
    std::sort(intervals.begin(), intervals.end());
    int32_t covered = 0;
    int32_t current_start = -1;
    int32_t current_end = -1;
    for (const auto& interval : intervals) {
        const int32_t start = std::clamp(interval.first, 0, length);
        const int32_t end = std::clamp(interval.second, 0, length);
        if (end <= start) {
            continue;
        }
        if (current_start < 0) {
            current_start = start;
            current_end = end;
            continue;
        }
        if (start <= current_end) {
            current_end = std::max(current_end, end);
            continue;
        }
        covered += current_end - current_start;
        current_start = start;
        current_end = end;
    }
    if (current_start >= 0) {
        covered += current_end - current_start;
    }
    return std::clamp(
        static_cast<double>(covered) / static_cast<double>(length),
        0.0,
        1.0);
}

std::vector<BlastSubjectHit> collapse_blast_hsps(
    const std::vector<BlastHsp>& hsps,
    int32_t query_len) {
    struct Accumulator {
        std::string subject_id;
        std::vector<std::pair<int32_t, int32_t>> query_intervals;
        double identity_weighted = 0.0;
        int32_t aligned_bases = 0;
        double bitscore = 0.0;
        double best_evalue = 1.0;
        int32_t query_start = -1;
        int32_t query_end = -1;
        int32_t target_start = -1;
        int32_t target_end = -1;
    };

    std::unordered_map<std::string, Accumulator> by_subject;
    for (const BlastHsp& hsp : hsps) {
        if (hsp.query_length > 0) {
            query_len = std::max(query_len, hsp.query_length);
        }
        Accumulator& acc = by_subject[hsp.subject_id];
        if (acc.subject_id.empty()) {
            acc.subject_id = hsp.subject_id;
            acc.best_evalue = hsp.evalue;
        }
        acc.query_intervals.push_back({hsp.query_start, hsp.query_end});
        acc.identity_weighted +=
            hsp.identity * static_cast<double>(hsp.alignment_length);
        acc.aligned_bases += hsp.alignment_length;
        acc.bitscore += hsp.bitscore;
        acc.best_evalue = std::min(acc.best_evalue, hsp.evalue);
        if (acc.query_start < 0 || hsp.query_start < acc.query_start) {
            acc.query_start = hsp.query_start;
        }
        acc.query_end = std::max(acc.query_end, hsp.query_end);
        if (acc.target_start < 0 || hsp.target_start < acc.target_start) {
            acc.target_start = hsp.target_start;
        }
        acc.target_end = std::max(acc.target_end, hsp.target_end);
    }

    std::vector<BlastSubjectHit> hits;
    hits.reserve(by_subject.size());
    for (auto& kv : by_subject) {
        Accumulator& acc = kv.second;
        if (acc.aligned_bases <= 0) {
            continue;
        }
        BlastSubjectHit hit;
        hit.subject_id = acc.subject_id;
        hit.name_parts = parse_te_name_parts(acc.subject_id);
        hit.identity =
            std::clamp(acc.identity_weighted / static_cast<double>(acc.aligned_bases), 0.0, 1.0);
        hit.query_coverage = covered_fraction_from_intervals(
            query_len,
            std::move(acc.query_intervals));
        hit.score = hit.identity * hit.query_coverage;
        hit.bitscore = acc.bitscore;
        hit.best_evalue = acc.best_evalue;
        hit.query_start = acc.query_start;
        hit.query_end = acc.query_end;
        hit.target_start = acc.target_start;
        hit.target_end = acc.target_end;
        hits.push_back(std::move(hit));
    }

    std::sort(hits.begin(), hits.end(), [](const BlastSubjectHit& lhs, const BlastSubjectHit& rhs) {
        if (lhs.best_evalue != rhs.best_evalue) {
            return lhs.best_evalue < rhs.best_evalue;
        }
        if (lhs.bitscore != rhs.bitscore) {
            return lhs.bitscore > rhs.bitscore;
        }
        if (lhs.score != rhs.score) {
            return lhs.score > rhs.score;
        }
        return lhs.subject_id < rhs.subject_id;
    });
    return hits;
}

std::vector<BlastSubjectHit> run_blastn_against_te_library(
    const std::string& blastn_path,
    const std::string& blast_db_prefix,
    const std::string& insert_seq) {
    if (blastn_path.empty()) {
        throw std::runtime_error("BLAST+ blastn path is empty");
    }
    if (blast_db_prefix.empty()) {
        throw std::runtime_error("BLAST database prefix is empty");
    }

    const std::string query_path = write_blast_query_fasta(insert_seq);
    const std::filesystem::path output_path =
        std::filesystem::path(query_path).replace_extension(".blast.tsv");
    const std::string command =
        shell_single_quote(blastn_path) +
        " -query " + shell_single_quote(query_path) +
        " -db " + shell_single_quote(blast_db_prefix) +
        " -task blastn" +
        " -dust no -soft_masking false" +
        " -max_target_seqs 25" +
        " -outfmt " +
        shell_single_quote("6 qseqid sseqid pident length qlen qstart qend sstart send bitscore evalue") +
        " -out " + shell_single_quote(output_path.string());
    const int status = run_shell_command_checked(command);
    if (status != 0) {
        std::remove(query_path.c_str());
        std::remove(output_path.string().c_str());
        throw std::runtime_error(
            "blastn failed for insert consensus TE classification with executable '" +
            blastn_path + "'");
    }

    std::ifstream in(output_path);
    if (!in.is_open()) {
        std::remove(query_path.c_str());
        std::remove(output_path.string().c_str());
        throw std::runtime_error(
            "blastn did not create output file: " + output_path.string());
    }

    std::vector<BlastHsp> hsps;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        BlastHsp hsp;
        if (parse_blast_hsp_line(line, hsp)) {
            hsps.push_back(std::move(hsp));
        }
    }
    std::remove(query_path.c_str());
    std::remove(output_path.string().c_str());
    return collapse_blast_hsps(
        hsps,
        static_cast<int32_t>(insert_seq.size()));
}

std::unordered_map<std::string, std::vector<BlastSubjectHit>> run_blastn_batch_against_te_library(
    const std::string& blastn_path,
    const std::string& blast_db_prefix,
    const std::vector<std::pair<std::string, std::string>>& queries) {
    if (blastn_path.empty()) {
        throw std::runtime_error("BLAST+ blastn path is empty");
    }
    if (blast_db_prefix.empty()) {
        throw std::runtime_error("BLAST database prefix is empty");
    }

    std::unordered_map<std::string, std::vector<BlastSubjectHit>> out;
    if (queries.empty()) {
        return out;
    }

    std::unordered_map<std::string, int32_t> query_lengths;
    query_lengths.reserve(queries.size());
    for (const auto& query : queries) {
        query_lengths[query.first] = static_cast<int32_t>(query.second.size());
    }

    const std::string query_path = write_blast_batch_query_fasta(queries);
    const std::filesystem::path output_path =
        std::filesystem::path(query_path).replace_extension(".blast.tsv");
    const std::string command =
        shell_single_quote(blastn_path) +
        " -query " + shell_single_quote(query_path) +
        " -db " + shell_single_quote(blast_db_prefix) +
        " -task blastn" +
        " -dust no -soft_masking false" +
        " -max_target_seqs 25" +
        " -outfmt " +
        shell_single_quote("6 qseqid sseqid pident length qlen qstart qend sstart send bitscore evalue") +
        " -out " + shell_single_quote(output_path.string());
    const int status = run_shell_command_checked(command);
    if (status != 0) {
        std::remove(query_path.c_str());
        std::remove(output_path.string().c_str());
        throw std::runtime_error(
            "blastn failed for batched insert consensus TE classification with executable '" +
            blastn_path + "'");
    }

    std::ifstream in(output_path);
    if (!in.is_open()) {
        std::remove(query_path.c_str());
        std::remove(output_path.string().c_str());
        throw std::runtime_error(
            "blastn did not create output file: " + output_path.string());
    }

    std::unordered_map<std::string, std::vector<BlastHsp>> hsps_by_query;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        BlastHsp hsp;
        if (!parse_blast_hsp_line(line, hsp)) {
            continue;
        }
        if (query_lengths.find(hsp.query_id) == query_lengths.end()) {
            continue;
        }
        hsps_by_query[hsp.query_id].push_back(std::move(hsp));
    }

    for (const auto& query : queries) {
        const auto it = hsps_by_query.find(query.first);
        if (it == hsps_by_query.end()) {
            out.emplace(query.first, std::vector<BlastSubjectHit>{});
            continue;
        }
        out.emplace(
            query.first,
            collapse_blast_hsps(
                it->second,
                query_lengths[query.first]));
    }

    std::remove(query_path.c_str());
    std::remove(output_path.string().c_str());
    return out;
}

TEAlignmentEvidence build_insert_alignment_evidence_from_blast_hits(
    const std::string& insert_seq,
    bool has_blast_db,
    const std::vector<BlastSubjectHit>& hits) {
    TEAlignmentEvidence evidence;
    auto finalize_evidence = [&](double effective_query_coverage) {
        evidence.te_sequence_explanation = explain_te_sequence_structure(
            insert_seq,
            evidence.qc_reason,
            evidence.best_family,
            evidence.best_subfamily,
            evidence.best_identity,
            effective_query_coverage,
            evidence.annotation_residual_fraction,
            evidence.annotation_masked_fraction,
            evidence.cross_family_margin,
            evidence.second_score,
            evidence.sequence_model_label,
            evidence.sequence_model_score);
        return evidence;
    };
    if (!has_blast_db) {
        evidence.qc_reason = "TE_LIBRARY_UNAVAILABLE";
        return finalize_evidence(evidence.best_query_coverage);
    }

    if (insert_seq.empty()) {
        evidence.qc_reason = "EMPTY_INSERT_SEQUENCE";
        return finalize_evidence(evidence.best_query_coverage);
    }

    if (hits.empty()) {
        evidence.qc_reason = "NO_TE_ALIGNMENT_MATCH";
        evidence.sequence_model_label = "TE_MODEL_OUTLIER";
        evidence.sequence_model_score = -0.50;
        return finalize_evidence(evidence.best_query_coverage);
    }

    // Assign the family by aggregating BLAST hits per family on the discriminative
    // score (identity x query-coverage), keeping each family's best copy, then
    // ranking families by that best copy. TE libraries are highly redundant, so
    // ranking by a family's best copy is robust to a single noisy high-scoring
    // copy -- unlike taking the single top hit. e-value only governs whether there
    // is any TE hit at all (already checked above).
    std::vector<std::pair<std::string, const BlastSubjectHit*>> family_best;
    std::unordered_map<std::string, size_t> family_slot;
    for (const BlastSubjectHit& hit : hits) {
        auto slot = family_slot.find(hit.name_parts.family);
        if (slot == family_slot.end()) {
            family_slot.emplace(hit.name_parts.family, family_best.size());
            family_best.emplace_back(hit.name_parts.family, &hit);
        } else if (hit.score > family_best[slot->second].second->score) {
            family_best[slot->second].second = &hit;
        }
    }
    std::sort(family_best.begin(), family_best.end(),
              [](const auto& lhs, const auto& rhs) {
                  if (lhs.second->score != rhs.second->score) {
                      return lhs.second->score > rhs.second->score;
                  }
                  return lhs.second->best_evalue < rhs.second->best_evalue;
              });

    const BlastSubjectHit& best_hit = *family_best.front().second;
    evidence.best_family = family_best.front().first;
    if (family_best.size() > 1) {
        evidence.second_family = family_best[1].first;
        evidence.second_score = family_best[1].second->score;
    } else {
        evidence.second_family = "NA";
        evidence.second_score = 0.0;
    }

    evidence.best_subfamily = best_hit.name_parts.subfamily;
    evidence.best_identity = best_hit.identity;
    evidence.best_query_coverage = best_hit.query_coverage;
    evidence.best_score = best_hit.score;
    evidence.te_consensus_start = best_hit.target_start;
    evidence.te_consensus_end = best_hit.target_end;
    const double effective_query_coverage = evidence.best_query_coverage;
    evidence.annotation_class = best_hit.name_parts.class_label;
    evidence.annotation_order = best_hit.name_parts.order_label;
    evidence.annotation_masked_fraction = 0.0;
    evidence.annotation_residual_fraction =
        std::clamp(1.0 - effective_query_coverage, 0.0, 1.0);
    if (best_hit.query_start >= 0 && best_hit.query_end > best_hit.query_start) {
        std::ostringstream interval;
        interval << "q=" << best_hit.query_start << "-" << best_hit.query_end
                 << ",t=" << best_hit.target_start << "-" << best_hit.target_end
                 << ",id=" << best_hit.identity
                 << ",cov=" << best_hit.query_coverage;
        evidence.annotation_intervals = interval.str();
    }
    evidence.coarse_prefilter_score = best_hit.score;
    evidence.coarse_chain_coverage = evidence.best_query_coverage;
    evidence.cross_family_margin = std::max(
        0.0,
        evidence.best_score - evidence.second_score);
    evidence.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
    evidence.sequence_model_score = evidence.best_score;

    if (evidence.best_family.empty() || evidence.best_family == "NA") {
        evidence.best_family = "UNKNOWN";
        evidence.best_subfamily = "UNKNOWN";
        evidence.pass = true;
        evidence.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";
        evidence.annotation_confidence =
            confidence_from_qc_reason(evidence.qc_reason);
        return finalize_evidence(effective_query_coverage);
    }

    const bool subfamily_ambiguous =
        std::any_of(hits.begin(), hits.end(), [&](const BlastSubjectHit& hit) {
            return &hit != &best_hit &&
                   hit.name_parts.family == evidence.best_family &&
                   hit.name_parts.subfamily != evidence.best_subfamily &&
                   hit.score >= evidence.best_score;
        });
    if (subfamily_ambiguous || evidence.best_subfamily.empty() || evidence.best_subfamily == "NA") {
        evidence.best_subfamily.clear();
        evidence.pass = true;
        evidence.qc_reason = "PASS_INSERT_TE_ALIGNMENT_FAMILY_ONLY";
        evidence.annotation_confidence =
            confidence_from_qc_reason(evidence.qc_reason);
        return finalize_evidence(effective_query_coverage);
    }

    evidence.pass = true;
    evidence.qc_reason = "PASS_INSERT_TE_ALIGNMENT";
    evidence.annotation_confidence = confidence_from_qc_reason(evidence.qc_reason);
    return finalize_evidence(effective_query_coverage);
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

struct TEKmerQuickClassifierModule::AlignmentEvidenceCache {
    mutable std::mutex mutex;
    std::condition_variable cv;
    std::unordered_map<std::string, TEAlignmentEvidence> by_insert_seq;
    std::unordered_set<std::string> in_flight;

    bool lookup(const std::string& insert_seq, TEAlignmentEvidence& out) const {
        std::lock_guard<std::mutex> lock(mutex);
        const auto it = by_insert_seq.find(insert_seq);
        if (it == by_insert_seq.end()) {
            return false;
        }
        out = it->second;
        return true;
    }

    bool reserve_or_wait(const std::string& insert_seq, TEAlignmentEvidence& out) {
        std::unique_lock<std::mutex> lock(mutex);
        for (;;) {
            const auto cached = by_insert_seq.find(insert_seq);
            if (cached != by_insert_seq.end()) {
                out = cached->second;
                return false;
            }
            if (in_flight.insert(insert_seq).second) {
                return true;
            }
            cv.wait(lock);
        }
    }

    void store(const std::string& insert_seq, const TEAlignmentEvidence& evidence) {
        {
            std::lock_guard<std::mutex> lock(mutex);
            by_insert_seq[insert_seq] = evidence;
            in_flight.erase(insert_seq);
        }
        cv.notify_all();
    }

    void cancel(const std::vector<std::string>& insert_seqs) {
        {
            std::lock_guard<std::mutex> lock(mutex);
            for (const std::string& insert_seq : insert_seqs) {
                in_flight.erase(insert_seq);
            }
        }
        cv.notify_all();
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
        std::shared_ptr<const std::string> blast_db_prefix;
        std::shared_ptr<AlignmentEvidenceCache> alignment_evidence_cache;
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
            blast_db_prefix_ = cached->blast_db_prefix;
            alignment_evidence_cache_ = cached->alignment_evidence_cache;
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

    const std::string blast_db_prefix = ensure_te_blast_db(
        config_.te_fasta_path,
        config_.te_makeblastdb_path,
        cache_key);

    std::shared_ptr<const AlignmentShortlistDb> alignment_shortlist_db;
    std::shared_ptr<const TemplateSeedDb> cached_template_seed_db;
    TeFamilyCacheBundle family_cache;
    std::shared_ptr<const Index> best_primary;
    if (!indices.empty()) {
        auto shortlist_db = std::make_shared<AlignmentShortlistDb>(
            std::max(7, std::min(31, config_.te_kmer_size)));
        if (shortlist_db->build_from_entries(entries)) {
            alignment_shortlist_db = std::move(shortlist_db);
        }
        auto template_seed_db = std::make_shared<TemplateSeedDb>(
            std::max(7, std::min(31, config_.te_kmer_size)));
        if (template_seed_db->build_from_entries(entries)) {
            cached_template_seed_db = std::move(template_seed_db);
        }
        family_cache = build_te_family_cache(
            *te_names,
            *te_sequences,
            *te_reverse_complement_sequences,
            config_.te_family_representatives,
            config_.te_kmer_size);

        best_primary = indices.front();
        for (const auto& idx : indices) {
            if (std::abs(idx->k - config_.te_kmer_size) < std::abs(best_primary->k - config_.te_kmer_size)) {
                best_primary = idx;
            }
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
    cached->blast_db_prefix = std::make_shared<const std::string>(blast_db_prefix);
    cached->alignment_evidence_cache = std::make_shared<AlignmentEvidenceCache>();
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
    blast_db_prefix_ = cached->blast_db_prefix;
    alignment_evidence_cache_ = cached->alignment_evidence_cache;
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
    const std::vector<TEAlignmentEvidence> batch =
        align_insert_sequences({raw_insert_seq});
    return batch.empty() ? TEAlignmentEvidence{} : batch.front();
}

std::vector<TEAlignmentEvidence> TEKmerQuickClassifierModule::align_insert_sequences(
    const std::vector<std::string>& raw_insert_seqs) const {
    last_exact_alignments_ = 0;
    TEAlignmentBatchStats stats;
    std::vector<TEAlignmentEvidence> out(raw_insert_seqs.size());
    std::vector<std::string> normalized;
    normalized.reserve(raw_insert_seqs.size());
    std::unordered_map<std::string, std::string> query_id_by_sequence;
    std::vector<std::pair<std::string, std::string>> queries;
    for (size_t i = 0; i < raw_insert_seqs.size(); ++i) {
        const std::string& raw_insert_seq = raw_insert_seqs[i];
        const std::string insert_seq = upper_acgt(raw_insert_seq);
        normalized.push_back(insert_seq);
        if (insert_seq.empty()) {
            continue;
        }
        if (!blast_db_prefix_ || blast_db_prefix_->empty()) {
            continue;
        }
        if (query_id_by_sequence.find(insert_seq) != query_id_by_sequence.end()) {
            stats.cache_misses += 1;
            stats.deduplicated_queries += 1;
            continue;
        }
        if (alignment_evidence_cache_ &&
            !alignment_evidence_cache_->reserve_or_wait(insert_seq, out[i])) {
            stats.cache_hits += 1;
            continue;
        }
        stats.cache_misses += 1;
        const std::string query_id = "q" + std::to_string(query_id_by_sequence.size());
        query_id_by_sequence.emplace(insert_seq, query_id);
        queries.push_back({query_id, insert_seq});
    }

    std::unordered_map<std::string, std::vector<BlastSubjectHit>> hits_by_query;
    if (!queries.empty()) {
        stats.blast_batches = 1;
        stats.blast_queries = static_cast<int64_t>(queries.size());
        try {
            hits_by_query = run_blastn_batch_against_te_library(
                config_.te_blastn_path,
                *blast_db_prefix_,
                queries);
        } catch (...) {
            if (alignment_evidence_cache_) {
                std::vector<std::string> reserved_insert_seqs;
                reserved_insert_seqs.reserve(queries.size());
                for (const auto& query : queries) {
                    reserved_insert_seqs.push_back(query.second);
                }
                alignment_evidence_cache_->cancel(reserved_insert_seqs);
            }
            throw;
        }
    }

    for (size_t i = 0; i < normalized.size(); ++i) {
        const std::string& insert_seq = normalized[i];
        if (!insert_seq.empty() &&
            alignment_evidence_cache_ &&
            alignment_evidence_cache_->lookup(insert_seq, out[i])) {
            continue;
        }
        std::vector<BlastSubjectHit> hits;
        const auto query_it = query_id_by_sequence.find(insert_seq);
        if (query_it != query_id_by_sequence.end()) {
            const auto hit_it = hits_by_query.find(query_it->second);
            if (hit_it != hits_by_query.end()) {
                hits = hit_it->second;
            }
        }
        out[i] = build_insert_alignment_evidence_from_blast_hits(
            insert_seq,
            static_cast<bool>(blast_db_prefix_) && !blast_db_prefix_->empty(),
            hits);
        if (!insert_seq.empty() && alignment_evidence_cache_) {
            alignment_evidence_cache_->store(insert_seq, out[i]);
        }
    }
    last_alignment_batch_stats_ = stats;
    return out;
}

TEAlignmentBatchStats TEKmerQuickClassifierModule::last_alignment_batch_stats() const {
    return last_alignment_batch_stats_;
}

}  // namespace placer
