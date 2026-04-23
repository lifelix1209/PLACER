#include "te_family_alignment.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

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

template <typename Fn>
void for_each_valid_kmer(const std::string& seq, int32_t k, Fn&& fn) {
    if (k <= 0 || static_cast<int32_t>(seq.size()) < k) {
        return;
    }

    const uint64_t mask = (k >= 32)
        ? ~uint64_t{0}
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

std::string upper_ascii(std::string s) {
    for (char& c : s) {
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    }
    return s;
}

struct ParsedTeName {
    std::string exact_name = "NA";
    std::string family = "NA";
    std::string family_key = "NA";
};

ParsedTeName parse_te_name(const std::string& te_name) {
    ParsedTeName parts;
    const std::string token = take_header_token(te_name);
    if (token.empty()) {
        return parts;
    }

    std::string family;
    std::string exact = token;
    const size_t hash = token.find('#');
    if (hash != std::string::npos) {
        const std::string left = token.substr(0, hash);
        const std::string right = token.substr(hash + 1);
        exact = left.empty() ? token : left;
        family = right.empty() ? exact : right;
        const size_t slash = family.rfind('/');
        if (slash != std::string::npos && (slash + 1) < family.size()) {
            family = family.substr(slash + 1);
        }
    } else {
        const size_t colon = token.find(':');
        if (colon != std::string::npos) {
            family = token.substr(0, colon);
            exact = ((colon + 1) < token.size()) ? token.substr(colon + 1) : token;
        } else {
            family = token;
        }
    }

    const std::string family_key = upper_ascii(family);
    if (family_key.rfind("ALU", 0) == 0) {
        parts.family = "ALU";
        parts.family_key = "ALU";
    } else if (family_key.rfind("L1", 0) == 0) {
        parts.family = "L1";
        parts.family_key = "L1";
    } else if (family_key.rfind("SVA", 0) == 0) {
        parts.family = "SVA";
        parts.family_key = "SVA";
    } else if (family_key.rfind("HERV", 0) == 0) {
        parts.family = "HERV";
        parts.family_key = "HERV";
    } else {
        parts.family = family.empty() ? exact : family;
        parts.family_key = family_key.empty() ? upper_ascii(exact) : family_key;
    }

    parts.exact_name = exact.empty() ? token : exact;
    return parts;
}

struct TeRepresentativeSketch {
    int32_t te_id = -1;
    std::vector<uint64_t> mins;
};

std::vector<uint64_t> build_family_sketch(
    const std::string& seq,
    int32_t k) {
    std::vector<uint64_t> mins;
    mins.reserve(64);
    for_each_valid_kmer(seq, k, [&](int32_t, uint64_t key) {
        mins.push_back(key);
    });
    std::sort(mins.begin(), mins.end());
    mins.erase(std::unique(mins.begin(), mins.end()), mins.end());
    if (mins.size() > 64) {
        mins.resize(64);
    }
    return mins;
}

double sketch_distance(
    const std::vector<uint64_t>& lhs,
    const std::vector<uint64_t>& rhs) {
    size_t i = 0;
    size_t j = 0;
    int32_t intersect = 0;
    while (i < lhs.size() && j < rhs.size()) {
        if (lhs[i] == rhs[j]) {
            ++intersect;
            ++i;
            ++j;
        } else if (lhs[i] < rhs[j]) {
            ++i;
        } else {
            ++j;
        }
    }
    const int32_t denom = std::max<int32_t>(
        1,
        static_cast<int32_t>(lhs.size() + rhs.size() - intersect));
    return 1.0 - (static_cast<double>(intersect) / static_cast<double>(denom));
}

std::vector<int32_t> pick_family_representative_ids(
    const std::vector<TeRepresentativeSketch>& sketches,
    int32_t rep_limit) {
    std::vector<int32_t> out;
    if (sketches.empty() || rep_limit <= 0) {
        return out;
    }

    size_t medoid_index = 0;
    double best_total_distance = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < sketches.size(); ++i) {
        double total_distance = 0.0;
        for (size_t j = 0; j < sketches.size(); ++j) {
            if (i == j) {
                continue;
            }
            total_distance += sketch_distance(sketches[i].mins, sketches[j].mins);
        }
        if (total_distance < best_total_distance) {
            best_total_distance = total_distance;
            medoid_index = i;
        }
    }
    out.push_back(sketches[medoid_index].te_id);

    while (static_cast<int32_t>(out.size()) < rep_limit &&
           static_cast<int32_t>(out.size()) < static_cast<int32_t>(sketches.size())) {
        int32_t best_next_id = -1;
        double best_min_distance = -1.0;
        for (const auto& sketch : sketches) {
            if (std::find(out.begin(), out.end(), sketch.te_id) != out.end()) {
                continue;
            }

            double min_distance_to_selected = std::numeric_limits<double>::infinity();
            for (int32_t selected_id : out) {
                const auto selected_it = std::find_if(
                    sketches.begin(),
                    sketches.end(),
                    [&](const TeRepresentativeSketch& item) {
                        return item.te_id == selected_id;
                    });
                if (selected_it == sketches.end()) {
                    continue;
                }
                min_distance_to_selected = std::min(
                    min_distance_to_selected,
                    sketch_distance(sketch.mins, selected_it->mins));
            }

            if (min_distance_to_selected > best_min_distance) {
                best_min_distance = min_distance_to_selected;
                best_next_id = sketch.te_id;
            }
        }

        if (best_next_id < 0) {
            break;
        }
        out.push_back(best_next_id);
    }

    return out;
}

double weight_from_family_count(int32_t family_count) {
    if (family_count <= 1) {
        return 1.0;
    }
    return 1.0 / static_cast<double>(family_count);
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

struct SeedHit {
    int32_t query_pos = 0;
    int32_t target_pos = 0;
    double weight = 0.0;
};

struct ChainSummary {
    double score = 0.0;
    double query_coverage = 0.0;
};

ChainSummary best_chain_summary(
    std::vector<SeedHit> hits,
    int32_t query_len,
    int32_t kmer_size) {
    ChainSummary summary;
    if (hits.empty() || query_len <= 0 || kmer_size <= 0) {
        return summary;
    }

    std::sort(hits.begin(), hits.end(), [](const SeedHit& lhs, const SeedHit& rhs) {
        if (lhs.query_pos != rhs.query_pos) {
            return lhs.query_pos < rhs.query_pos;
        }
        return lhs.target_pos < rhs.target_pos;
    });

    std::vector<double> dp(hits.size(), 0.0);
    std::vector<int32_t> prev(hits.size(), -1);
    size_t best_index = 0;
    double best_score = 0.0;

    for (size_t i = 0; i < hits.size(); ++i) {
        dp[i] = hits[i].weight;
        for (size_t j = 0; j < i; ++j) {
            const int32_t q_gap = hits[i].query_pos - hits[j].query_pos;
            const int32_t t_gap = hits[i].target_pos - hits[j].target_pos;
            if (q_gap <= 0 || t_gap <= 0) {
                continue;
            }
            const int32_t diag_delta = std::abs(
                (hits[i].target_pos - hits[i].query_pos) -
                (hits[j].target_pos - hits[j].query_pos));
            const double transition_penalty =
                (0.01 * static_cast<double>(std::abs(q_gap - t_gap))) +
                (0.005 * static_cast<double>(diag_delta));
            const double candidate = dp[j] + hits[i].weight - transition_penalty;
            if (candidate > dp[i]) {
                dp[i] = candidate;
                prev[i] = static_cast<int32_t>(j);
            }
        }
        if (dp[i] > best_score) {
            best_score = dp[i];
            best_index = i;
        }
    }

    std::vector<std::pair<int32_t, int32_t>> chain_intervals;
    for (int32_t idx = static_cast<int32_t>(best_index); idx >= 0; idx = prev[static_cast<size_t>(idx)]) {
        chain_intervals.push_back({
            hits[static_cast<size_t>(idx)].query_pos,
            hits[static_cast<size_t>(idx)].query_pos + kmer_size,
        });
    }
    std::sort(chain_intervals.begin(), chain_intervals.end());

    int32_t covered = 0;
    int32_t cur_start = -1;
    int32_t cur_end = -1;
    for (const auto& interval : chain_intervals) {
        if (cur_start < 0) {
            cur_start = interval.first;
            cur_end = interval.second;
            continue;
        }
        if (interval.first > cur_end) {
            covered += (cur_end - cur_start);
            cur_start = interval.first;
            cur_end = interval.second;
            continue;
        }
        cur_end = std::max(cur_end, interval.second);
    }
    if (cur_start >= 0) {
        covered += (cur_end - cur_start);
    }

    summary.score = std::max(0.0, best_score);
    summary.query_coverage = std::clamp(
        static_cast<double>(covered) / static_cast<double>(std::max(1, query_len)),
        0.0,
        1.0);
    return summary;
}

}  // namespace

TeFamilyCacheBundle build_te_family_cache(
    const std::vector<std::string>& te_names,
    const std::vector<std::string>& te_sequences,
    const std::vector<std::string>& te_reverse_complement_sequences,
    int32_t family_representatives,
    int32_t kmer_size) {
    TeFamilyCacheBundle bundle;
    auto alignment_index = std::make_shared<TeFamilyAlignmentIndex>();
    auto rep_groups = std::make_shared<TeFamilyGroupCache>();
    const int32_t rep_limit = std::max(1, family_representatives);
    alignment_index->kmer_size = std::max(7, kmer_size);

    std::unordered_map<std::string, size_t> family_to_group;
    family_to_group.reserve(te_names.size());
    std::vector<ParsedTeName> parsed_names;
    parsed_names.reserve(te_names.size());
    for (const auto& name : te_names) {
        parsed_names.push_back(parse_te_name(name));
    }

    for (size_t i = 0; i < te_names.size(); ++i) {
        const ParsedTeName& parsed = parsed_names[i];
        auto it = family_to_group.find(parsed.family_key);
        if (it == family_to_group.end()) {
            TeFamilyGroup group;
            group.family = parsed.family;
            group.family_key = parsed.family_key;
            alignment_index->groups.push_back(std::move(group));
            const size_t group_index = alignment_index->groups.size() - 1;
            family_to_group.emplace(parsed.family_key, group_index);
            alignment_index->family_to_group.emplace(parsed.family_key, group_index);
            it = family_to_group.find(parsed.family_key);
        }

        TeFamilyGroup& group = alignment_index->groups[it->second];
        group.member_te_ids.push_back(static_cast<int32_t>(i));
    }

    const int32_t sketch_k = std::max(7, alignment_index->kmer_size);
    for (auto& group : alignment_index->groups) {
        std::vector<TeRepresentativeSketch> sketches;
        sketches.reserve(group.member_te_ids.size());
        for (int32_t te_id : group.member_te_ids) {
            if (te_id < 0 ||
                te_id >= static_cast<int32_t>(te_sequences.size()) ||
                te_id >= static_cast<int32_t>(te_reverse_complement_sequences.size())) {
                continue;
            }
            sketches.push_back({
                te_id,
                build_family_sketch(te_sequences[static_cast<size_t>(te_id)], sketch_k),
            });
        }

        const std::vector<int32_t> representative_ids =
            pick_family_representative_ids(sketches, rep_limit);
        group.representatives.clear();
        group.representatives.reserve(representative_ids.size());
        for (int32_t te_id : representative_ids) {
            if (te_id < 0 ||
                te_id >= static_cast<int32_t>(te_sequences.size()) ||
                te_id >= static_cast<int32_t>(te_reverse_complement_sequences.size()) ||
                te_id >= static_cast<int32_t>(parsed_names.size())) {
                continue;
            }
            const ParsedTeName& parsed = parsed_names[static_cast<size_t>(te_id)];
            TeFamilyRepresentative rep;
            rep.te_id = te_id;
            rep.family = parsed.family;
            rep.family_key = parsed.family_key;
            rep.exact_name = parsed.exact_name;
            rep.sequence = te_sequences[static_cast<size_t>(te_id)];
            rep.reverse_complement_sequence =
                te_reverse_complement_sequences[static_cast<size_t>(te_id)];
            populate_kmer_positions(
                rep.sequence,
                alignment_index->kmer_size,
                rep.forward_kmer_positions);
            populate_kmer_positions(
                rep.reverse_complement_sequence,
                alignment_index->kmer_size,
                rep.reverse_kmer_positions);
            group.representatives.push_back(std::move(rep));
        }
    }

    for (const auto& group : alignment_index->groups) {
        std::unordered_set<uint64_t> family_keys;
        for (const auto& rep : group.representatives) {
            for (const auto& kv : rep.forward_kmer_positions) {
                family_keys.insert(kv.first);
            }
            for (const auto& kv : rep.reverse_kmer_positions) {
                family_keys.insert(kv.first);
            }
        }
        for (uint64_t key : family_keys) {
            alignment_index->kmer_family_counts[key] += 1;
        }
    }

    std::vector<std::pair<uint64_t, int32_t>> weighted_keys;
    weighted_keys.reserve(alignment_index->kmer_family_counts.size());
    for (const auto& kv : alignment_index->kmer_family_counts) {
        alignment_index->kmer_family_weights[kv.first] = weight_from_family_count(kv.second);
        weighted_keys.push_back(kv);
    }
    std::sort(weighted_keys.begin(), weighted_keys.end(), [](const auto& lhs, const auto& rhs) {
        if (lhs.second != rhs.second) {
            return lhs.second > rhs.second;
        }
        return lhs.first < rhs.first;
    });
    alignment_index->representative_seed_keys.reserve(weighted_keys.size());
    for (const auto& kv : weighted_keys) {
        alignment_index->representative_seed_keys.push_back(kv.first);
    }

    rep_groups->groups = alignment_index->groups;
    bundle.alignment_index = std::move(alignment_index);
    bundle.rep_groups = std::move(rep_groups);
    return bundle;
}

std::vector<TeFamilyChainScore> score_family_prefilter(
    const TeFamilyAlignmentIndex& index,
    const std::string& query,
    int32_t kmer_size,
    int32_t family_topn) {
    std::vector<TeFamilyChainScore> out;
    const int32_t k = std::max(7, kmer_size);
    if (query.empty() || index.groups.empty() || static_cast<int32_t>(query.size()) < k) {
        return out;
    }

    std::unordered_map<uint64_t, std::vector<int32_t>> query_positions;
    query_positions.reserve(query.size());
    double total_query_weight = 0.0;
    for_each_valid_kmer(query, k, [&](int32_t start, uint64_t key) {
        query_positions[key].push_back(start);
        const auto weight_it = index.kmer_family_weights.find(key);
        total_query_weight +=
            (weight_it == index.kmer_family_weights.end()) ? 1.0 : weight_it->second;
    });
    if (query_positions.empty() || total_query_weight <= 0.0) {
        return out;
    }

    out.reserve(index.groups.size());
    for (size_t group_index = 0; group_index < index.groups.size(); ++group_index) {
        const auto& group = index.groups[group_index];
        std::vector<std::tuple<double, double, int32_t, bool>> rep_scores;
        rep_scores.reserve(group.representatives.size() * 2);

        for (size_t rep_index = 0; rep_index < group.representatives.size(); ++rep_index) {
            const auto& rep = group.representatives[rep_index];

            auto score_orientation = [&](const auto& kmer_positions, bool reverse) {
                std::vector<SeedHit> hits;
                for (const auto& qkv : query_positions) {
                    const auto it = kmer_positions.find(qkv.first);
                    if (it == kmer_positions.end()) {
                        continue;
                    }
                    const auto weight_it = index.kmer_family_weights.find(qkv.first);
                    const double weight =
                        (weight_it == index.kmer_family_weights.end()) ? 1.0 : weight_it->second;
                    for (int32_t qpos : qkv.second) {
                        for (int32_t tpos : it->second) {
                            hits.push_back({qpos, tpos, weight});
                        }
                    }
                }
                const ChainSummary summary = best_chain_summary(
                    std::move(hits),
                    static_cast<int32_t>(query.size()),
                    k);
                const double score_norm = summary.score / total_query_weight;
                rep_scores.push_back({score_norm, summary.query_coverage, static_cast<int32_t>(rep_index), reverse});
            };

            score_orientation(rep.forward_kmer_positions, false);
            score_orientation(rep.reverse_kmer_positions, true);
        }

        if (rep_scores.empty()) {
            continue;
        }

        std::sort(rep_scores.begin(), rep_scores.end(), [](const auto& lhs, const auto& rhs) {
            if (std::get<0>(lhs) != std::get<0>(rhs)) {
                return std::get<0>(lhs) > std::get<0>(rhs);
            }
            return std::get<1>(lhs) > std::get<1>(rhs);
        });

        TeFamilyChainScore score;
        score.group_index = static_cast<int32_t>(group_index);
        score.representative_index = std::get<2>(rep_scores.front());
        score.reverse = std::get<3>(rep_scores.front());
        score.best_rep_chain_norm = std::get<0>(rep_scores.front());
        score.chained_query_coverage = std::get<1>(rep_scores.front());
        if (rep_scores.size() >= 2) {
            score.second_rep_chain_norm = std::get<0>(rep_scores[1]);
        }
        score.family_prefilter_score =
            (0.60 * score.best_rep_chain_norm) +
            (0.20 * score.second_rep_chain_norm) +
            (0.20 * score.chained_query_coverage);
        out.push_back(score);
    }

    std::sort(out.begin(), out.end(), [](const TeFamilyChainScore& lhs, const TeFamilyChainScore& rhs) {
        if (lhs.family_prefilter_score != rhs.family_prefilter_score) {
            return lhs.family_prefilter_score > rhs.family_prefilter_score;
        }
        return lhs.group_index < rhs.group_index;
    });

    const int32_t keep = std::max(1, family_topn);
    if (static_cast<int32_t>(out.size()) > keep) {
        out.resize(static_cast<size_t>(keep));
    }
    return out;
}

TemplateBandEstimate estimate_template_band_from_positions(
    const KmerPositionMap& query_positions,
    int32_t query_len,
    const KmerPositionMap& target_positions,
    int32_t target_len,
    int32_t kmer_size) {
    TemplateBandEstimate estimate;
    estimate.band_center = target_len - query_len;
    estimate.band_width = std::max(32, std::abs(target_len - query_len) + 32);

    const int32_t k = std::max(7, kmer_size);
    if (query_len < k || target_len < k ||
        query_positions.empty() || target_positions.empty()) {
        return estimate;
    }

    struct SeedDiagonalHit {
        int32_t diag = 0;
        int32_t query_start = 0;
        int32_t query_end = 0;
    };

    std::vector<SeedDiagonalHit> hits;
    hits.reserve(static_cast<size_t>(query_len));
    constexpr size_t kMaxSeedMultiplicity = 4;
    for (const auto& qkv : query_positions) {
        const auto it = target_positions.find(qkv.first);
        if (it == target_positions.end()) {
            continue;
        }
        if (qkv.second.size() > kMaxSeedMultiplicity ||
            it->second.size() > kMaxSeedMultiplicity) {
            continue;
        }
        for (int32_t qpos : qkv.second) {
            for (int32_t tpos : it->second) {
                hits.push_back({tpos - qpos, qpos, qpos + k});
            }
        }
    }
    if (hits.empty()) {
        return estimate;
    }

    constexpr int32_t kDiagBinSize = 8;
    std::unordered_map<int32_t, int32_t> bin_counts;
    bin_counts.reserve(hits.size());
    int32_t best_bin = 0;
    int32_t best_bin_count = 0;
    for (const auto& hit : hits) {
        const int32_t bin = hit.diag / kDiagBinSize;
        const int32_t count = ++bin_counts[bin];
        if (count > best_bin_count) {
            best_bin_count = count;
            best_bin = bin;
        }
    }

    std::vector<SeedDiagonalHit> cluster_hits;
    cluster_hits.reserve(hits.size());
    const int32_t best_bin_center = best_bin * kDiagBinSize;
    for (const auto& hit : hits) {
        if (std::abs(hit.diag - best_bin_center) <= kDiagBinSize) {
            cluster_hits.push_back(hit);
        }
    }
    if (cluster_hits.empty()) {
        cluster_hits = hits;
    }

    std::vector<int32_t> diags;
    diags.reserve(cluster_hits.size());
    std::vector<std::pair<int32_t, int32_t>> intervals;
    intervals.reserve(cluster_hits.size());
    int32_t min_diag = std::numeric_limits<int32_t>::max();
    int32_t max_diag = std::numeric_limits<int32_t>::min();
    for (const auto& hit : cluster_hits) {
        diags.push_back(hit.diag);
        intervals.push_back({hit.query_start, hit.query_end});
        min_diag = std::min(min_diag, hit.diag);
        max_diag = std::max(max_diag, hit.diag);
    }
    std::sort(diags.begin(), diags.end());
    estimate.band_center = diags[diags.size() / 2];

    std::sort(intervals.begin(), intervals.end());
    int32_t covered = 0;
    int32_t cur_start = -1;
    int32_t cur_end = -1;
    for (const auto& interval : intervals) {
        if (cur_start < 0) {
            cur_start = interval.first;
            cur_end = interval.second;
            continue;
        }
        if (interval.first > cur_end) {
            covered += (cur_end - cur_start);
            cur_start = interval.first;
            cur_end = interval.second;
            continue;
        }
        cur_end = std::max(cur_end, interval.second);
    }
    if (cur_start >= 0) {
        covered += (cur_end - cur_start);
    }

    estimate.query_coverage = std::clamp(
        static_cast<double>(covered) / static_cast<double>(std::max(1, query_len)),
        0.0,
        1.0);
    estimate.seed_score = estimate.query_coverage;
    estimate.band_width = std::max({
        24,
        std::abs(max_diag - min_diag) + k + 8,
        std::abs(target_len - query_len) + 16,
    });
    estimate.band_width = std::min(
        std::max(query_len, target_len),
        estimate.band_width);
    return estimate;
}

TemplateBandEstimate estimate_template_band(
    const std::string& query,
    const std::string& target,
    int32_t kmer_size) {
    const int32_t query_len = static_cast<int32_t>(query.size());
    const int32_t target_len = static_cast<int32_t>(target.size());
    const int32_t k = std::max(7, kmer_size);

    KmerPositionMap query_positions;
    KmerPositionMap target_positions;
    if (query_len >= k) {
        query_positions.reserve(query.size());
        populate_kmer_positions(query, k, query_positions);
    }
    if (target_len >= k) {
        target_positions.reserve(target.size());
        populate_kmer_positions(target, k, target_positions);
    }

    return estimate_template_band_from_positions(
        query_positions,
        query_len,
        target_positions,
        target_len,
        k);
}

namespace {

struct AffineAlignmentCell {
    int32_t score = std::numeric_limits<int32_t>::min() / 4;
    int32_t matches = 0;
    int32_t alignment_len = 0;
    int32_t query_bases = 0;
    bool valid = false;
};

bool better_affine_alignment_cell(
    const AffineAlignmentCell& lhs,
    const AffineAlignmentCell& rhs,
    int32_t query_len) {
    if (lhs.valid != rhs.valid) {
        return lhs.valid;
    }
    if (!lhs.valid) {
        return false;
    }
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

AffineAlignmentCell extend_diag(
    const AffineAlignmentCell& prev,
    bool match) {
    if (!prev.valid) {
        return AffineAlignmentCell{};
    }
    AffineAlignmentCell next = prev;
    next.score += match ? 2 : -2;
    next.alignment_len += 1;
    next.query_bases += 1;
    if (match) {
        next.matches += 1;
    }
    return next;
}

AffineAlignmentCell extend_up(
    const AffineAlignmentCell& prev,
    int32_t penalty) {
    if (!prev.valid) {
        return AffineAlignmentCell{};
    }
    AffineAlignmentCell next = prev;
    next.score -= penalty;
    next.alignment_len += 1;
    next.query_bases += 1;
    return next;
}

AffineAlignmentCell extend_left(
    const AffineAlignmentCell& prev,
    int32_t penalty) {
    if (!prev.valid) {
        return AffineAlignmentCell{};
    }
    AffineAlignmentCell next = prev;
    next.score -= penalty;
    next.alignment_len += 1;
    return next;
}

AffineAlignmentCell start_local_diag(bool match) {
    if (!match) {
        return AffineAlignmentCell{};
    }

    AffineAlignmentCell cell;
    cell.score = 2;
    cell.matches = 1;
    cell.alignment_len = 1;
    cell.query_bases = 1;
    cell.valid = true;
    return cell;
}

AffineAlignmentCell reset_if_nonpositive(AffineAlignmentCell cell) {
    if (!cell.valid || cell.score <= 0) {
        return AffineAlignmentCell{};
    }
    return cell;
}

}  // namespace

ExactAlignmentSummary banded_semiglobal_affine_align(
    const std::string& query,
    const std::string& target,
    int32_t band_center,
    int32_t band_width) {
    ExactAlignmentSummary summary;
    const int32_t n = static_cast<int32_t>(query.size());
    const int32_t m = static_cast<int32_t>(target.size());
    if (n <= 0 || m <= 0) {
        return summary;
    }

    const int32_t width = std::max(8, band_width);
    constexpr int32_t kGapOpenPenalty = 3;
    constexpr int32_t kGapExtendPenalty = 1;

    std::vector<AffineAlignmentCell> match_prev(static_cast<size_t>(m + 1));
    std::vector<AffineAlignmentCell> ins_prev(static_cast<size_t>(m + 1));
    std::vector<AffineAlignmentCell> del_prev(static_cast<size_t>(m + 1));
    std::vector<AffineAlignmentCell> match_curr(static_cast<size_t>(m + 1));
    std::vector<AffineAlignmentCell> ins_curr(static_cast<size_t>(m + 1));
    std::vector<AffineAlignmentCell> del_curr(static_cast<size_t>(m + 1));

    AffineAlignmentCell best;

    for (int32_t i = 1; i <= n; ++i) {
        std::fill(match_curr.begin(), match_curr.end(), AffineAlignmentCell{});
        std::fill(ins_curr.begin(), ins_curr.end(), AffineAlignmentCell{});
        std::fill(del_curr.begin(), del_curr.end(), AffineAlignmentCell{});

        const int32_t row_center = i + band_center;
        const int32_t j_min = std::max(0, row_center - width);
        const int32_t j_max = std::min(m, row_center + width);
        for (int32_t j = j_min; j <= j_max; ++j) {
            const size_t idx = static_cast<size_t>(j);

            if (j == 0) {
                continue;
            }

            const size_t prev_idx = static_cast<size_t>(j - 1);
            const bool match =
                query[static_cast<size_t>(i - 1)] == target[static_cast<size_t>(j - 1)];
            AffineAlignmentCell best_match = start_local_diag(match);
            const AffineAlignmentCell match_from_match = reset_if_nonpositive(extend_diag(
                match_prev[prev_idx],
                match));
            if (better_affine_alignment_cell(match_from_match, best_match, n)) {
                best_match = match_from_match;
            }
            const AffineAlignmentCell match_from_ins = reset_if_nonpositive(extend_diag(
                ins_prev[prev_idx],
                match));
            if (better_affine_alignment_cell(match_from_ins, best_match, n)) {
                best_match = match_from_ins;
            }
            const AffineAlignmentCell match_from_del = reset_if_nonpositive(extend_diag(
                del_prev[prev_idx],
                match));
            if (better_affine_alignment_cell(match_from_del, best_match, n)) {
                best_match = match_from_del;
            }
            match_curr[idx] = best_match;

            AffineAlignmentCell best_ins = reset_if_nonpositive(
                extend_up(match_prev[idx], kGapOpenPenalty));
            const AffineAlignmentCell ins_from_ins = reset_if_nonpositive(
                extend_up(ins_prev[idx], kGapExtendPenalty));
            if (better_affine_alignment_cell(ins_from_ins, best_ins, n)) {
                best_ins = ins_from_ins;
            }
            const AffineAlignmentCell ins_from_del = reset_if_nonpositive(
                extend_up(del_prev[idx], kGapOpenPenalty));
            if (better_affine_alignment_cell(ins_from_del, best_ins, n)) {
                best_ins = ins_from_del;
            }
            ins_curr[idx] = best_ins;

            AffineAlignmentCell best_del = reset_if_nonpositive(
                extend_left(match_curr[prev_idx], kGapOpenPenalty));
            const AffineAlignmentCell del_from_del = reset_if_nonpositive(
                extend_left(del_curr[prev_idx], kGapExtendPenalty));
            if (better_affine_alignment_cell(del_from_del, best_del, n)) {
                best_del = del_from_del;
            }
            const AffineAlignmentCell del_from_ins = reset_if_nonpositive(
                extend_left(ins_curr[prev_idx], kGapOpenPenalty));
            if (better_affine_alignment_cell(del_from_ins, best_del, n)) {
                best_del = del_from_ins;
            }
            del_curr[idx] = best_del;

            if (better_affine_alignment_cell(match_curr[idx], best, n)) {
                best = match_curr[idx];
            }
            if (better_affine_alignment_cell(ins_curr[idx], best, n)) {
                best = ins_curr[idx];
            }
            if (better_affine_alignment_cell(del_curr[idx], best, n)) {
                best = del_curr[idx];
            }
        }

        std::swap(match_prev, match_curr);
        std::swap(ins_prev, ins_curr);
        std::swap(del_prev, del_curr);
    }

    if (!best.valid || best.alignment_len <= 0 || best.query_bases <= 0) {
        return summary;
    }
    summary.identity = static_cast<double>(best.matches) /
        static_cast<double>(best.alignment_len);
    summary.query_coverage = static_cast<double>(best.query_bases) /
        static_cast<double>(std::max(1, n));
    summary.score = summary.identity * summary.query_coverage;
    return summary;
}

}  // namespace placer
