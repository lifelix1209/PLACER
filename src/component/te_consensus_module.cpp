#include "pipeline.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace placer {
namespace {

constexpr int32_t kInf = std::numeric_limits<int32_t>::max() / 8;
constexpr double kMadScale = 1.4826;

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

std::string upper_acgt(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        out.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
    }
    return out;
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

double median(std::vector<double> values) {
    if (values.empty()) {
        return 0.0;
    }
    std::sort(values.begin(), values.end());
    const size_t mid = values.size() / 2;
    return (values.size() % 2 == 1)
        ? values[mid]
        : 0.5 * (values[mid - 1] + values[mid]);
}

double weighted_median(std::vector<std::pair<double, double>> value_weights) {
    value_weights.erase(
        std::remove_if(value_weights.begin(), value_weights.end(), [](const auto& row) {
            return row.second <= 0.0;
        }),
        value_weights.end());
    if (value_weights.empty()) {
        return 0.0;
    }
    std::sort(value_weights.begin(), value_weights.end(), [](const auto& a, const auto& b) {
        if (a.first != b.first) {
            return a.first < b.first;
        }
        return a.second > b.second;
    });

    double total = 0.0;
    for (const auto& row : value_weights) {
        total += row.second;
    }
    if (total <= 0.0) {
        return value_weights.front().first;
    }

    const double half = 0.5 * total;
    double acc = 0.0;
    for (const auto& row : value_weights) {
        acc += row.second;
        if (acc >= half) {
            return row.first;
        }
    }
    return value_weights.back().first;
}

double mad_scaled(const std::vector<double>& values, double center) {
    if (values.empty()) {
        return 0.0;
    }
    std::vector<double> abs_dev;
    abs_dev.reserve(values.size());
    for (double v : values) {
        abs_dev.push_back(std::abs(v - center));
    }
    return kMadScale * median(std::move(abs_dev));
}

double relative_diff(double a, double b) {
    const double denom = std::max({1e-9, std::abs(a), std::abs(b)});
    return std::abs(a - b) / denom;
}

int32_t safe_floor_to_i32(double v) {
    if (v < static_cast<double>(std::numeric_limits<int32_t>::min())) {
        return std::numeric_limits<int32_t>::min();
    }
    if (v > static_cast<double>(std::numeric_limits<int32_t>::max())) {
        return std::numeric_limits<int32_t>::max();
    }
    return static_cast<int32_t>(std::floor(v));
}

int32_t safe_ceil_to_i32(double v) {
    if (v < static_cast<double>(std::numeric_limits<int32_t>::min())) {
        return std::numeric_limits<int32_t>::min();
    }
    if (v > static_cast<double>(std::numeric_limits<int32_t>::max())) {
        return std::numeric_limits<int32_t>::max();
    }
    return static_cast<int32_t>(std::ceil(v));
}

struct TeRecord {
    std::string name;
    std::string sequence;
    std::unordered_map<uint64_t, std::vector<int32_t>> kmer_to_positions;
};

struct Placement {
    int32_t te_start = 0;
    int32_t te_end = 0;
    double norm_cost = 1.0;
    double span_start = 0.0;
};

struct FragmentPlacementState {
    const InsertionFragment* fragment = nullptr;
    const FragmentTEHit* hit = nullptr;
    std::vector<Placement> placements;
    size_t best_idx = 0;
    size_t chosen_idx = 0;
    double delta_tpl = 0.0;
    double start_span_best = 0.0;
    bool core_candidate = false;
    bool split_sa = false;
};

struct SeedStarts {
    std::vector<int32_t> starts;
};

struct AffineAlignResult {
    bool ok = false;
    int32_t mismatch_cnt = 0;
    int32_t gap_open_cnt = 0;
    int32_t gap_extend_cnt = 0;
    int32_t aligned_cols = 0;
    int32_t target_consumed = 0;
    double norm_cost = 1.0;
};

inline size_t mat_index(int32_t i, int32_t j, int32_t cols) {
    return static_cast<size_t>(i) * static_cast<size_t>(cols) + static_cast<size_t>(j);
}

AffineAlignResult semiglobal_affine_align_suffix(
    const std::string& query,
    const std::string& target_suffix) {
    AffineAlignResult out;
    const int32_t n = static_cast<int32_t>(query.size());
    const int32_t m = static_cast<int32_t>(target_suffix.size());
    if (n <= 0 || m <= 0) {
        return out;
    }

    const int32_t rows = n + 1;
    const int32_t cols = m + 1;
    const size_t total = static_cast<size_t>(rows) * static_cast<size_t>(cols);

    std::vector<int32_t> m_cost(total, kInf);
    std::vector<int32_t> i_cost(total, kInf);
    std::vector<int32_t> d_cost(total, kInf);

    std::vector<uint8_t> m_prev(total, 255);
    std::vector<uint8_t> i_prev(total, 255);
    std::vector<uint8_t> d_prev(total, 255);

    const int32_t gap_open = 2;
    const int32_t gap_extend = 1;

    m_cost[mat_index(0, 0, cols)] = 0;

    for (int32_t i = 1; i <= n; ++i) {
        const size_t idx = mat_index(i, 0, cols);
        if (i == 1) {
            i_cost[idx] = gap_open + gap_extend;
            i_prev[idx] = 0;
        } else {
            const size_t up = mat_index(i - 1, 0, cols);
            if (i_cost[up] < kInf) {
                i_cost[idx] = i_cost[up] + gap_extend;
                i_prev[idx] = 1;
            }
        }
    }

    for (int32_t i = 1; i <= n; ++i) {
        for (int32_t j = 1; j <= m; ++j) {
            const size_t idx = mat_index(i, j, cols);
            const size_t diag = mat_index(i - 1, j - 1, cols);

            const int32_t sub = (query[static_cast<size_t>(i - 1)] == target_suffix[static_cast<size_t>(j - 1)]) ? 0 : 1;

            int32_t best_m = kInf;
            uint8_t best_m_prev = 255;
            if (m_cost[diag] < best_m) {
                best_m = m_cost[diag];
                best_m_prev = 0;
            }
            if (i_cost[diag] < best_m) {
                best_m = i_cost[diag];
                best_m_prev = 1;
            }
            if (d_cost[diag] < best_m) {
                best_m = d_cost[diag];
                best_m_prev = 2;
            }
            if (best_m < kInf) {
                m_cost[idx] = best_m + sub;
                m_prev[idx] = best_m_prev;
            }

            const size_t up = mat_index(i - 1, j, cols);
            int32_t best_i = kInf;
            uint8_t best_i_prev = 255;

            if (m_cost[up] < kInf && m_cost[up] + gap_open + gap_extend < best_i) {
                best_i = m_cost[up] + gap_open + gap_extend;
                best_i_prev = 0;
            }
            if (i_cost[up] < kInf && i_cost[up] + gap_extend < best_i) {
                best_i = i_cost[up] + gap_extend;
                best_i_prev = 1;
            }
            if (d_cost[up] < kInf && d_cost[up] + gap_open + gap_extend < best_i) {
                best_i = d_cost[up] + gap_open + gap_extend;
                best_i_prev = 2;
            }
            if (best_i < kInf) {
                i_cost[idx] = best_i;
                i_prev[idx] = best_i_prev;
            }

            const size_t left = mat_index(i, j - 1, cols);
            int32_t best_d = kInf;
            uint8_t best_d_prev = 255;
            if (m_cost[left] < kInf && m_cost[left] + gap_open + gap_extend < best_d) {
                best_d = m_cost[left] + gap_open + gap_extend;
                best_d_prev = 0;
            }
            if (d_cost[left] < kInf && d_cost[left] + gap_extend < best_d) {
                best_d = d_cost[left] + gap_extend;
                best_d_prev = 2;
            }
            if (i_cost[left] < kInf && i_cost[left] + gap_open + gap_extend < best_d) {
                best_d = i_cost[left] + gap_open + gap_extend;
                best_d_prev = 1;
            }
            if (best_d < kInf) {
                d_cost[idx] = best_d;
                d_prev[idx] = best_d_prev;
            }
        }
    }

    int32_t best_j = 0;
    uint8_t best_state = 0;
    int32_t best_cost = kInf;
    for (int32_t j = 0; j <= m; ++j) {
        const size_t idx = mat_index(n, j, cols);
        const int32_t c_m = m_cost[idx];
        const int32_t c_i = i_cost[idx];
        const int32_t c_d = d_cost[idx];

        auto better = [&](int32_t cost, uint8_t state) {
            if (cost > best_cost) {
                return false;
            }
            if (cost < best_cost) {
                return true;
            }
            if (j > best_j) {
                return true;
            }
            if (j < best_j) {
                return false;
            }
            return state < best_state;
        };

        if (better(c_m, 0)) {
            best_cost = c_m;
            best_j = j;
            best_state = 0;
        }
        if (better(c_i, 1)) {
            best_cost = c_i;
            best_j = j;
            best_state = 1;
        }
        if (better(c_d, 2)) {
            best_cost = c_d;
            best_j = j;
            best_state = 2;
        }
    }

    if (best_cost >= kInf) {
        return out;
    }

    int32_t i = n;
    int32_t j = best_j;
    uint8_t state = best_state;

    int32_t mismatch_cnt = 0;
    int32_t gap_open_cnt = 0;
    int32_t gap_extend_cnt = 0;
    int32_t aligned_cols = 0;
    int32_t target_consumed = 0;

    while (i > 0 || j > 0) {
        const size_t idx = mat_index(i, j, cols);
        if (state == 0) {
            if (i <= 0 || j <= 0 || m_prev[idx] == 255) {
                return out;
            }
            const uint8_t prev = m_prev[idx];
            if (query[static_cast<size_t>(i - 1)] != target_suffix[static_cast<size_t>(j - 1)]) {
                ++mismatch_cnt;
            }
            ++aligned_cols;
            ++target_consumed;
            --i;
            --j;
            state = prev;
        } else if (state == 1) {
            if (i <= 0 || i_prev[idx] == 255) {
                return out;
            }
            const uint8_t prev = i_prev[idx];
            ++aligned_cols;
            if (prev == 1) {
                ++gap_extend_cnt;
            } else {
                ++gap_open_cnt;
            }
            --i;
            state = prev;
        } else if (state == 2) {
            if (j <= 0 || d_prev[idx] == 255) {
                return out;
            }
            const uint8_t prev = d_prev[idx];
            ++aligned_cols;
            ++target_consumed;
            if (prev == 2) {
                ++gap_extend_cnt;
            } else {
                ++gap_open_cnt;
            }
            --j;
            state = prev;
        } else {
            return out;
        }
    }

    const int32_t raw_cost = mismatch_cnt + (2 * gap_open_cnt) + gap_extend_cnt;
    out.ok = true;
    out.mismatch_cnt = mismatch_cnt;
    out.gap_open_cnt = gap_open_cnt;
    out.gap_extend_cnt = gap_extend_cnt;
    out.aligned_cols = aligned_cols;
    out.target_consumed = target_consumed;
    out.norm_cost = static_cast<double>(raw_cost) / static_cast<double>(std::max(1, aligned_cols));
    return out;
}

double phi_value(const Placement& p, ReferenceSide side, InsertionTheta theta) {
    if (theta == InsertionTheta::kFwd) {
        if (side == ReferenceSide::kRefLeft) {
            return static_cast<double>(p.te_start);
        }
        if (side == ReferenceSide::kRefRight) {
            return static_cast<double>(p.te_end);
        }
        return static_cast<double>(p.te_start);
    }
    if (theta == InsertionTheta::kRev) {
        if (side == ReferenceSide::kRefLeft) {
            return static_cast<double>(p.te_end);
        }
        if (side == ReferenceSide::kRefRight) {
            return static_cast<double>(p.te_start);
        }
        return static_cast<double>(p.te_start);
    }
    return static_cast<double>(p.te_start);
}

}  // namespace

struct DeterministicAnchorLockedModule::TemplateDb {
    int32_t k = 13;
    std::unordered_map<std::string, TeRecord> records_by_name;

    bool build_from_fasta(const std::string& fasta_path, int32_t seed_k) {
        k = std::max(5, seed_k);

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

            TeRecord rec;
            rec.name = take_header_token(header);
            rec.sequence = upper_acgt(seq);
            if (static_cast<int32_t>(rec.sequence.size()) >= k) {
                for (int32_t i = 0; i + k <= static_cast<int32_t>(rec.sequence.size()); ++i) {
                    uint64_t key;
                    if (!build_kmer(rec.sequence, i, k, key)) {
                        continue;
                    }
                    rec.kmer_to_positions[key].push_back(i);
                }
            }

            if (!rec.name.empty() && !rec.sequence.empty()) {
                records_by_name.emplace(rec.name, std::move(rec));
            }

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

        return !records_by_name.empty();
    }
};

DeterministicAnchorLockedModule::DeterministicAnchorLockedModule(PipelineConfig config)
    : config_(std::move(config)) {
    if (!config_.te_consensus_enable || config_.te_fasta_path.empty()) {
        return;
    }

    auto db = std::make_shared<TemplateDb>();
    if (!db->build_from_fasta(config_.te_fasta_path, config_.te_consensus_seed_k)) {
        return;
    }
    template_db_ = std::move(db);
}

bool DeterministicAnchorLockedModule::is_enabled() const {
    return static_cast<bool>(template_db_);
}

AnchorLockedReport DeterministicAnchorLockedModule::resolve(
    const ComponentCall& component,
    const std::vector<InsertionFragment>& fragments,
    const std::vector<FragmentTEHit>& hits,
    const ClusterTECall& te_call) const {
    AnchorLockedReport report;
    report.enabled = is_enabled();
    report.total_fragments = static_cast<int32_t>(fragments.size());
    if (!is_enabled() || !template_db_) {
        return report;
    }

    std::unordered_map<std::string, const FragmentTEHit*> hit_by_fragment;
    hit_by_fragment.reserve(hits.size());
    for (const auto& h : hits) {
        hit_by_fragment[h.fragment_id] = &h;
    }

    std::string target_te_name;
    if (te_call.passed && !te_call.te_name.empty()) {
        target_te_name = te_call.te_name;
    } else {
        std::unordered_map<std::string, int32_t> votes;
        for (const auto& h : hits) {
            if (!h.te_name.empty()) {
                votes[h.te_name] += 1;
            }
        }
        int32_t best_votes = 0;
        for (const auto& kv : votes) {
            if (kv.second > best_votes ||
                (kv.second == best_votes && (target_te_name.empty() || kv.first < target_te_name))) {
                target_te_name = kv.first;
                best_votes = kv.second;
            }
        }
    }
    report.te_name = target_te_name;

    if (target_te_name.empty()) {
        return report;
    }

    const auto rec_it = template_db_->records_by_name.find(target_te_name);
    if (rec_it == template_db_->records_by_name.end()) {
        return report;
    }
    const TeRecord& te = rec_it->second;
    if (te.sequence.empty()) {
        return report;
    }

    auto enumerate_starts = [&](const std::string& query) -> SeedStarts {
        SeedStarts out;
        const int32_t k = template_db_->k;
        if (k <= 0 || static_cast<int32_t>(query.size()) < k || te.kmer_to_positions.empty()) {
            out.starts.push_back(0);
            return out;
        }

        std::unordered_map<int32_t, int32_t> counts;
        for (int32_t q = 0; q + k <= static_cast<int32_t>(query.size()); ++q) {
            uint64_t key;
            if (!build_kmer(query, q, k, key)) {
                continue;
            }
            const auto it = te.kmer_to_positions.find(key);
            if (it == te.kmer_to_positions.end()) {
                continue;
            }
            const auto& positions = it->second;
            const size_t max_positions = std::min<size_t>(positions.size(), 128);
            for (size_t pi = 0; pi < max_positions; ++pi) {
                const int32_t tpos = positions[pi];
                const int32_t start = tpos - q;
                if (start < 0 || start >= static_cast<int32_t>(te.sequence.size())) {
                    continue;
                }
                counts[start] += 1;
            }
        }

        if (counts.empty()) {
            out.starts.push_back(0);
            return out;
        }

        std::vector<std::pair<int32_t, int32_t>> rows;
        rows.reserve(counts.size());
        for (const auto& kv : counts) {
            rows.push_back(kv);
        }
        std::sort(rows.begin(), rows.end(), [](const auto& a, const auto& b) {
            if (a.second != b.second) {
                return a.second > b.second;
            }
            return a.first < b.first;
        });

        const int32_t keep = std::max(1, config_.te_consensus_max_start_candidates);
        for (int32_t i = 0; i < keep && i < static_cast<int32_t>(rows.size()); ++i) {
            out.starts.push_back(rows[static_cast<size_t>(i)].first);
        }

        if (out.starts.empty()) {
            out.starts.push_back(0);
        }
        std::sort(out.starts.begin(), out.starts.end());
        out.starts.erase(std::unique(out.starts.begin(), out.starts.end()), out.starts.end());
        return out;
    };

    std::vector<FragmentPlacementState> states;
    states.reserve(fragments.size());

    for (const auto& frag : fragments) {
        const auto hit_it = hit_by_fragment.find(frag.fragment_id);
        if (hit_it == hit_by_fragment.end() || !hit_it->second) {
            continue;
        }
        const FragmentTEHit* hit = hit_it->second;
        if (hit->te_name != target_te_name) {
            continue;
        }
        if (frag.sequence.empty()) {
            continue;
        }

        FragmentPlacementState state;
        state.fragment = &frag;
        state.hit = hit;
        state.split_sa = (frag.source == InsertionFragmentSource::kSplitSa);

        const std::string query = upper_acgt(frag.sequence);
        const SeedStarts starts = enumerate_starts(query);
        state.placements.reserve(starts.starts.size());

        for (int32_t start : starts.starts) {
            if (start < 0 || start >= static_cast<int32_t>(te.sequence.size())) {
                continue;
            }
            const std::string target_suffix = te.sequence.substr(static_cast<size_t>(start));
            if (target_suffix.empty()) {
                continue;
            }

            const AffineAlignResult aln = semiglobal_affine_align_suffix(query, target_suffix);
            if (!aln.ok) {
                continue;
            }
            Placement p;
            p.te_start = start;
            p.te_end = start + std::max(0, aln.target_consumed);
            p.norm_cost = aln.norm_cost;
            state.placements.push_back(p);
        }

        if (state.placements.empty()) {
            continue;
        }

        std::sort(state.placements.begin(), state.placements.end(), [](const Placement& a, const Placement& b) {
            if (a.norm_cost != b.norm_cost) {
                return a.norm_cost < b.norm_cost;
            }
            return a.te_start < b.te_start;
        });

        state.best_idx = 0;
        state.chosen_idx = 0;
        const double best_cost = state.placements.front().norm_cost;
        state.delta_tpl = (state.placements.size() >= 2)
            ? (state.placements[1].norm_cost - state.placements[0].norm_cost)
            : 1.0;

        int32_t span_min = state.placements.front().te_start;
        int32_t span_max = state.placements.front().te_start;
        for (const auto& p : state.placements) {
            if (p.norm_cost <= best_cost + 0.05) {
                span_min = std::min(span_min, p.te_start);
                span_max = std::max(span_max, p.te_start);
            }
        }
        state.start_span_best = static_cast<double>(span_max - span_min);

        const double span_s_max = config_.te_consensus_start_span_bias +
            (config_.te_consensus_start_span_sqrt_scale *
             std::sqrt(static_cast<double>(std::max(1, frag.length))));

        bool split_side_ok = true;
        if (state.split_sa && frag.ref_junc_pos >= 0) {
            split_side_ok =
                std::abs(frag.ref_junc_pos - component.anchor_pos) <= config_.te_consensus_w_side_max;
        }

        state.core_candidate =
            frag.ref_side != ReferenceSide::kUnknown &&
            frag.anchor_len >= config_.te_consensus_anchor_min &&
            state.start_span_best <= span_s_max &&
            state.delta_tpl >= config_.te_consensus_delta_tpl_min &&
            split_side_ok;

        states.push_back(std::move(state));
    }

    report.fragments_with_placements = static_cast<int32_t>(states.size());
    if (states.empty()) {
        return report;
    }

    auto weight_for = [&](const FragmentPlacementState& s, size_t placement_idx) {
        const Placement& p = s.placements[placement_idx];
        const double t = std::max(1e-6, config_.te_consensus_temp_t);
        return static_cast<double>(std::max(1, s.fragment->anchor_len)) * std::exp(-p.norm_cost / t);
    };

    std::vector<std::pair<double, double>> fwd_phi_w;
    std::vector<std::pair<double, double>> rev_phi_w;
    std::vector<double> fwd_values;
    std::vector<double> rev_values;
    for (const auto& s : states) {
        if (!s.core_candidate) {
            continue;
        }
        const Placement& p = s.placements[s.best_idx];
        const double w = weight_for(s, s.best_idx);
        const double phi_f = phi_value(p, s.fragment->ref_side, InsertionTheta::kFwd);
        const double phi_r = phi_value(p, s.fragment->ref_side, InsertionTheta::kRev);
        fwd_phi_w.push_back({phi_f, w});
        rev_phi_w.push_back({phi_r, w});
        fwd_values.push_back(phi_f);
        rev_values.push_back(phi_r);
        report.core_candidate_count += 1;
        if (s.split_sa) {
            report.split_sa_core_count += 1;
        }
        if (s.fragment->ref_junc_pos >= 0) {
            if (report.ref_junc_pos_min < 0) {
                report.ref_junc_pos_min = s.fragment->ref_junc_pos;
                report.ref_junc_pos_max = s.fragment->ref_junc_pos;
            } else {
                report.ref_junc_pos_min = std::min(report.ref_junc_pos_min, s.fragment->ref_junc_pos);
                report.ref_junc_pos_max = std::max(report.ref_junc_pos_max, s.fragment->ref_junc_pos);
            }
        }
    }

    if (report.core_candidate_count <= 0) {
        return report;
    }
    report.split_sa_core_frac = static_cast<double>(report.split_sa_core_count) /
        static_cast<double>(report.core_candidate_count);

    double sum_w = 0.0;
    for (const auto& row : fwd_phi_w) {
        sum_w += row.second;
    }
    report.sum_w_fwd = sum_w;
    report.sum_w_rev = sum_w;

    const double c_fwd = weighted_median(fwd_phi_w);
    const double c_rev = weighted_median(rev_phi_w);
    const double mad_fwd = mad_scaled(fwd_values, c_fwd);
    const double mad_rev = mad_scaled(rev_values, c_rev);
    report.mad_fwd = mad_fwd;
    report.mad_rev = mad_rev;

    const double penalty = config_.te_consensus_beta * (1.0 / std::max(1e-9, sum_w));
    const double score_fwd = mad_fwd + penalty;
    const double score_rev = mad_rev + penalty;

    InsertionTheta theta_for_phi = InsertionTheta::kFwd;
    if (score_rev < score_fwd) {
        theta_for_phi = InsertionTheta::kRev;
    } else if (score_rev == score_fwd && mad_rev < mad_fwd) {
        theta_for_phi = InsertionTheta::kRev;
    }

    report.theta0 = theta_for_phi;
    report.fail_theta_uncertain =
        std::abs(mad_fwd - mad_rev) < config_.te_consensus_eps_theta &&
        relative_diff(report.sum_w_fwd, report.sum_w_rev) < config_.te_consensus_rel_sumw_eps;
    if (report.fail_theta_uncertain) {
        report.theta0 = InsertionTheta::kUnknown;
    }

    std::vector<std::pair<double, double>> c0_values;
    c0_values.reserve(states.size());
    for (const auto& s : states) {
        if (!s.core_candidate) {
            continue;
        }
        const Placement& p = s.placements[s.best_idx];
        c0_values.push_back({phi_value(p, s.fragment->ref_side, theta_for_phi), weight_for(s, s.best_idx)});
    }
    const double c0 = weighted_median(c0_values);
    report.center_c0 = c0;

    for (auto& s : states) {
        if (s.placements.empty()) {
            continue;
        }
        size_t best_idx = 0;
        double best_obj = std::numeric_limits<double>::infinity();
        for (size_t pi = 0; pi < s.placements.size(); ++pi) {
            const Placement& p = s.placements[pi];
            const double phi = phi_value(p, s.fragment->ref_side, theta_for_phi);
            const double span_start = std::abs(static_cast<double>(p.te_start - s.placements[s.best_idx].te_start));
            const double obj =
                p.norm_cost +
                (config_.te_consensus_lambda * std::abs(phi - c0)) +
                (config_.te_consensus_mu * span_start);
            if (obj < best_obj ||
                (obj == best_obj && p.norm_cost < s.placements[best_idx].norm_cost) ||
                (obj == best_obj && p.norm_cost == s.placements[best_idx].norm_cost &&
                 p.te_start < s.placements[best_idx].te_start)) {
                best_obj = obj;
                best_idx = pi;
            }
        }
        s.chosen_idx = best_idx;
    }

    std::vector<std::pair<double, double>> c1_values;
    std::vector<double> c1_phi_values;
    c1_values.reserve(states.size());
    c1_phi_values.reserve(states.size());
    for (const auto& s : states) {
        if (!s.core_candidate) {
            continue;
        }
        const Placement& p = s.placements[s.chosen_idx];
        const double phi = phi_value(p, s.fragment->ref_side, theta_for_phi);
        const double w = weight_for(s, s.chosen_idx);
        c1_values.push_back({phi, w});
        c1_phi_values.push_back(phi);
    }

    const double c1 = weighted_median(c1_values);
    report.center_c1 = c1;

    std::vector<double> core_set_phi;
    for (const auto& s : states) {
        if (!s.core_candidate) {
            continue;
        }
        const Placement& p = s.placements[s.chosen_idx];
        const double phi = phi_value(p, s.fragment->ref_side, theta_for_phi);
        if (std::abs(phi - c1) <= config_.te_consensus_w_core_gate) {
            core_set_phi.push_back(phi);
            report.core_set_count += 1;
        }
    }

    if (core_set_phi.empty()) {
        return report;
    }

    const double core = median(core_set_phi);
    const double span = mad_scaled(core_set_phi, core);
    const double w = std::max(config_.te_consensus_w_min, config_.te_consensus_gamma * span);

    report.te_breakpoint_core = safe_floor_to_i32(core);
    report.core_span = span;
    report.te_breakpoint_window_start = safe_floor_to_i32(core - w);
    report.te_breakpoint_window_end = safe_ceil_to_i32(core + w);
    report.has_result = report.core_set_count >= config_.te_consensus_min_core_set;
    return report;
}

}  // namespace placer
