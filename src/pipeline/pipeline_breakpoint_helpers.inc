constexpr int32_t kBreakpointClusterGapBp = 75;
constexpr size_t kBreakpointHypothesisMaxClustersPerSide = 6;
constexpr size_t kBreakpointHypothesisMaxSingleClusters = 2;
constexpr int32_t kBreakpointPairCompatibilityBp = 250;

struct BreakpointPositionCluster {
    int32_t center = -1;
    int32_t start = -1;
    int32_t end = -1;
    int32_t support = 0;
};

std::vector<BreakpointPositionCluster> collect_breakpoint_position_clusters(
    std::vector<int32_t> positions) {
    std::vector<BreakpointPositionCluster> clusters;
    if (positions.empty()) {
        return clusters;
    }
    std::sort(positions.begin(), positions.end());

    size_t cluster_begin = 0;
    for (size_t i = 1; i <= positions.size(); ++i) {
        const bool split_cluster =
            i == positions.size() ||
            (positions[i] - positions[i - 1]) > kBreakpointClusterGapBp;
        if (!split_cluster) {
            continue;
        }

        BreakpointPositionCluster cluster;
        cluster.start = positions[cluster_begin];
        cluster.end = positions[i - 1];
        cluster.support = static_cast<int32_t>(i - cluster_begin);
        cluster.center = median_i32(std::vector<int32_t>(
            positions.begin() + static_cast<std::ptrdiff_t>(cluster_begin),
            positions.begin() + static_cast<std::ptrdiff_t>(i)));
        clusters.push_back(cluster);
        cluster_begin = i;
    }
    return clusters;
}

int32_t breakpoint_hypothesis_support_weight(int32_t priority) {
    switch (priority) {
        case 0: return 8;  // fragment split
        case 1: return 8;  // fragment indel
        case 2: return 7;  // raw split
        case 3: return 1;  // fragment clip pair
        case 4: return 6;  // raw indel
        case 5: return 1;  // raw clip pair
        default: return 1;
    }
}

struct LocalBreakpointHypothesis {
    bool valid = false;
    int32_t left = -1;
    int32_t right = -1;
    int32_t center = -1;
    int32_t support = 0;
    int32_t priority = std::numeric_limits<int32_t>::max();
};

double mad_from_positions(const std::vector<int32_t>& positions) {
    if (positions.empty()) {
        return 0.0;
    }
    const int32_t center = median_i32(positions);
    std::vector<int32_t> abs_dev;
    abs_dev.reserve(positions.size());
    for (int32_t pos : positions) {
        abs_dev.push_back(std::abs(pos - center));
    }
    return static_cast<double>(median_i32(std::move(abs_dev)));
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
    int32_t at_count = 0;
    int32_t total = 0;
    for (char c : seq) {
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
            continue;
        }
        ++total;
        if (c == 'A' || c == 'T') {
            ++at_count;
        }
    }
    if (total <= 0) {
        return 0.0;
    }
    return static_cast<double>(at_count) / static_cast<double>(total);
}

double shannon_entropy_acgt(const std::string& seq) {
    std::array<int32_t, 4> counts = {0, 0, 0, 0};
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
    for (int32_t count : counts) {
        if (count <= 0) {
            continue;
        }
        const double p = static_cast<double>(count) / static_cast<double>(total);
        entropy -= p * std::log2(p);
    }
    return entropy;
}

uint8_t char_to_2bit(char c);
bool build_kmer(const std::string& s, int32_t start, int32_t k, uint64_t& out);
double kmer_uniqueness_ratio(const std::string& seq, int32_t k);

bool is_low_complexity_softclip(
    const InsertionFragment& fragment,
    const PipelineConfig& config) {
    if (!is_softclip_source(fragment.source)) {
        return false;
    }
    if (fragment.sequence.empty()) {
        return true;
    }
    const double at_frac = at_fraction(fragment.sequence);
    const int32_t homopolymer = max_homopolymer_run(fragment.sequence);
    const double entropy = shannon_entropy_acgt(fragment.sequence);
    const double kmer_uniqueness = kmer_uniqueness_ratio(fragment.sequence, 5);
    return at_frac >= config.te_softclip_low_complexity_at_frac_min ||
           homopolymer >= config.te_softclip_low_complexity_homopolymer_min ||
           entropy < std::max(0.0, config.te_softclip_entropy_min) ||
           kmer_uniqueness < std::clamp(config.te_softclip_kmer_uniqueness_min, 0.0, 1.0);
}

bool is_low_quality_softclip_anchor(
    const InsertionFragment& fragment,
    const PipelineConfig& config) {
    if (!is_softclip_source(fragment.source)) {
        return false;
    }
    if (fragment.anchor_len > 0 &&
        fragment.anchor_len < std::max(1, config.te_softclip_min_anchor_len)) {
        return true;
    }
    if (fragment.anchor_len > 0 && fragment.nm >= 0) {
        const double nm_per_bp =
            static_cast<double>(fragment.nm) / static_cast<double>(std::max(1, fragment.anchor_len));
        if (nm_per_bp > std::max(0.0, config.te_softclip_max_nm_per_bp)) {
            return true;
        }
    }
    return false;
}

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

std::unordered_set<uint64_t> kmer_set(const std::string& seq, int32_t k) {
    std::unordered_set<uint64_t> out;
    if (k <= 0 || static_cast<int32_t>(seq.size()) < k) {
        return out;
    }
    out.reserve(static_cast<size_t>(seq.size()));
    for (int32_t i = 0; i + k <= static_cast<int32_t>(seq.size()); ++i) {
        uint64_t key = 0;
        if (!build_kmer(seq, i, k, key)) {
            continue;
        }
        out.insert(key);
    }
    return out;
}

double kmer_uniqueness_ratio(const std::string& seq, int32_t k) {
    if (k <= 0 || static_cast<int32_t>(seq.size()) < k) {
        return 0.0;
    }
    int32_t total = 0;
    std::unordered_set<uint64_t> uniq;
    uniq.reserve(static_cast<size_t>(seq.size()));
    for (int32_t i = 0; i + k <= static_cast<int32_t>(seq.size()); ++i) {
        uint64_t key = 0;
        if (!build_kmer(seq, i, k, key)) {
            continue;
        }
        uniq.insert(key);
        total += 1;
    }
    if (total <= 0) {
        return 0.0;
    }
    return static_cast<double>(uniq.size()) / static_cast<double>(total);
}

std::string upper_acgt(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        out.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
    }
    return out;
}

char complement_acgt(char c) {
    switch (static_cast<char>(std::toupper(static_cast<unsigned char>(c)))) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return 'N';
    }
}

std::string reverse_complement_acgt(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (auto it = s.rbegin(); it != s.rend(); ++it) {
        out.push_back(complement_acgt(*it));
    }
    return out;
}

struct EditDistanceWorkspace {
    std::vector<int32_t> prev;
    std::vector<int32_t> curr;

    void ensure_columns(int32_t m) {
        const size_t columns = static_cast<size_t>(m + 1);
        if (prev.size() < columns) {
            prev.resize(columns, 0);
        }
        if (curr.size() < columns) {
            curr.resize(columns, 0);
        }
    }
};

int32_t max_edits_for_identity_threshold(
    int32_t lhs_len,
    int32_t rhs_len,
    double min_identity) {
    const int32_t denom = std::max(lhs_len, rhs_len);
    if (denom <= 0) {
        return 0;
    }
    const double clamped = std::clamp(min_identity, 0.0, 1.0);
    return std::max(
        0,
        static_cast<int32_t>(
            std::floor(((1.0 - clamped) * static_cast<double>(denom)) + 1e-9)));
}

bool edit_identity_if_at_least(
    std::string_view lhs,
    std::string_view rhs,
    int32_t max_edits,
    EditDistanceWorkspace& workspace,
    double& identity_out) {
    identity_out = 0.0;
    const int32_t n = static_cast<int32_t>(lhs.size());
    const int32_t m = static_cast<int32_t>(rhs.size());
    if (n <= 0 || m <= 0 || max_edits < 0) {
        return false;
    }
    if (std::abs(n - m) > max_edits) {
        return false;
    }

    workspace.ensure_columns(m);
    auto& prev = workspace.prev;
    auto& curr = workspace.curr;
    const size_t columns = static_cast<size_t>(m + 1);
    const int32_t inf = max_edits + 1;
    std::fill(prev.begin(), prev.begin() + columns, inf);
    std::fill(curr.begin(), curr.begin() + columns, inf);
    for (int32_t j = 0; j <= std::min(m, max_edits); ++j) {
        prev[static_cast<size_t>(j)] = j;
    }

    for (int32_t i = 1; i <= n; ++i) {
        std::fill(curr.begin(), curr.begin() + columns, inf);
        if (i <= max_edits) {
            curr[0] = i;
        }

        const int32_t j_lo = std::max(1, i - max_edits);
        const int32_t j_hi = std::min(m, i + max_edits);
        if (j_lo > j_hi) {
            return false;
        }

        for (int32_t j = j_lo; j <= j_hi; ++j) {
            const int32_t sub_cost =
                prev[static_cast<size_t>(j - 1)] +
                ((lhs[static_cast<size_t>(i - 1)] == rhs[static_cast<size_t>(j - 1)]) ? 0 : 1);
            const int32_t del_cost = prev[static_cast<size_t>(j)] + 1;
            const int32_t ins_cost = curr[static_cast<size_t>(j - 1)] + 1;
            curr[static_cast<size_t>(j)] = std::min({sub_cost, del_cost, ins_cost});
        }
        std::swap(prev, curr);
    }

    const int32_t dist = prev[static_cast<size_t>(m)];
    if (dist > max_edits) {
        return false;
    }

    const double denom = static_cast<double>(std::max(n, m));
    identity_out = std::clamp(1.0 - (static_cast<double>(dist) / denom), 0.0, 1.0);
    return true;
}

bool edit_identity_if_at_least_fast_impl(
    std::string_view lhs,
    std::string_view rhs,
    int32_t max_edits,
    EditDistanceWorkspace& workspace,
    double& identity_out) {
    identity_out = 0.0;
    const int32_t n = static_cast<int32_t>(lhs.size());
    const int32_t m = static_cast<int32_t>(rhs.size());
    if (n <= 0 || m <= 0 || max_edits < 0) {
        return false;
    }
    if (std::abs(n - m) > max_edits) {
        return false;
    }

    const int32_t neg = std::numeric_limits<int32_t>::min() / 4;
    const int32_t offset = max_edits + 1;
    const size_t width = static_cast<size_t>((2 * max_edits) + 3);
    if (workspace.prev.size() < width) {
        workspace.prev.resize(width, neg);
    }
    if (workspace.curr.size() < width) {
        workspace.curr.resize(width, neg);
    }
    auto& prev = workspace.prev;
    auto& curr = workspace.curr;
    const auto index_for_diag = [offset](int32_t diag) {
        return static_cast<size_t>(diag + offset);
    };
    const auto extend_match = [&](int32_t x, int32_t diag) {
        int32_t y = x - diag;
        if (x < 0 || y < 0 || x > n || y > m) {
            return neg;
        }
        while (x < n && y < m &&
               lhs[static_cast<size_t>(x)] == rhs[static_cast<size_t>(y)]) {
            ++x;
            ++y;
        }
        return x;
    };
    const auto reached_end = [m](int32_t x, int32_t diag) {
        return x >= 0 && (x - diag) == m;
    };
    const auto maybe_finish = [&](int32_t dist, int32_t x, int32_t diag) {
        if (x == n && reached_end(x, diag)) {
            const double denom = static_cast<double>(std::max(n, m));
            identity_out = std::clamp(
                1.0 - (static_cast<double>(dist) / denom),
                0.0,
                1.0);
            return true;
        }
        return false;
    };

    std::fill(prev.begin(), prev.begin() + width, neg);
    std::fill(curr.begin(), curr.begin() + width, neg);
    prev[index_for_diag(0)] = extend_match(0, 0);
    if (maybe_finish(0, prev[index_for_diag(0)], 0)) {
        return true;
    }

    for (int32_t edits = 1; edits <= max_edits; ++edits) {
        std::fill(curr.begin(), curr.begin() + width, neg);
        const int32_t prev_lo = -(edits - 1);
        const int32_t prev_hi = edits - 1;
        for (int32_t diag = -edits; diag <= edits; ++diag) {
            int32_t best_x = neg;

            if (diag - 1 >= prev_lo && diag - 1 <= prev_hi) {
                const int32_t x0 = prev[index_for_diag(diag - 1)];
                if (x0 >= 0 && x0 < n) {
                    best_x = std::max(best_x, extend_match(x0 + 1, diag));
                }
            }
            if (diag + 1 >= prev_lo && diag + 1 <= prev_hi) {
                const int32_t x0 = prev[index_for_diag(diag + 1)];
                const int32_t y0 = x0 - (diag + 1);
                if (x0 >= 0 && y0 >= 0 && y0 < m) {
                    best_x = std::max(best_x, extend_match(x0, diag));
                }
            }
            if (diag >= prev_lo && diag <= prev_hi) {
                const int32_t x0 = prev[index_for_diag(diag)];
                const int32_t y0 = x0 - diag;
                if (x0 >= 0 && x0 < n && y0 >= 0 && y0 < m) {
                    best_x = std::max(best_x, extend_match(x0 + 1, diag));
                }
            }

            curr[index_for_diag(diag)] = best_x;
            if (maybe_finish(edits, best_x, diag)) {
                return true;
            }
        }
        std::swap(prev, curr);
    }

    return false;
}

double fixed_window_edit_identity(
    const std::string& query,
    const std::string& target) {
    if (query.empty() || target.empty()) {
        return 0.0;
    }
    const int32_t max_edits = static_cast<int32_t>(
        std::max(query.size(), target.size()));
    EditDistanceWorkspace workspace;
    double identity = 0.0;
    if (!edit_identity_if_at_least(
            query,
            target,
            max_edits,
            workspace,
            identity)) {
        return 0.0;
    }
    return identity;
}

uint8_t nt_to_abpoa_code(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return 4;
    }
}

char abpoa_code_to_nt(uint8_t c) {
    switch (c) {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        default: return 'N';
    }
}

std::string abpoa_consensus(const std::vector<std::string>& seqs) {
    if (seqs.empty()) {
        return "";
    }

    abpoa_t* ab = abpoa_init();
    abpoa_para_t* abpt = abpoa_init_para();
    if (!ab || !abpt) {
        if (ab) {
            abpoa_free(ab);
        }
        if (abpt) {
            abpoa_free_para(abpt);
        }
        return "";
    }

    abpt->out_cons = 1;
    abpt->out_msa = 0;
    abpt->out_gfa = 0;
    abpt->max_n_cons = 1;
    abpoa_post_set_para(abpt);

    const int32_t n_seq = static_cast<int32_t>(seqs.size());
    std::vector<int> seq_lens(static_cast<size_t>(n_seq), 0);
    std::vector<std::vector<uint8_t>> packed(static_cast<size_t>(n_seq));
    std::vector<uint8_t*> seq_ptrs(static_cast<size_t>(n_seq), nullptr);

    for (int32_t i = 0; i < n_seq; ++i) {
        const std::string seq = upper_acgt(seqs[static_cast<size_t>(i)]);
        seq_lens[static_cast<size_t>(i)] = static_cast<int>(seq.size());
        packed[static_cast<size_t>(i)].resize(seq.size(), 4);
        for (size_t j = 0; j < seq.size(); ++j) {
            packed[static_cast<size_t>(i)][j] = nt_to_abpoa_code(seq[j]);
        }
        seq_ptrs[static_cast<size_t>(i)] = packed[static_cast<size_t>(i)].data();
    }

    const int ret = abpoa_msa(
        ab,
        abpt,
        n_seq,
        nullptr,
        seq_lens.data(),
        seq_ptrs.data(),
        nullptr,
        nullptr);

    std::string consensus;
    if (ret == 0 && ab->abc && ab->abc->n_cons > 0) {
        const int cons_len = ab->abc->cons_len[0];
        if (cons_len > 0 && ab->abc->cons_base && ab->abc->cons_base[0]) {
            consensus.reserve(static_cast<size_t>(cons_len));
            for (int i = 0; i < cons_len; ++i) {
                consensus.push_back(abpoa_code_to_nt(ab->abc->cons_base[0][i]));
            }
        }
    }

    abpoa_free(ab);
    abpoa_free_para(abpt);
    return consensus;
}

bool has_prefix(const std::string& value, const char* prefix) {
    return value.rfind(prefix, 0) == 0;
}

bool has_qc_token(const std::string& qc, const char* token) {
    if (!token || *token == '\0') {
        return false;
    }
    size_t start = 0;
    while (start <= qc.size()) {
        const size_t end = qc.find('|', start);
        const size_t len = (end == std::string::npos) ? (qc.size() - start) : (end - start);
        if (len == std::strlen(token) && qc.compare(start, len, token) == 0) {
            return true;
        }
        if (end == std::string::npos) {
            break;
        }
        start = end + 1;
    }
    return false;
}

void append_qc_token(std::string& qc, const char* token) {
    if (!token || *token == '\0' || has_qc_token(qc, token)) {
        return;
    }
    if (!qc.empty()) {
        qc += "|";
    }
    qc += token;
}

void update_best_bp_by_read(
    std::unordered_map<size_t, int32_t>& bp_by_read,
    size_t read_index,
    int32_t pos,
    int32_t anchor_pos) {
    if (pos < 0) {
        return;
    }
    const auto it = bp_by_read.find(read_index);
    if (it == bp_by_read.end()) {
        bp_by_read.emplace(read_index, pos);
        return;
    }
    const int32_t cur_dist = std::abs(it->second - anchor_pos);
    const int32_t new_dist = std::abs(pos - anchor_pos);
    if (new_dist < cur_dist || (new_dist == cur_dist && pos < it->second)) {
        it->second = pos;
    }
}

std::pair<int32_t, int32_t> infer_component_breakpoint_bounds(
    const ComponentCall& component) {
    int32_t bp_min = std::numeric_limits<int32_t>::max();
    int32_t bp_max = std::numeric_limits<int32_t>::min();

    for (const auto& bp : component.breakpoint_candidates) {
        if (bp.pos < 0) {
            continue;
        }
        bp_min = std::min(bp_min, bp.pos);
        bp_max = std::max(bp_max, bp.pos);
    }

    if (bp_min <= bp_max) {
        return {bp_min, bp_max};
    }

    const int32_t fallback = std::max(0, component.anchor_pos);
    return {fallback, fallback};
}

std::vector<LocalBreakpointHypothesis> collect_breakpoint_hypotheses(
    const ComponentCall& component,
    const std::vector<const bam1_t*>& local_records,
    const std::vector<InsertionFragment>& fragments,
    int32_t seed_left,
    int32_t seed_right,
    size_t top_k) {
    constexpr int32_t kEventBreakpointSearchSlackBp = 200;
    constexpr int32_t kPreciseFragmentSearchSlackBp = 200;
    constexpr int32_t kPreciseFragmentComponentLinkBp = kLocalEventFetchSlackBp;

    int32_t search_start = std::max(
        0,
        std::min(seed_left, seed_right) - kEventBreakpointSearchSlackBp);
    int32_t search_end = std::max(
        search_start + 1,
        std::max(seed_left, seed_right) + kEventBreakpointSearchSlackBp);
    const auto expand_search_to_precise_pos = [&](int32_t pos) {
        if (pos < 0) {
            return;
        }
        search_start = std::max(
            0,
            std::min(search_start, pos - kPreciseFragmentSearchSlackBp));
        search_end = std::max(
            search_end,
            pos + kPreciseFragmentSearchSlackBp);
    };
    const auto linked_to_component_seed = [&](int32_t pos) {
        if (pos < 0) {
            return false;
        }
        const int32_t seed_min = std::min(seed_left, seed_right);
        const int32_t seed_max = std::max(seed_left, seed_right);
        const int32_t distance = (pos < seed_min)
            ? (seed_min - pos)
            : ((pos > seed_max) ? (pos - seed_max) : 0);
        return distance <= kPreciseFragmentComponentLinkBp;
    };
    for (const auto& candidate : component.breakpoint_candidates) {
        if ((candidate.class_mask & kCandidateLongInsertion) != 0 ||
            (candidate.class_mask & kCandidateSplitSaSupplementary) != 0 ||
            candidate.ins_len >= kLongInsertionSignalMin) {
            expand_search_to_precise_pos(candidate.pos);
        }
    }
    for (const auto& fragment : fragments) {
        if ((fragment.source == InsertionFragmentSource::kCigarInsertion ||
             fragment.source == InsertionFragmentSource::kSplitSa) &&
            linked_to_component_seed(fragment.ref_junc_pos)) {
            expand_search_to_precise_pos(fragment.ref_junc_pos);
        }
    }

    std::vector<int32_t> fragment_split_left_positions;
    std::vector<int32_t> fragment_split_right_positions;
    std::vector<int32_t> fragment_indel_positions;
    std::vector<int32_t> fragment_clip_left_positions;
    std::vector<int32_t> fragment_clip_right_positions;
    std::vector<int32_t> raw_split_left_positions;
    std::vector<int32_t> raw_split_right_positions;
    std::vector<int32_t> raw_indel_positions;
    std::vector<int32_t> raw_clip_left_positions;
    std::vector<int32_t> raw_clip_right_positions;

    auto add_position = [&](std::vector<int32_t>& out, int32_t pos) {
        if (pos >= search_start && pos <= search_end) {
            out.push_back(pos);
        }
    };

    for (const bam1_t* record : local_records) {
        if (!record) {
            continue;
        }
        ReadView read(record);
        if (read.tid() != component.tid) {
            continue;
        }

        const LocalEventSignal signal = classify_local_event_signal(
            record,
            component.chrom,
            search_start,
            search_end);
        const bool has_raw_split_left = signal.split_left_pos >= 0;
        const bool has_raw_split_right = signal.split_right_pos >= 0;
        if (has_raw_split_left && has_raw_split_right) {
            const int32_t left_dist = std::abs(signal.split_left_pos - component.anchor_pos);
            const int32_t right_dist = std::abs(signal.split_right_pos - component.anchor_pos);
            if (left_dist < right_dist ||
                (left_dist == right_dist && signal.split_left_pos <= signal.split_right_pos)) {
                add_position(raw_split_left_positions, signal.split_left_pos);
            } else {
                add_position(raw_split_right_positions, signal.split_right_pos);
            }
        } else {
            add_position(raw_split_left_positions, signal.split_left_pos);
            add_position(raw_split_right_positions, signal.split_right_pos);
        }
        add_position(raw_clip_left_positions, signal.left_clip_pos);
        add_position(raw_clip_right_positions, signal.right_clip_pos);
        add_position(raw_indel_positions, signal.indel_pos);
    }

    for (const auto& fragment : fragments) {
        if (fragment.ref_junc_pos < search_start || fragment.ref_junc_pos > search_end) {
            continue;
        }
        switch (fragment.source) {
            case InsertionFragmentSource::kSplitSa:
                if (fragment.ref_side == ReferenceSide::kRefLeft) {
                    fragment_split_left_positions.push_back(fragment.ref_junc_pos);
                } else if (fragment.ref_side == ReferenceSide::kRefRight) {
                    fragment_split_right_positions.push_back(fragment.ref_junc_pos);
                }
                break;
            case InsertionFragmentSource::kCigarInsertion:
                fragment_indel_positions.push_back(fragment.ref_junc_pos);
                break;
            case InsertionFragmentSource::kClipRefLeft:
                fragment_clip_left_positions.push_back(fragment.ref_junc_pos);
                break;
            case InsertionFragmentSource::kClipRefRight:
                fragment_clip_right_positions.push_back(fragment.ref_junc_pos);
                break;
            default:
                break;
        }
    }

    const auto sort_clusters = [&](std::vector<BreakpointPositionCluster>& clusters) {
        std::sort(clusters.begin(), clusters.end(), [&](const auto& lhs, const auto& rhs) {
            if (lhs.support != rhs.support) {
                return lhs.support > rhs.support;
            }
            const int32_t lhs_dist = std::abs(lhs.center - component.anchor_pos);
            const int32_t rhs_dist = std::abs(rhs.center - component.anchor_pos);
            if (lhs_dist != rhs_dist) {
                return lhs_dist < rhs_dist;
            }
            return lhs.center < rhs.center;
        });
    };

    const auto make_single_breakpoint_hypotheses = [&](std::vector<int32_t> positions,
                                                       int32_t priority) {
        std::vector<LocalBreakpointHypothesis> out;
        auto clusters = collect_breakpoint_position_clusters(std::move(positions));
        sort_clusters(clusters);
        if (clusters.size() > kBreakpointHypothesisMaxSingleClusters) {
            clusters.resize(kBreakpointHypothesisMaxSingleClusters);
        }

        out.reserve(clusters.size());
        for (const auto& cluster : clusters) {
            LocalBreakpointHypothesis hypothesis;
            hypothesis.valid = cluster.center >= 0;
            hypothesis.left = cluster.center;
            hypothesis.right = cluster.center;
            hypothesis.center = cluster.center;
            hypothesis.support = cluster.support;
            hypothesis.priority = priority;
            out.push_back(hypothesis);
        }
        return out;
    };

    const auto make_paired_breakpoint_hypotheses = [&](std::vector<int32_t> left_positions,
                                                       std::vector<int32_t> right_positions,
                                                       int32_t priority) {
        std::vector<LocalBreakpointHypothesis> out;
        auto left_clusters = collect_breakpoint_position_clusters(std::move(left_positions));
        auto right_clusters = collect_breakpoint_position_clusters(std::move(right_positions));
        sort_clusters(left_clusters);
        sort_clusters(right_clusters);

        if (left_clusters.size() > kBreakpointHypothesisMaxClustersPerSide) {
            left_clusters.resize(kBreakpointHypothesisMaxClustersPerSide);
        }
        if (right_clusters.size() > kBreakpointHypothesisMaxClustersPerSide) {
            right_clusters.resize(kBreakpointHypothesisMaxClustersPerSide);
        }

        if (left_clusters.empty() && right_clusters.empty()) {
            return out;
        }

        if (left_clusters.empty()) {
            out.reserve(right_clusters.size());
            for (const auto& right_cluster : right_clusters) {
                LocalBreakpointHypothesis hypothesis;
                hypothesis.valid = right_cluster.center >= 0;
                hypothesis.left = right_cluster.center;
                hypothesis.right = right_cluster.center;
                hypothesis.center = right_cluster.center;
                hypothesis.support = right_cluster.support;
                hypothesis.priority = priority;
                out.push_back(hypothesis);
            }
            return out;
        }

        if (right_clusters.empty()) {
            out.reserve(left_clusters.size());
            for (const auto& left_cluster : left_clusters) {
                LocalBreakpointHypothesis hypothesis;
                hypothesis.valid = left_cluster.center >= 0;
                hypothesis.left = left_cluster.center;
                hypothesis.right = left_cluster.center;
                hypothesis.center = left_cluster.center;
                hypothesis.support = left_cluster.support;
                hypothesis.priority = priority;
                out.push_back(hypothesis);
            }
            return out;
        }

        out.reserve(left_clusters.size() * right_clusters.size());
        for (const auto& left_cluster : left_clusters) {
            for (const auto& right_cluster : right_clusters) {
                const int32_t left = std::min(left_cluster.center, right_cluster.center);
                const int32_t right = std::max(left_cluster.center, right_cluster.center);
                if ((right - left) > kBreakpointPairCompatibilityBp) {
                    continue;
                }

                LocalBreakpointHypothesis hypothesis;
                hypothesis.valid = true;
                hypothesis.left = left;
                hypothesis.right = right;
                hypothesis.center =
                    static_cast<int32_t>((static_cast<int64_t>(left) + right) / 2);
                hypothesis.support = left_cluster.support + right_cluster.support;
                hypothesis.priority = priority;
                out.push_back(hypothesis);
            }
        }
        return out;
    };

    std::vector<LocalBreakpointHypothesis> hypotheses;
    hypotheses.reserve(24);
    const auto append_unique = [&](const std::vector<LocalBreakpointHypothesis>& batch) {
        for (const auto& hypothesis : batch) {
            if (!hypothesis.valid) {
                continue;
            }
            bool duplicate = false;
            for (const auto& existing : hypotheses) {
                if (existing.left == hypothesis.left &&
                    existing.right == hypothesis.right) {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate) {
                hypotheses.push_back(hypothesis);
            }
        }
    };

    append_unique(make_paired_breakpoint_hypotheses(
        fragment_split_left_positions,
        fragment_split_right_positions,
        0));
    append_unique(make_single_breakpoint_hypotheses(
        fragment_indel_positions,
        1));
    append_unique(make_paired_breakpoint_hypotheses(
        raw_split_left_positions,
        raw_split_right_positions,
        2));
    append_unique(make_paired_breakpoint_hypotheses(
        fragment_clip_left_positions,
        fragment_clip_right_positions,
        3));
    append_unique(make_single_breakpoint_hypotheses(
        raw_indel_positions,
        4));
    append_unique(make_paired_breakpoint_hypotheses(
        raw_clip_left_positions,
        raw_clip_right_positions,
        5));

    if (hypotheses.empty()) {
        LocalBreakpointHypothesis fallback;
        fallback.valid = true;
        fallback.left = std::max(0, component.anchor_pos);
        fallback.right = fallback.left;
        fallback.center = fallback.left;
        fallback.support = 1;
        fallback.priority = std::numeric_limits<int32_t>::max();
        hypotheses.push_back(fallback);
    }

    std::sort(
        hypotheses.begin(),
        hypotheses.end(),
        [&](const LocalBreakpointHypothesis& lhs, const LocalBreakpointHypothesis& rhs) {
            const int32_t lhs_score =
                lhs.support * breakpoint_hypothesis_support_weight(lhs.priority);
            const int32_t rhs_score =
                rhs.support * breakpoint_hypothesis_support_weight(rhs.priority);
            if (lhs_score != rhs_score) {
                return lhs_score > rhs_score;
            }
            const int32_t lhs_dist = std::abs(lhs.center - component.anchor_pos);
            const int32_t rhs_dist = std::abs(rhs.center - component.anchor_pos);
            if (lhs_dist != rhs_dist) {
                return lhs_dist < rhs_dist;
            }
            if (lhs.support != rhs.support) {
                return lhs.support > rhs.support;
            }
            if (lhs.priority != rhs.priority) {
                return lhs.priority < rhs.priority;
            }
            if (lhs.left != rhs.left) {
                return lhs.left < rhs.left;
            }
            return lhs.right < rhs.right;
        });
    if (top_k > 0 && hypotheses.size() > top_k) {
        hypotheses.resize(top_k);
    }
    return hypotheses;
}
