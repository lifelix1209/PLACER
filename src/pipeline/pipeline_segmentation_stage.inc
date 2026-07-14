bool Pipeline::edit_identity_if_at_least_fast(
    std::string_view lhs,
    std::string_view rhs,
    int32_t max_edits,
    double& identity_out) {
    EditDistanceWorkspace workspace;
    return edit_identity_if_at_least_fast_impl(
        lhs,
        rhs,
        max_edits,
        workspace,
        identity_out);
}

std::string Pipeline::pre_segmentation_gate_reason(
    const EventReadEvidence& event_evidence,
    const EventConsensus& event_consensus) const {
    if (event_consensus.left_anchor_input_reads <= 0 ||
        event_consensus.right_anchor_input_reads <= 0) {
        return "PRESEG_NO_BILATERAL_ANCHOR";
    }
    const bool has_precise_or_full_context =
        event_consensus.full_context_input_reads > 0 ||
        (event_evidence.alt_split_reads + event_evidence.alt_indel_reads) > 0;
    const bool has_identifiable_bilateral_partial_context =
        has_identifiable_bilateral_partial_context_support(
            event_consensus.left_anchor_input_reads,
            event_consensus.right_anchor_input_reads,
            event_consensus.partial_context_input_reads,
            event_consensus.input_event_reads);
    if (!has_precise_or_full_context && !has_identifiable_bilateral_partial_context) {
        return "PRESEG_NO_PRECISE_OR_FULL_CONTEXT";
    }
    return {};
}

std::vector<Pipeline::AnchorSeedBin> Pipeline::collect_anchor_seed_bins(
    const std::string& query,
    int32_t breakpoint,
    int32_t ref_window_start,
    const std::string& ref_window,
    bool is_left) const {
    std::unordered_map<uint64_t, std::vector<int32_t>> ref_hits;
    ref_hits.reserve(ref_window.size());
    for_each_valid_kmer(ref_window, kEventSegmentationSeedK, [&](int32_t pos, uint64_t key) {
        ref_hits[key].push_back(pos);
    });

    std::unordered_map<int32_t, AnchorSeedBin> bins_by_start;
    bins_by_start.reserve(ref_window.size());
    for_each_valid_kmer(query, kEventSegmentationSeedK, [&](int32_t query_pos, uint64_t key) {
        const auto it = ref_hits.find(key);
        if (it == ref_hits.end()) {
            return;
        }
        for (int32_t ref_pos : it->second) {
            const int32_t ref_start = ref_window_start + ref_pos - query_pos;
            const int32_t ref_bin_start =
                (ref_start / kEventSegmentationSeedBinBp) * kEventSegmentationSeedBinBp;
            auto& bin = bins_by_start[ref_bin_start];
            bin.ref_bin_start = ref_bin_start;
            ++bin.support;

            const int32_t candidate_breakpoint = is_left
                ? (ref_start + static_cast<int32_t>(query.size()))
                : ref_start;
            bin.best_breakpoint_delta = std::min(
                bin.best_breakpoint_delta,
                std::abs(candidate_breakpoint - breakpoint));
        }
    });

    std::vector<AnchorSeedBin> bins;
    bins.reserve(bins_by_start.size());
    for (const auto& entry : bins_by_start) {
        bins.push_back(entry.second);
    }
    std::sort(bins.begin(), bins.end(), [](const AnchorSeedBin& lhs, const AnchorSeedBin& rhs) {
        if (lhs.support != rhs.support) {
            return lhs.support > rhs.support;
        }
        if (lhs.best_breakpoint_delta != rhs.best_breakpoint_delta) {
            return lhs.best_breakpoint_delta < rhs.best_breakpoint_delta;
        }
        return lhs.ref_bin_start < rhs.ref_bin_start;
    });
    if (static_cast<int32_t>(bins.size()) > kEventSegmentationSeedTopBins) {
        bins.resize(static_cast<size_t>(kEventSegmentationSeedTopBins));
    }
    return bins;
}

EventSegmentation Pipeline::segment_event_consensus(
    const ComponentCall& component,
    const EventReadEvidence& event_evidence,
    const EventConsensus& event_consensus) const {
    EventSegmentation segmentation;
    EventSegmentationSearchStats stats;
    last_segmentation_stats_ = {};
    const auto finish = [&](EventSegmentation result) {
        last_segmentation_stats_ = stats;
        return result;
    };

    if (!event_consensus.qc_pass || event_consensus.consensus_seq.empty()) {
        segmentation.qc_reason = "NO_EVENT_CONSENSUS_TO_SEGMENT";
        return finish(segmentation);
    }
    if (!tsd_detector_.can_fetch_reference()) {
        segmentation.qc_reason = "REFERENCE_FETCH_UNAVAILABLE";
        return finish(segmentation);
    }

    const std::string& consensus = event_consensus.consensus_seq;
    const int32_t consensus_len = static_cast<int32_t>(consensus.size());
    const auto retry_reverse_complement_segmentation = [&]() {
        EventSegmentation retry_segmentation;
        if (event_consensus.qc_reason == "PASS_EVENT_CONSENSUS_REVCOMP_RETRY") {
            return retry_segmentation;
        }
        EventConsensus reverse_consensus = event_consensus;
        reverse_consensus.consensus_seq = reverse_complement_acgt(consensus);
        reverse_consensus.consensus_len =
            static_cast<int32_t>(reverse_consensus.consensus_seq.size());
        if (reverse_consensus.consensus_seq.empty() ||
            reverse_consensus.consensus_seq == consensus) {
            return retry_segmentation;
        }
        reverse_consensus.qc_reason = "PASS_EVENT_CONSENSUS_REVCOMP_RETRY";
        stats.reverse_complement_retries += 1;
        retry_segmentation = segment_event_consensus(
            component,
            event_evidence,
            reverse_consensus);
        const EventSegmentationSearchStats retry_stats = last_segmentation_stats_;
        stats.edit_distance_calls += retry_stats.edit_distance_calls;
        stats.edit_distance_cache_hits += retry_stats.edit_distance_cache_hits;
        stats.edit_distance_cache_misses += retry_stats.edit_distance_cache_misses;
        stats.seed_bins_total += retry_stats.seed_bins_total;
        stats.paired_searches += retry_stats.paired_searches;
        stats.one_sided_searches += retry_stats.one_sided_searches;
        stats.endpoint_slack_searches += retry_stats.endpoint_slack_searches;
        stats.reverse_complement_retries += retry_stats.reverse_complement_retries;
        stats.segmentation_cache_hits += retry_stats.segmentation_cache_hits;
        stats.segmentation_cache_misses += retry_stats.segmentation_cache_misses;
        last_segmentation_stats_ = stats;
        return retry_segmentation.pass ? retry_segmentation : EventSegmentation{};
    };
    const int32_t min_one_sided_len =
        kEventSegmentationMinFlankAlignBp + kEventSegmentationMinInsertBp;
    if (consensus_len < min_one_sided_len) {
        segmentation.qc_reason = "EVENT_CONSENSUS_TOO_SHORT";
        return finish(segmentation);
    }

    const int32_t paired_max_flank_query_len = std::min(
        kEventSegmentationMaxFlankQueryBp,
        consensus_len - kEventSegmentationMinFlankAlignBp - kEventSegmentationMinInsertBp);
    const int32_t one_sided_max_flank_query_len = std::min(
        kEventSegmentationMaxFlankQueryBp,
        consensus_len - kEventSegmentationMinInsertBp);

    const int32_t bp_left = std::min(event_evidence.bp_left, event_evidence.bp_right);
    const int32_t bp_right = std::max(event_evidence.bp_left, event_evidence.bp_right);
    if (component.chrom.empty() || bp_left < 0 || bp_right < 0) {
        segmentation.qc_reason = "INVALID_EVENT_BREAKPOINTS";
        return finish(segmentation);
    }

    const auto emit_unplaced_insert_segmentation = [&]() {
        EventSegmentation fallback;
        if (event_consensus.qc_reason == "PASS_EVENT_CONSENSUS_REVCOMP_RETRY") {
            return fallback;
        }
        const bool strong_unopposed_alt =
            event_evidence.ref_span_reads == 0 &&
            event_evidence.alt_struct_reads >= 8 &&
            event_evidence.alt_struct_reads >=
                (event_consensus.input_event_reads * 2);
        const bool anchored_consensus =
            event_consensus.qc_pass &&
            event_consensus.left_anchor_input_reads >= 2 &&
            event_consensus.right_anchor_input_reads >= 2 &&
            event_consensus.partial_context_input_reads >= 4;
        if (!strong_unopposed_alt || !anchored_consensus) {
            return fallback;
        }

        fallback.insert_seq = consensus;
        fallback.left_ref_start = bp_left;
        fallback.left_ref_end = bp_left;
        fallback.right_ref_start = bp_right;
        fallback.right_ref_end = bp_right;
        fallback.pass = true;
        fallback.qc_reason = "PASS_EVENT_SEGMENTATION_UNPLACED_INSERT";
        return fallback;
    };

    const int32_t left_window_start = std::max(
        0,
        bp_left - kEventSegmentationBreakpointSlackBp - one_sided_max_flank_query_len);
    const int32_t left_window_end = std::max(
        left_window_start + 1,
        bp_left + kEventSegmentationBreakpointSlackBp);
    const int32_t right_window_start = std::max(
        0,
        bp_right - kEventSegmentationBreakpointSlackBp);
    const int32_t right_window_end = std::max(
        right_window_start + 1,
        bp_right + kEventSegmentationBreakpointSlackBp + one_sided_max_flank_query_len);

    const std::string left_ref_window = tsd_detector_.fetch_window(
        component.chrom,
        left_window_start,
        left_window_end);
    const std::string right_ref_window = tsd_detector_.fetch_window(
        component.chrom,
        right_window_start,
        right_window_end);
    if (left_ref_window.empty() || right_ref_window.empty()) {
        segmentation.qc_reason = "REFERENCE_WINDOW_FETCH_FAILED";
        return finish(segmentation);
    }

    struct EditMemoKey {
        int32_t query_start = 0;
        int32_t ref_start = 0;
        int32_t align_len = 0;

        bool operator==(const EditMemoKey& other) const {
            return query_start == other.query_start &&
                   ref_start == other.ref_start &&
                   align_len == other.align_len;
        }
    };
    struct EditMemoKeyHash {
        size_t operator()(const EditMemoKey& key) const {
            size_t seed = static_cast<size_t>(key.query_start);
            seed ^= static_cast<size_t>(key.ref_start) + 0x9e3779b97f4a7c15ULL +
                    (seed << 6) + (seed >> 2);
            seed ^= static_cast<size_t>(key.align_len) + 0x9e3779b97f4a7c15ULL +
                    (seed << 6) + (seed >> 2);
            return seed;
        }
    };
    struct EditMemoValue {
        bool pass = false;
        double identity = 0.0;
    };

    EditDistanceWorkspace edit_workspace;
    std::unordered_map<EditMemoKey, EditMemoValue, EditMemoKeyHash> edit_memo;
    const auto collect_candidates = [&](
        bool is_left,
        int32_t breakpoint,
        int32_t ref_window_start,
        const std::string& ref_window,
        int32_t flank_query_len,
        int32_t query_endpoint_slack,
        bool is_one_sided_search) {
        EventFlankSearchResult result;
        if (ref_window.empty() || flank_query_len < kEventSegmentationMinFlankAlignBp) {
            return result;
        }
        if (is_one_sided_search) {
            stats.one_sided_searches += 1;
        } else {
            stats.paired_searches += 1;
        }
        if (query_endpoint_slack > 0) {
            stats.endpoint_slack_searches += 1;
        }

        const int32_t ref_window_end =
            ref_window_start + static_cast<int32_t>(ref_window.size());
        const int32_t search_lo = std::max(0, breakpoint - kEventSegmentationBreakpointSlackBp);
        const int32_t search_hi = breakpoint + kEventSegmentationBreakpointSlackBp;
        const int32_t query_seed_len = std::min(
            flank_query_len + query_endpoint_slack,
            consensus_len);
        const int32_t query_seed_start = is_left ? 0 : (consensus_len - query_seed_len);
        const std::string seed_query = is_left
            ? consensus.substr(0, static_cast<size_t>(query_seed_len))
            : consensus.substr(static_cast<size_t>(query_seed_start));
        const std::vector<AnchorSeedBin> seed_bins = collect_anchor_seed_bins(
            seed_query,
            breakpoint,
            ref_window_start,
            ref_window,
            is_left);
        stats.seed_bins_total += static_cast<int64_t>(seed_bins.size());
        if (seed_bins.empty()) {
            return result;
        }

        for (int32_t align_len = flank_query_len;
             align_len >= kEventSegmentationMinFlankAlignBp;
             --align_len) {
            const int32_t max_query_start = consensus_len - align_len;
            if (max_query_start < 0) {
                continue;
            }
            const int32_t max_edits = max_edits_for_identity_threshold(
                align_len,
                align_len,
                kEventSegmentationMinFlankIdentity);

            const int32_t endpoint_slack = std::max(0, std::min(
                query_endpoint_slack,
                max_query_start));
            const int32_t query_lo = is_left
                ? 0
                : std::max(0, max_query_start - endpoint_slack);
            const int32_t query_hi = is_left
                ? std::min(max_query_start, endpoint_slack)
                : max_query_start;

            std::vector<EventFlankPlacement> placements;
            placements.reserve(
                seed_bins.size() *
                static_cast<size_t>(kEventSegmentationSeedBinBp) *
                static_cast<size_t>(std::max(1, query_hi - query_lo + 1)));

            for (int32_t query_start = query_lo; query_start <= query_hi; ++query_start) {
                const std::string_view query(
                    consensus.data() + query_start,
                    static_cast<size_t>(align_len));
                const int32_t endpoint_offset = is_left
                    ? query_start
                    : (consensus_len - (query_start + align_len));
                const int32_t seed_query_offset = query_start - query_seed_start;

                for (const auto& seed_bin : seed_bins) {
                    for (int32_t bin_offset = 0;
                         bin_offset < kEventSegmentationSeedBinBp;
                         ++bin_offset) {
                        const int32_t ref_start =
                            seed_bin.ref_bin_start + seed_query_offset + bin_offset;
                        const int32_t ref_end = ref_start + align_len;
                        const int32_t candidate_breakpoint = is_left ? ref_end : ref_start;
                        if (candidate_breakpoint < search_lo || candidate_breakpoint > search_hi) {
                            continue;
                        }
                        if (ref_start < ref_window_start || ref_end > ref_window_end) {
                            continue;
                        }

                        const size_t offset = static_cast<size_t>(ref_start - ref_window_start);
                        const std::string_view ref_seq(
                            ref_window.data() + offset,
                            static_cast<size_t>(align_len));
                        double identity = 0.0;
                        const EditMemoKey memo_key{query_start, ref_start, align_len};
                        bool pass_identity = false;
                        const auto memo_it = edit_memo.find(memo_key);
                        if (memo_it != edit_memo.end()) {
                            stats.edit_distance_cache_hits += 1;
                            pass_identity = memo_it->second.pass;
                            identity = memo_it->second.identity;
                        } else {
                            stats.edit_distance_cache_misses += 1;
                            stats.edit_distance_calls += 1;
                            pass_identity = edit_identity_if_at_least_fast_impl(
                                query,
                                ref_seq,
                                max_edits,
                                edit_workspace,
                                identity);
                            edit_memo.emplace(
                                memo_key,
                                EditMemoValue{pass_identity, identity});
                        }
                        if (!pass_identity) {
                            continue;
                        }

                        EventFlankPlacement placement;
                        placement.query_start = query_start;
                        placement.query_end = query_start + align_len;
                        placement.ref_start = ref_start;
                        placement.ref_end = ref_end;
                        placement.align_len = align_len;
                        placement.identity = identity;
                        placement.breakpoint_delta = std::abs(candidate_breakpoint - breakpoint);
                        placement.endpoint_offset = endpoint_offset;
                        placements.push_back(placement);
                    }
                }
            }

            if (placements.empty()) {
                continue;
            }

            std::sort(
                placements.begin(),
                placements.end(),
                better_event_flank_placement);
            const double best_identity = placements.front().identity;
            for (const auto& placement : placements) {
                if ((placement.identity + kEventSegmentationUniquenessMargin) < best_identity) {
                    break;
                }
                result.candidates.push_back(placement);
            }
        }
        std::sort(
            result.candidates.begin(),
            result.candidates.end(),
            better_event_flank_placement);
        return result;
    };

    EventFlankSearchResult left_search = collect_candidates(
        true,
        bp_left,
        left_window_start,
        left_ref_window,
        paired_max_flank_query_len,
        0,
        false);
    if (left_search.candidates.empty() && kEventSegmentationEndpointSlackBp > 0) {
        left_search = collect_candidates(
            true,
            bp_left,
            left_window_start,
            left_ref_window,
            paired_max_flank_query_len,
            kEventSegmentationEndpointSlackBp,
            false);
    }

    EventFlankSearchResult right_search = collect_candidates(
        false,
        bp_right,
        right_window_start,
        right_ref_window,
        paired_max_flank_query_len,
        0,
        false);
    if (right_search.candidates.empty() && kEventSegmentationEndpointSlackBp > 0) {
        right_search = collect_candidates(
            false,
            bp_right,
            right_window_start,
            right_ref_window,
            paired_max_flank_query_len,
            kEventSegmentationEndpointSlackBp,
            false);
    }

    const auto emit_left_sided_segmentation = [&](const EventFlankPlacement& best_left) {
        const int32_t insert_start = best_left.query_end;
        const int32_t insert_len = consensus_len - insert_start;
        if (insert_len < kEventSegmentationMinInsertBp) {
            segmentation.qc_reason = "EMPTY_EVENT_INSERT_SEGMENT";
            return;
        }

        segmentation.left_flank_seq = consensus.substr(
            static_cast<size_t>(best_left.query_start),
            static_cast<size_t>(best_left.align_len));
        segmentation.insert_seq = consensus.substr(
            static_cast<size_t>(insert_start),
            static_cast<size_t>(insert_len));
        segmentation.left_ref_start = best_left.ref_start;
        segmentation.left_ref_end = best_left.ref_end;
        segmentation.right_ref_start = bp_right;
        segmentation.right_ref_end = bp_right;
        segmentation.left_flank_align_len = best_left.align_len;
        segmentation.right_flank_align_len = 0;
        segmentation.left_flank_identity = best_left.identity;
        segmentation.right_flank_identity = 0.0;
        segmentation.pass = true;
        segmentation.qc_reason = "PASS_EVENT_SEGMENTATION_ONE_SIDED_LEFT";
    };

    const auto emit_right_sided_segmentation = [&](const EventFlankPlacement& best_right) {
        const int32_t insert_len = best_right.query_start;
        if (insert_len < kEventSegmentationMinInsertBp) {
            segmentation.qc_reason = "EMPTY_EVENT_INSERT_SEGMENT";
            return;
        }

        segmentation.insert_seq = consensus.substr(
            0,
            static_cast<size_t>(insert_len));
        segmentation.right_flank_seq = consensus.substr(
            static_cast<size_t>(best_right.query_start),
            static_cast<size_t>(best_right.align_len));
        segmentation.left_ref_start = bp_left;
        segmentation.left_ref_end = bp_left;
        segmentation.right_ref_start = best_right.ref_start;
        segmentation.right_ref_end = best_right.ref_end;
        segmentation.left_flank_align_len = 0;
        segmentation.right_flank_align_len = best_right.align_len;
        segmentation.left_flank_identity = 0.0;
        segmentation.right_flank_identity = best_right.identity;
        segmentation.pass = true;
        segmentation.qc_reason = "PASS_EVENT_SEGMENTATION_ONE_SIDED_RIGHT";
    };

    const auto collect_one_sided_search = [&](bool is_left) {
        EventFlankSearchResult search = collect_candidates(
            is_left,
            is_left ? bp_left : bp_right,
            is_left ? left_window_start : right_window_start,
            is_left ? left_ref_window : right_ref_window,
            one_sided_max_flank_query_len,
            0,
            true);
        if (search.candidates.empty() && kEventSegmentationEndpointSlackBp > 0) {
            search = collect_candidates(
                is_left,
                is_left ? bp_left : bp_right,
                is_left ? left_window_start : right_window_start,
                is_left ? left_ref_window : right_ref_window,
                one_sided_max_flank_query_len,
                kEventSegmentationEndpointSlackBp,
                true);
        }
        return search;
    };

    const auto better_one_sided = [](
                                      const EventFlankPlacement& lhs,
                                      const EventFlankPlacement& rhs) {
        if (lhs.identity != rhs.identity) {
            return lhs.identity > rhs.identity;
        }
        if (lhs.breakpoint_delta != rhs.breakpoint_delta) {
            return lhs.breakpoint_delta < rhs.breakpoint_delta;
        }
        if (lhs.align_len != rhs.align_len) {
            return lhs.align_len > rhs.align_len;
        }
        if (lhs.endpoint_offset != rhs.endpoint_offset) {
            return lhs.endpoint_offset < rhs.endpoint_offset;
        }
        return lhs.ref_start < rhs.ref_start;
    };

    const auto better_pair = [&](const EventFlankPlacement& lhs_left,
                                 const EventFlankPlacement& lhs_right,
                                 const EventFlankPlacement& rhs_left,
                                 const EventFlankPlacement& rhs_right) {
        const double lhs_min_identity = std::min(lhs_left.identity, lhs_right.identity);
        const double rhs_min_identity = std::min(rhs_left.identity, rhs_right.identity);
        if (lhs_min_identity != rhs_min_identity) {
            return lhs_min_identity > rhs_min_identity;
        }

        const double lhs_mean_identity = 0.5 * (lhs_left.identity + lhs_right.identity);
        const double rhs_mean_identity = 0.5 * (rhs_left.identity + rhs_right.identity);
        if (lhs_mean_identity != rhs_mean_identity) {
            return lhs_mean_identity > rhs_mean_identity;
        }

        const int32_t lhs_delta = lhs_left.breakpoint_delta + lhs_right.breakpoint_delta;
        const int32_t rhs_delta = rhs_left.breakpoint_delta + rhs_right.breakpoint_delta;
        if (lhs_delta != rhs_delta) {
            return lhs_delta < rhs_delta;
        }

        const int32_t lhs_total = lhs_left.align_len + lhs_right.align_len;
        const int32_t rhs_total = rhs_left.align_len + rhs_right.align_len;
        if (lhs_total != rhs_total) {
            return lhs_total > rhs_total;
        }

        if (lhs_left.ref_start != rhs_left.ref_start) {
            return lhs_left.ref_start < rhs_left.ref_start;
        }
        return lhs_right.ref_start < rhs_right.ref_start;
    };

    const auto canonical_pair = [&](const EventFlankPlacement& left,
                                    const EventFlankPlacement& right) {
        if (left.query_end + kEventSegmentationMinInsertBp > right.query_start) {
            return false;
        }
        FinalBoundaryInput input;
        input.left_ref_start = left.ref_start;
        input.left_ref_end = left.ref_end;
        input.right_ref_start = right.ref_start;
        input.right_ref_end = right.ref_end;
        input.tsd_min_len = std::max(1, config_.tsd_min_len);
        input.tsd_max_len = std::max(input.tsd_min_len, config_.tsd_max_len);
        return check_boundary_consistency(input).pass;
    };

    const auto collect_tethered_opposite_flank = [&](bool missing_is_left,
                                                      int32_t anchored_breakpoint) {
        EventFlankSearchResult search;
        if (anchored_breakpoint < 0) {
            return search;
        }
        const int32_t ref_window_start = missing_is_left
            ? std::max(
                  0,
                  anchored_breakpoint - kEventSegmentationBreakpointSlackBp -
                      one_sided_max_flank_query_len)
            : std::max(0, anchored_breakpoint - kEventSegmentationBreakpointSlackBp);
        const int32_t ref_window_end = missing_is_left
            ? std::max(
                  ref_window_start + 1,
                  anchored_breakpoint + kEventSegmentationBreakpointSlackBp)
            : std::max(
                  ref_window_start + 1,
                  anchored_breakpoint + kEventSegmentationBreakpointSlackBp +
                      one_sided_max_flank_query_len);
        const std::string ref_window = tsd_detector_.fetch_window(
            component.chrom,
            ref_window_start,
            ref_window_end);
        if (ref_window.empty()) {
            return search;
        }

        search = collect_candidates(
            missing_is_left,
            anchored_breakpoint,
            ref_window_start,
            ref_window,
            one_sided_max_flank_query_len,
            0,
            false);
        if (search.candidates.empty() && kEventSegmentationEndpointSlackBp > 0) {
            search = collect_candidates(
                missing_is_left,
                anchored_breakpoint,
                ref_window_start,
                ref_window,
                one_sided_max_flank_query_len,
                kEventSegmentationEndpointSlackBp,
                false);
        }
        return search;
    };

    const auto try_tethered_pair = [&](const EventFlankSearchResult& anchored_search,
                                       bool anchor_is_left,
                                       EventFlankPlacement& rescued_left,
                                       EventFlankPlacement& rescued_right) {
        constexpr int32_t kPlacementDiagonalToleranceBp = 5;
        if (anchored_search.candidates.empty()) {
            return false;
        }

        const EventFlankPlacement& anchor = anchored_search.candidates.front();
        const int32_t anchor_boundary = anchor_is_left ? anchor.ref_end : anchor.ref_start;
        const int32_t anchor_diagonal = anchor.ref_start - anchor.query_start;
        for (const auto& alternative : anchored_search.candidates) {
            if (alternative.identity + kEventSegmentationUniquenessMargin <
                anchor.identity) {
                break;
            }
            const int32_t alternative_diagonal =
                alternative.ref_start - alternative.query_start;
            if (std::abs(alternative_diagonal - anchor_diagonal) >
                kPlacementDiagonalToleranceBp) {
                return false;
            }
        }

        const bool missing_is_left = !anchor_is_left;
        const EventFlankSearchResult opposite_search =
            collect_tethered_opposite_flank(missing_is_left, anchor_boundary);
        bool found = false;
        EventFlankPlacement best_left;
        EventFlankPlacement best_right;
        for (const auto& opposite : opposite_search.candidates) {
            const EventFlankPlacement& left = anchor_is_left ? anchor : opposite;
            const EventFlankPlacement& right = anchor_is_left ? opposite : anchor;
            if (!canonical_pair(left, right)) {
                continue;
            }
            if (!found || better_pair(left, right, best_left, best_right)) {
                best_left = left;
                best_right = right;
                found = true;
            }
        }
        if (!found) {
            return false;
        }

        const int32_t best_opposite_diagonal = anchor_is_left
            ? (best_right.ref_start - best_right.query_start)
            : (best_left.ref_start - best_left.query_start);
        const double best_opposite_identity = anchor_is_left
            ? best_right.identity
            : best_left.identity;
        for (const auto& opposite : opposite_search.candidates) {
            if (opposite.identity + kEventSegmentationUniquenessMargin <
                best_opposite_identity) {
                break;
            }
            const EventFlankPlacement& left = anchor_is_left ? anchor : opposite;
            const EventFlankPlacement& right = anchor_is_left ? opposite : anchor;
            if (!canonical_pair(left, right)) {
                continue;
            }
            const EventFlankPlacement& candidate = anchor_is_left ? right : left;
            const int32_t diagonal = candidate.ref_start - candidate.query_start;
            if (std::abs(diagonal - best_opposite_diagonal) >
                kPlacementDiagonalToleranceBp) {
                return false;
            }
        }

        rescued_left = best_left;
        rescued_right = best_right;
        return true;
    };

    const auto emit_best_one_sided = [&](bool allow_left, bool allow_right) {
        EventFlankSearchResult one_sided_left;
        EventFlankSearchResult one_sided_right;
        if (allow_left) {
            one_sided_left = collect_one_sided_search(true);
        }
        if (allow_right) {
            one_sided_right = collect_one_sided_search(false);
        }

        const bool have_left = !one_sided_left.candidates.empty();
        const bool have_right = !one_sided_right.candidates.empty();
        if (!have_left && !have_right) {
            return false;
        }
        if (have_left && !have_right) {
            emit_left_sided_segmentation(one_sided_left.candidates.front());
            return segmentation.pass;
        }
        if (!have_left && have_right) {
            emit_right_sided_segmentation(one_sided_right.candidates.front());
            return segmentation.pass;
        }
        if (better_one_sided(
                one_sided_left.candidates.front(),
                one_sided_right.candidates.front())) {
            emit_left_sided_segmentation(one_sided_left.candidates.front());
        } else {
            emit_right_sided_segmentation(one_sided_right.candidates.front());
        }
        return segmentation.pass;
    };

    EventFlankPlacement best_left;
    EventFlankPlacement best_right;
    bool have_pair = false;

    for (const auto& left_candidate : left_search.candidates) {
        for (const auto& right_candidate : right_search.candidates) {
            if (left_candidate.query_end + kEventSegmentationMinInsertBp >
                right_candidate.query_start) {
                continue;
            }

            if (!have_pair ||
                better_pair(left_candidate, right_candidate, best_left, best_right)) {
                best_left = left_candidate;
                best_right = right_candidate;
                have_pair = true;
            }
        }
    }

    bool tethered_pair = false;
    if (!have_pair) {
        EventFlankPlacement rescued_left;
        EventFlankPlacement rescued_right;
        if (try_tethered_pair(
                left_search,
                true,
                rescued_left,
                rescued_right)) {
            best_left = rescued_left;
            best_right = rescued_right;
            have_pair = true;
            tethered_pair = true;
        }

        EventFlankPlacement reverse_rescued_left;
        EventFlankPlacement reverse_rescued_right;
        if (try_tethered_pair(
                right_search,
                false,
                reverse_rescued_left,
                reverse_rescued_right) &&
            (!have_pair || better_pair(
                reverse_rescued_left,
                reverse_rescued_right,
                best_left,
                best_right))) {
            best_left = reverse_rescued_left;
            best_right = reverse_rescued_right;
            have_pair = true;
            tethered_pair = true;
        }
    }

    if (!have_pair) {
        if (emit_best_one_sided(true, true)) {
            return finish(segmentation);
        }
        if (const EventSegmentation retry = retry_reverse_complement_segmentation();
            retry.pass) {
            return finish(retry);
        }
        if (const EventSegmentation fallback = emit_unplaced_insert_segmentation();
            fallback.pass) {
            return finish(fallback);
        }
        segmentation.qc_reason = "NO_TRIPARTITE_EVENT_SEGMENTATION";
        return finish(segmentation);
    }

    const int32_t insert_start = best_left.query_end;
    const int32_t insert_len = best_right.query_start - best_left.query_end;
    if (insert_len < kEventSegmentationMinInsertBp) {
        segmentation.qc_reason = "EMPTY_EVENT_INSERT_SEGMENT";
        return finish(segmentation);
    }

    segmentation.left_flank_seq = consensus.substr(
        static_cast<size_t>(best_left.query_start),
        static_cast<size_t>(best_left.align_len));
    segmentation.insert_seq = consensus.substr(
        static_cast<size_t>(insert_start),
        static_cast<size_t>(insert_len));
    segmentation.right_flank_seq = consensus.substr(
        static_cast<size_t>(best_right.query_start),
        static_cast<size_t>(best_right.align_len));
    segmentation.left_ref_start = best_left.ref_start;
    segmentation.left_ref_end = best_left.ref_end;
    segmentation.right_ref_start = best_right.ref_start;
    segmentation.right_ref_end = best_right.ref_end;
    segmentation.left_flank_align_len = best_left.align_len;
    segmentation.right_flank_align_len = best_right.align_len;
    segmentation.left_flank_identity = best_left.identity;
    segmentation.right_flank_identity = best_right.identity;
    segmentation.pass = true;
    segmentation.qc_reason = tethered_pair
        ? "PASS_EVENT_SEGMENTATION_TETHERED_FLANK"
        : "PASS_EVENT_SEGMENTATION";
    return finish(segmentation);
}

EventSegmentation Pipeline::segment_event_consensus_cached(
    const ComponentCall& component,
    const EventReadEvidence& event_evidence,
    const EventConsensus& event_consensus,
    std::unordered_map<std::string, EventSegmentation>& segmentation_cache) const {
    const std::string key = event_segmentation_cache_key(
        component,
        event_evidence,
        event_consensus);
    const auto cached = segmentation_cache.find(key);
    if (cached != segmentation_cache.end()) {
        last_segmentation_stats_ = {};
        last_segmentation_stats_.segmentation_cache_hits = 1;
        return cached->second;
    }
    EventSegmentation segmentation = segment_event_consensus(
        component,
        event_evidence,
        event_consensus);
    last_segmentation_stats_.segmentation_cache_misses += 1;
    segmentation_cache.emplace(key, segmentation);
    return segmentation;
}

EventSegmentationSearchStats Pipeline::last_segmentation_stats() const {
    return last_segmentation_stats_;
}

EventSegmentationEvidence Pipeline::analyze_event_segmentation(
    const EventConsensus& event_consensus,
    const EventSegmentation& event_segmentation) const {
    return analyze_event_segmentation_for_test(
        event_consensus.qc_pass,
        event_segmentation);
}

TEAlignmentEvidence Pipeline::align_insert_seq_to_te(
    const EventSegmentation& event_segmentation) const {
    if (!event_segmentation.pass) {
        TEAlignmentEvidence evidence;
        evidence.qc_reason = "NO_EVENT_SEGMENTATION_FOR_TE_ALIGNMENT";
        return evidence;
    }
    return te_classifier_module_.align_insert_sequence(event_segmentation.insert_seq);
}

TEAlignmentEvidence Pipeline::align_insert_seq_to_te_cached(
    const EventSegmentation& event_segmentation,
    std::unordered_map<std::string, TEAlignmentEvidence>& te_alignment_cache) const {
    if (!event_segmentation.pass || event_segmentation.insert_seq.empty()) {
        return align_insert_seq_to_te(event_segmentation);
    }
    const auto cached = te_alignment_cache.find(event_segmentation.insert_seq);
    if (cached != te_alignment_cache.end()) {
        return cached->second;
    }
    TEAlignmentEvidence evidence = align_insert_seq_to_te(event_segmentation);
    te_alignment_cache.emplace(event_segmentation.insert_seq, evidence);
    return evidence;
}

ClipInsertConcordanceEvidence Pipeline::analyze_clip_insert_concordance(
    const EventReadEvidence& event_evidence,
    const EventSegmentation& event_segmentation,
    const std::vector<InsertionFragment>& fragments) const {
    ClipInsertConcordanceEvidence evidence;
    if (!event_segmentation.pass || event_segmentation.insert_seq.empty()) {
        evidence.qc = "NO_EVENT_SEGMENTATION_FOR_CLIP_INSERT_CONCORDANCE";
        return evidence;
    }

    constexpr int32_t kMinClipMatchLen = 50;
    constexpr int32_t kMaxClipMatchLen = 300;
    constexpr double kMinClipIdentity = 0.85;
    const std::string insert_seq = upper_acgt(event_segmentation.insert_seq);
    if (static_cast<int32_t>(insert_seq.size()) < kMinClipMatchLen) {
        evidence.qc = "CLIP_INSERT_CONCORDANCE_INSERT_TOO_SHORT";
        return evidence;
    }

    const std::unordered_set<std::string> support_qnames(
        event_evidence.support_qnames.begin(),
        event_evidence.support_qnames.end());
    std::unordered_set<std::string> full_insert_qnames;
    std::unordered_set<std::string> left_clip_qnames;
    std::unordered_set<std::string> right_clip_qnames;

    for (const auto& fragment : fragments) {
        if (fragment.read_id.empty() ||
            support_qnames.find(fragment.read_id) == support_qnames.end() ||
            !fragment_is_local_to_event(
                fragment,
                event_evidence.bp_left,
                event_evidence.bp_right)) {
            continue;
        }

        if (fragment.source == InsertionFragmentSource::kCigarInsertion ||
            fragment.source == InsertionFragmentSource::kSplitSa) {
            if (static_cast<int32_t>(fragment.sequence.size()) >= kMinClipMatchLen) {
                full_insert_qnames.insert(fragment.read_id);
            }
            continue;
        }

        if (!is_softclip_source(fragment.source) ||
            is_low_complexity_softclip(fragment, config_) ||
            is_low_quality_softclip_anchor(fragment, config_)) {
            continue;
        }

        const std::string clip_seq = upper_acgt(fragment.sequence);
        const int32_t match_len = std::min({
            static_cast<int32_t>(clip_seq.size()),
            static_cast<int32_t>(insert_seq.size()),
            kMaxClipMatchLen});
        if (match_len < kMinClipMatchLen) {
            continue;
        }

        const std::string clipped_window = clip_seq.substr(0, static_cast<size_t>(match_len));
        std::string insert_window;
        if (fragment.source == InsertionFragmentSource::kClipRefLeft) {
            insert_window = insert_seq.substr(0, static_cast<size_t>(match_len));
        } else {
            insert_window = insert_seq.substr(
                insert_seq.size() - static_cast<size_t>(match_len),
                static_cast<size_t>(match_len));
        }

        const double identity = fixed_window_edit_identity(clipped_window, insert_window);
        if (identity + 1e-9 < kMinClipIdentity) {
            continue;
        }

        if (fragment.source == InsertionFragmentSource::kClipRefLeft) {
            left_clip_qnames.insert(fragment.read_id);
            evidence.max_left_identity = std::max(evidence.max_left_identity, identity);
        } else {
            right_clip_qnames.insert(fragment.read_id);
            evidence.max_right_identity = std::max(evidence.max_right_identity, identity);
        }
    }

    evidence.full_insert_reads = static_cast<int32_t>(full_insert_qnames.size());
    evidence.left_clip_reads = static_cast<int32_t>(left_clip_qnames.size());
    evidence.right_clip_reads = static_cast<int32_t>(right_clip_qnames.size());

    const bool bilateral_clip =
        evidence.left_clip_reads >= 1 && evidence.right_clip_reads >= 1;
    const bool strong_one_sided_clip =
        evidence.left_clip_reads >= 2 || evidence.right_clip_reads >= 2;
    evidence.pass =
        evidence.full_insert_reads >= 1 &&
        (bilateral_clip || strong_one_sided_clip);
    evidence.qc = evidence.pass
        ? "PASS_CLIP_INSERT_CONCORDANCE"
        : "NO_CLIP_INSERT_CONCORDANCE";
    return evidence;
}
