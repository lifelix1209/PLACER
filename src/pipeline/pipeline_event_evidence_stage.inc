EventReadEvidence Pipeline::collect_event_read_evidence(
    const ComponentCall& component,
    const std::vector<const bam1_t*>& local_records,
    const std::vector<ReadReferenceSpan>& read_spans,
    const std::vector<InsertionFragment>& fragments) const {
    const auto seed_bounds = infer_component_breakpoint_bounds(component);
    const auto bp_bounds = resolve_event_breakpoint_bounds(
        component,
        local_records,
        fragments,
        seed_bounds.first,
        seed_bounds.second);
    return collect_event_read_evidence_for_bounds(
        component,
        local_records,
        read_spans,
        fragments,
        seed_bounds.first,
        seed_bounds.second,
        bp_bounds.first,
        bp_bounds.second);
}

EventReadEvidence Pipeline::collect_event_read_evidence_for_bounds(
    const ComponentCall& component,
    const std::vector<const bam1_t*>& local_records,
    const std::vector<ReadReferenceSpan>& read_spans,
    const std::vector<InsertionFragment>& fragments,
    int32_t seed_left,
    int32_t seed_right,
    int32_t bp_left,
    int32_t bp_right) const {
    EventReadEvidence evidence;
    (void)seed_left;
    (void)seed_right;
    constexpr int32_t kAltSignalSlackBp = 25;
    constexpr int32_t kRefSignalSlackBp = 75;
    constexpr int32_t kCandidateClipSlackBp = 25;
    constexpr int32_t kCandidateClipSupplementMaxOffsetBp = 75;
    const int32_t left = std::min(bp_left, bp_right);
    const int32_t right = std::max(bp_left, bp_right);
    evidence.bp_left = left;
    evidence.bp_right = right;

    std::unordered_set<std::string> split_qnames;
    std::unordered_set<std::string> indel_qnames;
    std::unordered_set<std::string> left_clip_qnames;
    std::unordered_set<std::string> right_clip_qnames;
    std::unordered_set<std::string> raw_cigar_insert_qnames;
    std::unordered_set<std::string> nearby_candidate_left_clip_qnames;
    std::unordered_set<std::string> nearby_candidate_right_clip_qnames;

    // Alt support should stay hypothesis-specific, while clean ref-spanning
    // reads need a slightly wider exclusion band so nearby clipped noise does
    // not masquerade as reference support.
    const int32_t alt_signal_start = std::max(0, evidence.bp_left - kAltSignalSlackBp);
    const int32_t alt_signal_end = std::max(alt_signal_start, evidence.bp_right + kAltSignalSlackBp);
    const int32_t ref_signal_start = std::max(0, evidence.bp_left - kRefSignalSlackBp);
    const int32_t ref_signal_end = std::max(ref_signal_start, evidence.bp_right + kRefSignalSlackBp);
    const int32_t ref_span_start = alt_signal_start;
    const int32_t ref_span_end = alt_signal_end;

    std::vector<int32_t> nearby_component_clip_candidates;
    nearby_component_clip_candidates.reserve(component.breakpoint_candidates.size());
    for (const auto& candidate : component.breakpoint_candidates) {
        if (candidate.pos < 0) {
            continue;
        }
        const int32_t dist_to_hypothesis = (candidate.pos < left)
            ? (left - candidate.pos)
            : ((candidate.pos > right) ? (candidate.pos - right) : 0);
        if (dist_to_hypothesis <= kCandidateClipSupplementMaxOffsetBp) {
            nearby_component_clip_candidates.push_back(candidate.pos);
        }
    }
    std::sort(
        nearby_component_clip_candidates.begin(),
        nearby_component_clip_candidates.end());
    nearby_component_clip_candidates.erase(
        std::unique(
            nearby_component_clip_candidates.begin(),
            nearby_component_clip_candidates.end()),
        nearby_component_clip_candidates.end());

    const bool have_nearby_component_clip_candidates =
        !nearby_component_clip_candidates.empty();
    int32_t candidate_clip_start = 0;
    int32_t candidate_clip_end = 0;
    if (have_nearby_component_clip_candidates) {
        candidate_clip_start = std::max(
            0,
            nearby_component_clip_candidates.front() - kCandidateClipSlackBp);
        candidate_clip_end = std::max(
            candidate_clip_start,
            nearby_component_clip_candidates.back() + kCandidateClipSlackBp);
    }
    const auto matches_nearby_component_clip_candidate = [&](int32_t pos) {
        if (pos < 0) {
            return false;
        }
        for (const int32_t candidate_pos : nearby_component_clip_candidates) {
            if (std::abs(pos - candidate_pos) <= kCandidateClipSlackBp) {
                return true;
            }
        }
        return false;
    };

    for (size_t read_idx = 0; read_idx < local_records.size(); ++read_idx) {
        const bam1_t* record = local_records[read_idx];
        if (!record) {
            continue;
        }
        ReadView read(record);
        if (read.tid() != component.tid) {
            continue;
        }

        const std::string qname = record_qname(record);
        if (qname.empty()) {
            continue;
        }

        const LocalEventSignal signal = classify_local_event_signal(
            record,
            component.chrom,
            alt_signal_start,
            alt_signal_end);
        if (signal.split) {
            split_qnames.insert(qname);
        }
        if (signal.indel && read.mapq() == kInsertionCandidateRequiredMapq) {
            indel_qnames.insert(qname);
        }
        if (signal.max_raw_cigar_insert_len > 0) {
            raw_cigar_insert_qnames.insert(qname);
            evidence.max_raw_cigar_insert_len = std::max(
                evidence.max_raw_cigar_insert_len,
                signal.max_raw_cigar_insert_len);
        }
        if (signal.left_clip) {
            left_clip_qnames.insert(qname);
        }
        if (signal.right_clip) {
            right_clip_qnames.insert(qname);
        }

        if (have_nearby_component_clip_candidates) {
            const LocalEventSignal nearby_signal = classify_local_event_signal(
                record,
                component.chrom,
                candidate_clip_start,
                candidate_clip_end);
            if (nearby_signal.left_clip &&
                matches_nearby_component_clip_candidate(nearby_signal.left_clip_pos)) {
                nearby_candidate_left_clip_qnames.insert(qname);
            }
            if (nearby_signal.right_clip &&
                matches_nearby_component_clip_candidate(nearby_signal.right_clip_pos)) {
                nearby_candidate_right_clip_qnames.insert(qname);
            }
        }
    }

    for (const auto& fragment : fragments) {
        if (fragment.read_id.empty() ||
            fragment.ref_junc_pos < alt_signal_start ||
            fragment.ref_junc_pos > alt_signal_end) {
            continue;
        }
        switch (fragment.source) {
            case InsertionFragmentSource::kSplitSa:
                split_qnames.insert(fragment.read_id);
                break;
            case InsertionFragmentSource::kCigarInsertion:
                indel_qnames.insert(fragment.read_id);
                raw_cigar_insert_qnames.insert(fragment.read_id);
                evidence.max_raw_cigar_insert_len = std::max(
                    evidence.max_raw_cigar_insert_len,
                    fragment.length);
                break;
            case InsertionFragmentSource::kClipRefLeft:
                left_clip_qnames.insert(fragment.read_id);
                break;
            case InsertionFragmentSource::kClipRefRight:
                right_clip_qnames.insert(fragment.read_id);
                break;
            default:
                break;
        }

        if (have_nearby_component_clip_candidates &&
            matches_nearby_component_clip_candidate(fragment.ref_junc_pos)) {
            switch (fragment.source) {
                case InsertionFragmentSource::kClipRefLeft:
                    nearby_candidate_left_clip_qnames.insert(fragment.read_id);
                    break;
                case InsertionFragmentSource::kClipRefRight:
                    nearby_candidate_right_clip_qnames.insert(fragment.read_id);
                    break;
                default:
                    break;
            }
        }
    }

    const bool has_precise_nonclip_support =
        !split_qnames.empty() || !indel_qnames.empty();
    if (has_precise_nonclip_support &&
        !nearby_candidate_left_clip_qnames.empty() &&
        !nearby_candidate_right_clip_qnames.empty()) {
        left_clip_qnames.insert(
            nearby_candidate_left_clip_qnames.begin(),
            nearby_candidate_left_clip_qnames.end());
        right_clip_qnames.insert(
            nearby_candidate_right_clip_qnames.begin(),
            nearby_candidate_right_clip_qnames.end());
    }

    evidence.alt_split_reads = static_cast<int32_t>(split_qnames.size());
    evidence.alt_indel_reads = static_cast<int32_t>(indel_qnames.size());
    evidence.alt_left_clip_reads = static_cast<int32_t>(left_clip_qnames.size());
    evidence.alt_right_clip_reads = static_cast<int32_t>(right_clip_qnames.size());
    evidence.raw_cigar_insert_reads =
        static_cast<int32_t>(raw_cigar_insert_qnames.size());

    std::unordered_set<std::string> alt_qnames = split_qnames;
    alt_qnames.insert(indel_qnames.begin(), indel_qnames.end());
    const bool clip_support_has_partner =
        !alt_qnames.empty() ||
        (!left_clip_qnames.empty() && !right_clip_qnames.empty());
    if (clip_support_has_partner) {
        alt_qnames.insert(left_clip_qnames.begin(), left_clip_qnames.end());
        alt_qnames.insert(right_clip_qnames.begin(), right_clip_qnames.end());
    }

    evidence.alt_struct_reads = static_cast<int32_t>(alt_qnames.size());
    evidence.support_qnames = sorted_string_set(alt_qnames);

    std::unordered_set<std::string> ref_qnames;
    std::unordered_set<std::string> low_mapq_ref_qnames;
    for (size_t read_idx = 0; read_idx < read_spans.size(); ++read_idx) {
        const auto& span = read_spans[read_idx];
        if (!span.valid || span.tid != component.tid) {
            continue;
        }
        if (span.start > ref_span_start || span.end < ref_span_end) {
            continue;
        }
        if (read_idx >= local_records.size() || !local_records[read_idx]) {
            continue;
        }
        const bam1_t* record = local_records[read_idx];
        if (!is_primary_alignment(record)) {
            continue;
        }
        const std::string qname = record_qname(record);
        if (qname.empty() || alt_qnames.find(qname) != alt_qnames.end()) {
            continue;
        }
        if (read_has_local_event_signal(record, component.chrom, ref_signal_start, ref_signal_end)) {
            continue;
        }
        if (record->core.qual >= 20) {
            ref_qnames.insert(qname);
        } else {
            low_mapq_ref_qnames.insert(qname);
        }
    }

    evidence.ref_span_reads = static_cast<int32_t>(ref_qnames.size());
    evidence.low_mapq_ref_span_reads = static_cast<int32_t>(low_mapq_ref_qnames.size());
    evidence.ref_span_qnames = sorted_string_set(ref_qnames);
    return evidence;
}

Pipeline::ComponentSignalCache Pipeline::build_component_signal_cache(
    const std::vector<const bam1_t*>& local_records,
    const std::vector<ReadReferenceSpan>& read_spans,
    const std::string& chrom) const {
    ComponentSignalCache cache;
    cache.reads.reserve(local_records.size());
    for (size_t read_idx = 0; read_idx < local_records.size(); ++read_idx) {
        const bam1_t* record = local_records[read_idx];
        if (!record) {
            continue;
        }

        ReadView read(record);
        const uint32_t* cigar = read.cigar();
        const int32_t n_cigar = read.n_cigar();

        CachedLocalSignal cached;
        cached.qname = record_qname(record);
        cached.tid = read.tid();
        cached.mapq = read.mapq();
        cached.is_primary = is_primary_alignment(record);
        cached.has_sa_or_supp =
            read.has_sa_tag() || ((read.flag() & BAM_FSUPPLEMENTARY) != 0);
        cached.split_positions = robust_local_split_insertion_positions(read, chrom);
        if (read_idx < read_spans.size() && read_spans[read_idx].valid) {
            cached.span_valid = true;
            cached.span_start = read_spans[read_idx].start;
            cached.span_end = read_spans[read_idx].end;
        }

        if (cigar && n_cigar > 0) {
            const int first = find_first_non_hard_clip(cigar, n_cigar);
            const int last = find_last_non_hard_clip(cigar, n_cigar);
            int32_t ref_pos = read.pos();
            for (int32_t ci = 0; ci < n_cigar; ++ci) {
                const int op = bam_cigar_op(cigar[ci]);
                const int32_t len = static_cast<int32_t>(bam_cigar_oplen(cigar[ci]));

                if (ci == first && op == BAM_CSOFT_CLIP && len >= kSoftClipSignalMin) {
                    cached.left_clip_pos = ref_pos;
                }
                if (ci == last && op == BAM_CSOFT_CLIP && len >= kSoftClipSignalMin) {
                    cached.right_clip_pos = ref_pos;
                }
                if (op == BAM_CINS && len >= kLongInsertionSignalMin) {
                    cached.indel_positions.push_back(ref_pos);
                }
                if (op == BAM_CINS) {
                    cached.raw_cigar_insertions.emplace_back(ref_pos, len);
                }
                if ((bam_cigar_type(op) & 2) != 0) {
                    ref_pos += len;
                }
            }
        }

        cache.reads.push_back(std::move(cached));
    }
    return cache;
}

EventReadEvidence Pipeline::collect_event_read_evidence_for_bounds_cached(
    const ComponentCall& component,
    const ComponentSignalCache& cache,
    const std::vector<InsertionFragment>& fragments,
    int32_t seed_left,
    int32_t seed_right,
    int32_t bp_left,
    int32_t bp_right) const {
    EventReadEvidence evidence;
    (void)seed_left;
    (void)seed_right;
    constexpr int32_t kAltSignalSlackBp = 25;
    constexpr int32_t kRefSignalSlackBp = 75;
    constexpr int32_t kCandidateClipSlackBp = 25;
    constexpr int32_t kCandidateClipSupplementMaxOffsetBp = 75;
    const int32_t left = std::min(bp_left, bp_right);
    const int32_t right = std::max(bp_left, bp_right);
    evidence.bp_left = left;
    evidence.bp_right = right;

    const auto classify_cached_signal = [&](
        const CachedLocalSignal& cached,
        int32_t window_start,
        int32_t window_end) {
        LocalEventSignal signal;
        if (window_end < window_start) {
            return signal;
        }

        if (cached.left_clip_pos >= window_start && cached.left_clip_pos <= window_end) {
            signal.left_clip = true;
            signal.left_clip_pos = cached.left_clip_pos;
        }
        if (cached.right_clip_pos >= window_start && cached.right_clip_pos <= window_end) {
            signal.right_clip = true;
            signal.right_clip_pos = cached.right_clip_pos;
        }

        const int32_t window_center = window_start + ((window_end - window_start) / 2);
        int32_t best_indel_dist = std::numeric_limits<int32_t>::max();
        for (int32_t pos : cached.indel_positions) {
            if (pos < window_start || pos > window_end) {
                continue;
            }
            signal.indel = true;
            const int32_t dist = std::abs(pos - window_center);
            if (dist < best_indel_dist) {
                best_indel_dist = dist;
                signal.indel_pos = pos;
            }
        }
        for (const auto& raw_insert : cached.raw_cigar_insertions) {
            const int32_t pos = raw_insert.first;
            const int32_t len = raw_insert.second;
            if (pos < window_start || pos > window_end) {
                continue;
            }
            signal.max_raw_cigar_insert_len = std::max(
                signal.max_raw_cigar_insert_len,
                len);
        }

        for (int32_t pos : cached.split_positions) {
            if (pos < window_start || pos > window_end) {
                continue;
            }
            signal.split = true;
            if (signal.split_left_pos < 0 || pos < signal.split_left_pos) {
                signal.split_left_pos = pos;
            }
            if (signal.split_right_pos < 0 || pos > signal.split_right_pos) {
                signal.split_right_pos = pos;
            }
        }
        return signal;
    };

    std::unordered_set<std::string> split_qnames;
    std::unordered_set<std::string> indel_qnames;
    std::unordered_set<std::string> left_clip_qnames;
    std::unordered_set<std::string> right_clip_qnames;
    std::unordered_set<std::string> raw_cigar_insert_qnames;
    std::unordered_set<std::string> nearby_candidate_left_clip_qnames;
    std::unordered_set<std::string> nearby_candidate_right_clip_qnames;

    const int32_t alt_signal_start = std::max(0, evidence.bp_left - kAltSignalSlackBp);
    const int32_t alt_signal_end = std::max(alt_signal_start, evidence.bp_right + kAltSignalSlackBp);
    const int32_t ref_signal_start = std::max(0, evidence.bp_left - kRefSignalSlackBp);
    const int32_t ref_signal_end = std::max(ref_signal_start, evidence.bp_right + kRefSignalSlackBp);
    const int32_t ref_span_start = alt_signal_start;
    const int32_t ref_span_end = alt_signal_end;

    std::vector<int32_t> nearby_component_clip_candidates;
    nearby_component_clip_candidates.reserve(component.breakpoint_candidates.size());
    for (const auto& candidate : component.breakpoint_candidates) {
        if (candidate.pos < 0) {
            continue;
        }
        const int32_t dist_to_hypothesis = (candidate.pos < left)
            ? (left - candidate.pos)
            : ((candidate.pos > right) ? (candidate.pos - right) : 0);
        if (dist_to_hypothesis <= kCandidateClipSupplementMaxOffsetBp) {
            nearby_component_clip_candidates.push_back(candidate.pos);
        }
    }
    std::sort(
        nearby_component_clip_candidates.begin(),
        nearby_component_clip_candidates.end());
    nearby_component_clip_candidates.erase(
        std::unique(
            nearby_component_clip_candidates.begin(),
            nearby_component_clip_candidates.end()),
        nearby_component_clip_candidates.end());

    const bool have_nearby_component_clip_candidates =
        !nearby_component_clip_candidates.empty();
    int32_t candidate_clip_start = 0;
    int32_t candidate_clip_end = 0;
    if (have_nearby_component_clip_candidates) {
        candidate_clip_start = std::max(
            0,
            nearby_component_clip_candidates.front() - kCandidateClipSlackBp);
        candidate_clip_end = std::max(
            candidate_clip_start,
            nearby_component_clip_candidates.back() + kCandidateClipSlackBp);
    }
    const auto matches_nearby_component_clip_candidate = [&](int32_t pos) {
        if (pos < 0) {
            return false;
        }
        for (const int32_t candidate_pos : nearby_component_clip_candidates) {
            if (std::abs(pos - candidate_pos) <= kCandidateClipSlackBp) {
                return true;
            }
        }
        return false;
    };

    for (const auto& cached : cache.reads) {
        if (cached.tid != component.tid || cached.qname.empty()) {
            continue;
        }

        const LocalEventSignal signal = classify_cached_signal(
            cached,
            alt_signal_start,
            alt_signal_end);
        if (signal.split) {
            split_qnames.insert(cached.qname);
        }
        if (signal.indel && cached.mapq == kInsertionCandidateRequiredMapq) {
            indel_qnames.insert(cached.qname);
        }
        if (signal.max_raw_cigar_insert_len > 0) {
            raw_cigar_insert_qnames.insert(cached.qname);
            evidence.max_raw_cigar_insert_len = std::max(
                evidence.max_raw_cigar_insert_len,
                signal.max_raw_cigar_insert_len);
        }
        if (signal.left_clip) {
            left_clip_qnames.insert(cached.qname);
        }
        if (signal.right_clip) {
            right_clip_qnames.insert(cached.qname);
        }

        if (have_nearby_component_clip_candidates) {
            const LocalEventSignal nearby_signal = classify_cached_signal(
                cached,
                candidate_clip_start,
                candidate_clip_end);
            if (nearby_signal.left_clip &&
                matches_nearby_component_clip_candidate(nearby_signal.left_clip_pos)) {
                nearby_candidate_left_clip_qnames.insert(cached.qname);
            }
            if (nearby_signal.right_clip &&
                matches_nearby_component_clip_candidate(nearby_signal.right_clip_pos)) {
                nearby_candidate_right_clip_qnames.insert(cached.qname);
            }
        }
    }

    for (const auto& fragment : fragments) {
        if (fragment.read_id.empty() ||
            fragment.ref_junc_pos < alt_signal_start ||
            fragment.ref_junc_pos > alt_signal_end) {
            continue;
        }
        switch (fragment.source) {
            case InsertionFragmentSource::kSplitSa:
                split_qnames.insert(fragment.read_id);
                break;
            case InsertionFragmentSource::kCigarInsertion:
                indel_qnames.insert(fragment.read_id);
                raw_cigar_insert_qnames.insert(fragment.read_id);
                evidence.max_raw_cigar_insert_len = std::max(
                    evidence.max_raw_cigar_insert_len,
                    fragment.length);
                break;
            case InsertionFragmentSource::kClipRefLeft:
                left_clip_qnames.insert(fragment.read_id);
                break;
            case InsertionFragmentSource::kClipRefRight:
                right_clip_qnames.insert(fragment.read_id);
                break;
            default:
                break;
        }

        if (have_nearby_component_clip_candidates &&
            matches_nearby_component_clip_candidate(fragment.ref_junc_pos)) {
            switch (fragment.source) {
                case InsertionFragmentSource::kClipRefLeft:
                    nearby_candidate_left_clip_qnames.insert(fragment.read_id);
                    break;
                case InsertionFragmentSource::kClipRefRight:
                    nearby_candidate_right_clip_qnames.insert(fragment.read_id);
                    break;
                default:
                    break;
            }
        }
    }

    const bool has_precise_nonclip_support =
        !split_qnames.empty() || !indel_qnames.empty();
    if (has_precise_nonclip_support &&
        !nearby_candidate_left_clip_qnames.empty() &&
        !nearby_candidate_right_clip_qnames.empty()) {
        left_clip_qnames.insert(
            nearby_candidate_left_clip_qnames.begin(),
            nearby_candidate_left_clip_qnames.end());
        right_clip_qnames.insert(
            nearby_candidate_right_clip_qnames.begin(),
            nearby_candidate_right_clip_qnames.end());
    }

    evidence.alt_split_reads = static_cast<int32_t>(split_qnames.size());
    evidence.alt_indel_reads = static_cast<int32_t>(indel_qnames.size());
    evidence.alt_left_clip_reads = static_cast<int32_t>(left_clip_qnames.size());
    evidence.alt_right_clip_reads = static_cast<int32_t>(right_clip_qnames.size());
    evidence.raw_cigar_insert_reads =
        static_cast<int32_t>(raw_cigar_insert_qnames.size());

    std::unordered_set<std::string> alt_qnames = split_qnames;
    alt_qnames.insert(indel_qnames.begin(), indel_qnames.end());
    const bool clip_support_has_partner =
        !alt_qnames.empty() ||
        (!left_clip_qnames.empty() && !right_clip_qnames.empty());
    if (clip_support_has_partner) {
        alt_qnames.insert(left_clip_qnames.begin(), left_clip_qnames.end());
        alt_qnames.insert(right_clip_qnames.begin(), right_clip_qnames.end());
    }

    evidence.alt_struct_reads = static_cast<int32_t>(alt_qnames.size());
    evidence.support_qnames = sorted_string_set(alt_qnames);

    std::unordered_set<std::string> ref_qnames;
    std::unordered_set<std::string> low_mapq_ref_qnames;
    for (const auto& cached : cache.reads) {
        if (!cached.span_valid || cached.tid != component.tid) {
            continue;
        }
        if (cached.span_start > ref_span_start || cached.span_end < ref_span_end) {
            continue;
        }
        if (!cached.is_primary) {
            continue;
        }
        if (cached.qname.empty() || alt_qnames.find(cached.qname) != alt_qnames.end()) {
            continue;
        }
        if (classify_cached_signal(cached, ref_signal_start, ref_signal_end).any()) {
            continue;
        }
        if (cached.mapq >= 20) {
            ref_qnames.insert(cached.qname);
        } else {
            low_mapq_ref_qnames.insert(cached.qname);
        }
    }

    evidence.ref_span_reads = static_cast<int32_t>(ref_qnames.size());
    evidence.low_mapq_ref_span_reads = static_cast<int32_t>(low_mapq_ref_qnames.size());
    evidence.ref_span_qnames = sorted_string_set(ref_qnames);
    return evidence;
}
