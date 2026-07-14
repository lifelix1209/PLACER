namespace {

EvidenceLedgerRow make_joint_evidence_ledger_row(
    const ComponentCall& component,
    const EventReadEvidence& event_evidence,
    const EventConsensus& event_consensus,
    const EventSegmentation& event_segmentation,
    const TEAlignmentEvidence& te_alignment,
    const JointDecisionResult& joint,
    int32_t coverage_bin_start,
    int32_t coverage_bin_end,
    const std::string& retention_reason) {
    EvidenceLedgerRow row;
    row.chrom = component.chrom;
    row.tid = component.tid;
    row.bp_left = event_evidence.bp_left;
    row.bp_right = event_evidence.bp_right;
    row.coverage_left = coverage_bin_start >= 0
        ? coverage_bin_start
        : (component.bin_start >= 0 ? component.bin_start : row.bp_left);
    row.coverage_right = coverage_bin_end > coverage_bin_start
        ? coverage_bin_end - 1
        : (component.bin_end >= 0 ? component.bin_end : row.bp_right);
    row.owner_context_left = coverage_bin_start >= 0
        ? coverage_bin_start
        : (component.bin_start >= 0 ? component.bin_start : row.coverage_left);
    row.owner_context_right = coverage_bin_end > coverage_bin_start
        ? coverage_bin_end - 1
        : (component.bin_end >= 0 ? component.bin_end : row.coverage_right);
    if (row.bp_left >= 0) {
        row.coverage_left = row.coverage_left >= 0
            ? std::min(row.coverage_left, row.bp_left)
            : row.bp_left;
        row.coverage_right = row.coverage_right >= 0
            ? std::max(row.coverage_right, row.bp_left)
            : row.bp_left;
    }
    if (row.bp_right >= 0) {
        row.coverage_left = row.coverage_left >= 0
            ? std::min(row.coverage_left, row.bp_right)
            : row.bp_right;
        row.coverage_right = row.coverage_right >= 0
            ? std::max(row.coverage_right, row.bp_right)
            : row.bp_right;
    }
    row.pos = event_evidence.bp_left >= 0 && event_evidence.bp_right >= 0
        ? event_evidence.bp_left + ((event_evidence.bp_right - event_evidence.bp_left) / 2)
        : component.anchor_pos;
    row.family = te_alignment.best_family.empty() ? "NA" : te_alignment.best_family;
    row.subfamily = te_alignment.best_subfamily.empty() ? "NA" : te_alignment.best_subfamily;
    row.family_alignment_resolved =
        te_alignment.pass &&
        (te_alignment.qc_reason == "PASS_INSERT_TE_ALIGNMENT" ||
         te_alignment.qc_reason == "PASS_INSERT_TE_ALIGNMENT_FAMILY_ONLY") &&
        !te_alignment.best_family.empty();
    row.final_qc = joint.final_qc;
    if (row.final_qc.rfind("PASS_TE", 0) != 0) {
        if (row.bp_left >= 0 && row.bp_right >= 0) {
            row.coverage_left = std::min(row.bp_left, row.bp_right);
            row.coverage_right = std::max(row.bp_left, row.bp_right);
        } else if (row.pos >= 0) {
            row.coverage_left = row.pos;
            row.coverage_right = row.pos;
        }
    }
    row.posterior_qc = joint.posterior_qc;
    row.lfdr_qc = joint.lfdr_qc;
    row.candidate_retention_reason = retention_reason;
    row.alt_struct_reads = event_evidence.alt_struct_reads;
    row.alt_split_reads = event_evidence.alt_split_reads;
    row.alt_indel_reads = event_evidence.alt_indel_reads;
    row.alt_left_clip_reads = event_evidence.alt_left_clip_reads;
    row.alt_right_clip_reads = event_evidence.alt_right_clip_reads;
    row.raw_cigar_insert_reads = event_evidence.raw_cigar_insert_reads;
    row.max_raw_cigar_insert_len = event_evidence.max_raw_cigar_insert_len;
    row.full_context_input_reads = event_consensus.full_context_input_reads;
    row.partial_context_input_reads = event_consensus.partial_context_input_reads;
    row.left_anchor_input_reads = event_consensus.left_anchor_input_reads;
    row.right_anchor_input_reads = event_consensus.right_anchor_input_reads;
    row.input_event_reads = event_consensus.input_event_reads;
    row.event_consensus_len = event_consensus.consensus_len;
    row.left_flank_align_len = event_segmentation.left_flank_align_len;
    row.right_flank_align_len = event_segmentation.right_flank_align_len;
    row.insert_seq = event_segmentation.insert_seq;
    row.ref_span_reads = event_evidence.ref_span_reads;
    row.support_qnames = event_evidence.support_qnames;
    row.best_te_identity = te_alignment.best_identity;
    row.best_te_query_coverage = te_alignment.best_query_coverage;
    row.cross_family_margin = te_alignment.cross_family_margin;
    row.te_structure_path = joint.te_structure_path;
    row.te_structure_log_evidence = joint.te_structure_log_evidence;
    row.nonte_structure_log_evidence = joint.nonte_structure_log_evidence;
    row.artifact_structure_log_evidence = joint.artifact_structure_log_evidence;
    row.te_structure_path_confidence = joint.te_structure_path_confidence;
    row.polyA_posterior = joint.polyA_posterior;
    row.transduction_posterior = joint.transduction_posterior;
    row.te_core_coverage = joint.te_core_coverage;
    row.unexplained_high_complexity_bp = joint.unexplained_high_complexity_bp;
    row.te_posterior = joint.te_posterior;
    row.non_te_posterior = joint.non_te_posterior;
    row.artifact_posterior = joint.artifact_posterior;
    row.lfdr = joint.lfdr;
    row.worst_case_lfdr = joint.worst_case_lfdr;
    row.mechanistic_lower_log_bf_te_vs_artifact =
        joint.mechanistic_lower_log_bf_te_vs_artifact;
    row.mechanistic_lower_log_bf_te_vs_non_te =
        joint.mechanistic_lower_log_bf_te_vs_non_te;
    row.mechanistic_ref_conflict_signal =
        joint.mechanistic_ref_conflict_signal;
    row.mechanistic_ambiguity_width =
        joint.mechanistic_ambiguity_width;
    row.mechanistic_blocks = joint.mechanistic_blocks;
    row.robust_mechanistic_lfdr = joint.robust_mechanistic_lfdr;
    row.robust_mechanistic_worst_case_lfdr =
        joint.robust_mechanistic_worst_case_lfdr;
    row.robust_mechanistic_qc = joint.robust_mechanistic_qc;
    return row;
}

}  // namespace

void Pipeline::process_bin_records(
    std::vector<const bam1_t*>&& bin_records,
    int32_t tid,
    int32_t bin_index,
    int32_t owner_bin_start,
    int32_t owner_bin_end,
    PipelineResult& result) const {
    if (bin_records.empty()) {
        return;
    }
    const auto bin_wall_started = std::chrono::steady_clock::now();

    result.processed_bins += 1;

    int32_t bin_start = std::numeric_limits<int32_t>::max();
    int32_t bin_end = std::numeric_limits<int32_t>::min();
    for (const bam1_t* record : bin_records) {
        if (!record) {
            continue;
        }
        ReadView read(record);
        bin_start = std::min(bin_start, read.pos());
        bin_end = std::max(bin_end, compute_ref_end(read));
    }
    if (bin_start == std::numeric_limits<int32_t>::max()) {
        bin_start = std::max(0, bin_index);
    }
    if (bin_end == std::numeric_limits<int32_t>::min()) {
        bin_end = bin_start;
    }
    const std::string chrom = bam_reader_->chromosome_name(tid);
    const auto make_summary_ledger_row = [owner_bin_start, owner_bin_end](
                                            const ComponentCall& component,
                                            const HypothesisSummary& summary,
                                            const std::string& retention_reason) {
        EvidenceLedgerRow row;
        row.chrom = component.chrom;
        row.tid = component.tid;
        row.bp_left = summary.bp_left;
        row.bp_right = summary.bp_right;
        row.pos = summary.bp_left >= 0 && summary.bp_right >= 0
            ? summary.bp_left + ((summary.bp_right - summary.bp_left) / 2)
            : component.anchor_pos;
        row.owner_context_left = owner_bin_start < owner_bin_end
            ? owner_bin_start
            : component.bin_start;
        row.owner_context_right = owner_bin_start < owner_bin_end
            ? owner_bin_end - 1
            : component.bin_end;
        if (row.bp_left >= 0 && row.bp_right >= 0) {
            row.coverage_left = std::min(row.bp_left, row.bp_right);
            row.coverage_right = std::max(row.bp_left, row.bp_right);
        } else if (row.pos >= 0) {
            row.coverage_left = row.pos;
            row.coverage_right = row.pos;
        }
        row.final_qc = "NOT_EVALUATED_PRE_EXPENSIVE_STAGE";
        row.posterior_qc = "POSTERIOR_NOT_EVALUATED";
        row.lfdr_qc = "LFDR_NOT_EVALUATED";
        row.candidate_retention_reason = retention_reason;
        row.alt_struct_reads = summary.alt_struct_reads;
        row.alt_split_reads = summary.alt_split_reads;
        row.alt_indel_reads = summary.alt_indel_reads;
        row.alt_left_clip_reads = summary.alt_left_clip_reads;
        row.alt_right_clip_reads = summary.alt_right_clip_reads;
        row.raw_cigar_insert_reads = summary.raw_cigar_insert_reads;
        row.max_raw_cigar_insert_len = summary.max_raw_cigar_insert_len;
        row.ref_span_reads = summary.ref_span_reads;
        row.support_qnames = summary.support_qnames;
        return row;
    };

    const auto component_build_started = std::chrono::steady_clock::now();
    auto components = component_module_.build(bin_records, chrom, tid, bin_start, bin_end);
    if (owner_bin_start < owner_bin_end) {
        components.erase(
            std::remove_if(
                components.begin(),
                components.end(),
                [owner_bin_start, owner_bin_end](const ComponentCall& component) {
                    return component.anchor_pos < owner_bin_start ||
                        component.anchor_pos >= owner_bin_end;
                }),
            components.end());
    }
    const double component_build_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - component_build_started).count();
    result.built_components += static_cast<int64_t>(components.size());
    std::vector<std::pair<int32_t, int32_t>> component_seed_bounds;
    component_seed_bounds.reserve(components.size());
    std::vector<LocalIntervalRequest> requests;
    requests.reserve(components.size());
    for (size_t ci = 0; ci < components.size(); ++ci) {
        const auto& component = components[ci];
        const auto seed_bounds = infer_component_breakpoint_bounds(component);
        component_seed_bounds.push_back(seed_bounds);
        requests.push_back(LocalIntervalRequest{
            component.chrom,
            std::max(0, std::min(seed_bounds.first, seed_bounds.second) - kLocalEventFetchSlackBp),
            std::max(
                std::max(0, std::min(seed_bounds.first, seed_bounds.second) - kLocalEventFetchSlackBp) + 1,
                std::max(seed_bounds.first, seed_bounds.second) + kLocalEventFetchSlackBp),
            ci,
        });
    }
    LocalIntervalReuseStats reuse_stats;
    reuse_stats.request_count = requests.size();
    const std::vector<CanonicalLocalInterval> intervals =
        build_canonical_local_intervals(requests, 128);
    reuse_stats.canonical_interval_count = intervals.size();
    std::vector<LocalIntervalCacheEntry> interval_cache;
    interval_cache.reserve(intervals.size());
    const auto local_fetch_started = std::chrono::steady_clock::now();
    for (const auto& interval : intervals) {
        LocalIntervalCacheEntry entry;
        entry.interval = interval;
        const bool fetch_ok = bam_reader_->fetch(
            interval.chrom,
            interval.start,
            interval.end,
            [&](BamRecordPtr&& record) {
                if (!record) {
                    return true;
                }
                entry.owned_records.push_back(std::move(record));
                const bam1_t* ptr = entry.owned_records.back().get();
                entry.records.push_back(ptr);
                ReadView view(ptr);
                ReadReferenceSpan span;
                span.valid = true;
                span.tid = view.tid();
                span.start = view.pos();
                span.end = compute_ref_end(view);
                entry.read_spans.push_back(span);
                return true;
            });
        if (!fetch_ok) {
            entry.records.clear();
            entry.read_spans.clear();
        }
        interval_cache.push_back(std::move(entry));
    }
    const double local_fetch_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - local_fetch_started).count();

    struct BinStageStats {
        int64_t breakpoint_hypotheses_total = 0;
        int64_t breakpoint_hypotheses_after_collapse = 0;
        int64_t breakpoint_hypotheses_after_gate = 0;
        int64_t validator_candidates_total = 0;
        int64_t validator_feasible_total = 0;
        int64_t shortlist_empty = 0;
        int64_t shortlist_primary_only = 0;
        int64_t shortlist_with_challenger = 0;
        int64_t expensive_stage_hypotheses_total = 0;
        int64_t event_consensus_calls = 0;
        int64_t event_consensus_rejected = 0;
        int64_t pre_segmentation_gate_rejected = 0;
        int64_t pre_segmentation_gate_no_bilateral_anchor = 0;
        int64_t pre_segmentation_gate_no_precise_or_full_context = 0;
        int64_t event_segmentation_calls = 0;
        int64_t event_segmentation_rejected = 0;
        int64_t te_alignment_calls = 0;
        int64_t te_alignment_rejected = 0;
        int64_t genotype_calls = 0;
        int64_t final_calls = 0;
        double component_build_seconds = 0.0;
        double local_fetch_seconds = 0.0;
        double projection_seconds = 0.0;
        double local_component_refresh_seconds = 0.0;
        double fragment_extract_seconds = 0.0;
        double signal_cache_build_seconds = 0.0;
        double hypothesis_summary_seconds = 0.0;
        double validator_input_count_seconds = 0.0;
        double event_evidence_seconds = 0.0;
        ExpensiveStageTiming timing;
    } bin_stats;
    bin_stats.component_build_seconds = component_build_seconds;
    bin_stats.local_fetch_seconds = local_fetch_seconds;

    if (config_.log_stage_bins) {
        std::ostringstream oss;
        oss << "[Pipeline][bin][start]"
            << " chrom=" << chrom
            << " tid=" << tid
            << " start=" << bin_start
            << " end=" << bin_end
            << " records=" << bin_records.size()
            << " components=" << components.size();
        emit_pipeline_log_line(oss.str());
    }

    const auto log_component = [&](
        const char* outcome,
        const ComponentCall& component,
        const std::vector<InsertionFragment>& fragments,
        const GenotypeCall* genotype,
        const FinalCall* final_call) {
        if (!config_.log_stage_components) {
            return;
        }

        std::ostringstream oss;
        oss << "[Pipeline][component]"
            << " outcome=" << outcome
            << " chrom=" << component.chrom
            << " anchor=" << component.anchor_pos
            << " bin=" << component.bin_start << "-" << component.bin_end
            << " reads=" << component.read_indices.size()
            << " softclip_reads=" << component.soft_clip_read_indices.size()
            << " split_reads=" << component.split_sa_read_indices.size()
            << " indel_reads=" << component.insertion_read_indices.size()
            << " breakpoints=" << component.breakpoint_candidates.size()
            << " fragments=" << fragments.size();

        if (genotype != nullptr) {
            oss << " gt=" << display_or_na(genotype->genotype)
                << " af=" << genotype->af
                << " gq=" << genotype->gq;
        }

        if (final_call != nullptr) {
            oss << " final_te=" << display_or_na(final_call->te_name)
                << " final_support=" << final_call->support_reads
                << " final_qc=" << display_or_na(final_call->final_qc)
                << " te_qc=" << display_or_na(final_call->te_qc);
        }

        emit_pipeline_log_line(oss.str());
    };

    std::unordered_map<std::string, EventSegmentation> event_segmentation_cache;
    std::unordered_map<std::string, TEAlignmentEvidence> te_alignment_cache;

    for (size_t ci = 0; ci < components.size(); ++ci) {
        const auto& component = components[ci];
        const std::vector<InsertionFragment> empty_fragments;
        const auto seed_bounds = component_seed_bounds[ci];
        const auto projection_started = std::chrono::steady_clock::now();
        const LocalIntervalProjection projection =
            project_cached_interval_reads(requests[ci], interval_cache);
        bin_stats.projection_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - projection_started).count();
        if (projection.records.empty()) {
            log_component(
                "LOCAL_EVENT_RECOLLECTION_FAILED",
                component,
                empty_fragments,
                nullptr,
                nullptr);
            continue;
        }

        const auto local_component_refresh_started = std::chrono::steady_clock::now();
        const ComponentCall local_component = build_local_fragment_component(
            config_,
            component,
            projection.records);
        bin_stats.local_component_refresh_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - local_component_refresh_started).count();
        const auto fragment_extract_started = std::chrono::steady_clock::now();
        const std::vector<InsertionFragment> fragments = ins_fragment_module_.extract(
            local_component,
            projection.records);
        bin_stats.fragment_extract_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - fragment_extract_started).count();
        const auto signal_cache_build_started = std::chrono::steady_clock::now();
        const ComponentSignalCache signal_cache = build_component_signal_cache(
            projection.records,
            projection.read_spans,
            component.chrom);
        bin_stats.signal_cache_build_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - signal_cache_build_started).count();
        const auto all_breakpoint_hypotheses = enumerate_breakpoint_hypotheses(
            component,
            projection.records,
            fragments,
            seed_bounds.first,
            seed_bounds.second,
            0);

        if (config_.log_stage_components) {
            for (const auto& hypothesis : all_breakpoint_hypotheses) {
                std::ostringstream bp_oss;
                bp_oss << "[Pipeline][breakpoint_hypothesis]"
                       << " chrom=" << component.chrom
                       << " anchor=" << component.anchor_pos
                       << " bp=" << hypothesis.left << "-" << hypothesis.right
                       << " center=" << hypothesis.center
                       << " support=" << hypothesis.support
                       << " priority=" << hypothesis.priority
                       << " score="
                       << (hypothesis.support *
                           breakpoint_hypothesis_support_weight(hypothesis.priority));
                emit_pipeline_log_line(bp_oss.str());
            }
        }

        const auto breakpoint_hypotheses = select_diverse_breakpoint_hypotheses(
            all_breakpoint_hypotheses,
            0,
            component.anchor_pos);
        bin_stats.breakpoint_hypotheses_total += static_cast<int64_t>(breakpoint_hypotheses.size());
        std::vector<HypothesisSummary> hypothesis_summaries;
        hypothesis_summaries.reserve(breakpoint_hypotheses.size());
        const auto hypothesis_summary_started = std::chrono::steady_clock::now();
        for (size_t hi = 0; hi < breakpoint_hypotheses.size(); ++hi) {
            const auto& hypothesis = breakpoint_hypotheses[hi];
            hypothesis_summaries.push_back(build_hypothesis_summary(
                component,
                signal_cache,
                fragments,
                hi,
                hypothesis.left,
                hypothesis.right,
                hypothesis.support,
                hypothesis.priority));
        }
        bin_stats.hypothesis_summary_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - hypothesis_summary_started).count();
        const auto collapsed_summaries = collapse_hypothesis_summaries(hypothesis_summaries);
        bin_stats.breakpoint_hypotheses_after_collapse +=
            static_cast<int64_t>(collapsed_summaries.size());
        std::vector<HypothesisSummary> coarse_survivors;
        coarse_survivors.reserve(collapsed_summaries.size());
        for (size_t si = 0; si < collapsed_summaries.size(); ++si) {
            const bool keep_for_expensive_stage =
                should_keep_hypothesis_for_expensive_stage(
                    collapsed_summaries[si],
                    si == 0);
            if (keep_for_expensive_stage) {
                coarse_survivors.push_back(collapsed_summaries[si]);
            } else if (should_record_hypothesis_in_evidence_ledger(
                           collapsed_summaries[si])) {
                result.evidence_ledger.push_back(make_summary_ledger_row(
                    component,
                    collapsed_summaries[si],
                    "LEDGER_ONLY_PRE_EXPENSIVE_STAGE"));
            }
        }
        bin_stats.breakpoint_hypotheses_after_gate +=
            static_cast<int64_t>(coarse_survivors.size());

        std::vector<HypothesisValidatorEvidence> validator_candidates;
        validator_candidates.reserve(coarse_survivors.size());
        for (const auto& summary : coarse_survivors) {
            const auto validator_input_started = std::chrono::steady_clock::now();
            const ConsensusInputCounts inputs = collect_event_consensus_input_counts(
                projection.records,
                fragments,
                summary.event_evidence);
            bin_stats.validator_input_count_seconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - validator_input_started).count();
            validator_candidates.push_back(
                collect_hypothesis_validator_evidence(
                    summary,
                    inputs,
                    component.anchor_pos));
        }
        const auto shortlist = build_expensive_stage_shortlist(validator_candidates);
        bin_stats.validator_candidates_total += static_cast<int64_t>(validator_candidates.size());
        for (const auto& candidate : validator_candidates) {
            if (candidate.feasible_for_expensive_stage) {
                ++bin_stats.validator_feasible_total;
            }
        }
        if (shortlist.empty()) {
            ++bin_stats.shortlist_empty;
        } else if (shortlist.size() == 1) {
            ++bin_stats.shortlist_primary_only;
        } else {
            ++bin_stats.shortlist_with_challenger;
        }
        bin_stats.expensive_stage_hypotheses_total += static_cast<int64_t>(shortlist.size());

        if (config_.log_stage_components) {
            std::ostringstream shortlist_oss;
            shortlist_oss << "[Pipeline][shortlist]"
                          << " chrom=" << component.chrom
                          << " anchor=" << component.anchor_pos
                          << " collapsed_hypotheses=" << collapsed_summaries.size()
                          << " coarse_survivors=" << coarse_survivors.size()
                          << " validator_feasible=" << std::count_if(
                                 validator_candidates.begin(),
                                 validator_candidates.end(),
                                 [](const auto& row) { return row.feasible_for_expensive_stage; })
                          << " shortlist_size=" << shortlist.size();
            if (!shortlist.empty()) {
                shortlist_oss << " primary_bp="
                              << shortlist.front().validator.summary.bp_left
                              << "-" << shortlist.front().validator.summary.bp_right;
            }
            if (shortlist.size() > 1) {
                shortlist_oss << " challenger_bp="
                              << shortlist[1].validator.summary.bp_left
                              << "-" << shortlist[1].validator.summary.bp_right;
            }
            emit_pipeline_log_line(shortlist_oss.str());

            for (const auto& candidate : validator_candidates) {
                std::ostringstream candidate_oss;
                candidate_oss << "[Pipeline][shortlist_candidate]"
                              << " chrom=" << component.chrom
                              << " anchor=" << component.anchor_pos
                              << " bp=" << candidate.summary.bp_left
                              << "-" << candidate.summary.bp_right
                              << " precise=" << candidate.precise_support
                              << " alt_struct=" << candidate.summary.alt_struct_reads
                              << " ref_span=" << candidate.summary.ref_span_reads
                              << " full_context=" << candidate.full_context_input_reads
                              << " partial_context=" << candidate.partial_context_input_reads
                              << " input_event_reads=" << candidate.input_event_reads
                              << " left_anchor=" << candidate.left_anchor_input_reads
                              << " right_anchor=" << candidate.right_anchor_input_reads
                              << " feasible=" << bool_as_int(candidate.feasible_for_expensive_stage)
                              << " qc=" << candidate.qc_reason;
                emit_pipeline_log_line(candidate_oss.str());
            }
        }

        bool have_best = false;
        double best_total = -1e9;
        JointDecisionResult best_joint;
        EventReadEvidence best_event_evidence;
        EventConsensus best_event_consensus;
        EventSegmentation best_event_segmentation;
        TEAlignmentEvidence best_te_alignment;
        BoundaryEvidence best_boundary_evidence;
        ClipInsertConcordanceEvidence best_clip_insert_concordance;
        GenotypeCall best_genotype;

        struct ComponentFinalCallEvaluation {
            FinalCall call;
            GenotypeCall genotype;
            double score = 0.0;
            int32_t anchor_pos = -1;
            int32_t anchor_support = 0;
            int32_t anchor_ref_span_reads = 0;
            int32_t anchor_priority = std::numeric_limits<int32_t>::max();
            double anchor_hypothesis_score = 0.0;
            int32_t retether_direction = 0;
            bool evidence_te = false;
            bool one_sided_segmentation = false;
        };
        std::vector<ComponentFinalCallEvaluation> final_call_candidates;
        final_call_candidates.reserve(shortlist.size());
        std::vector<ComponentFinalCallCandidate> component_anchor_candidates;
        component_anchor_candidates.reserve(all_breakpoint_hypotheses.size());
        for (const auto& hypothesis : all_breakpoint_hypotheses) {
            if (!hypothesis.valid || hypothesis.center < 0) {
                continue;
            }
            ComponentFinalCallCandidate anchor;
            anchor.pos = hypothesis.center;
            anchor.emit_te = false;
            anchor.anchor_support = hypothesis.support;
            anchor.anchor_priority = hypothesis.priority;
            anchor.anchor_hypothesis_score =
                static_cast<double>(hypothesis.support) *
                static_cast<double>(breakpoint_hypothesis_support_weight(hypothesis.priority));
            anchor.component_anchor_pos = component.anchor_pos;
            component_anchor_candidates.push_back(anchor);
        }

        struct EvaluatedShortlistCandidate {
            ShortlistedHypothesis shortlisted;
            EventReadEvidence event_evidence;
            EventConsensus event_consensus;
            GenotypeCall genotype;
            EventExistenceEvidence existence;
            EventSegmentation event_segmentation;
            EventSegmentationEvidence seg_evidence;
            TEAlignmentEvidence te_alignment;
        };
        std::vector<EvaluatedShortlistCandidate> evaluated_candidates;
        evaluated_candidates.reserve(shortlist.size());
        std::vector<std::string> te_batch_insert_sequences;
        std::vector<std::vector<size_t>> te_batch_candidate_index_groups;
        std::unordered_map<std::string, size_t> te_batch_query_index_by_sequence;
        te_batch_insert_sequences.reserve(shortlist.size());
        te_batch_candidate_index_groups.reserve(shortlist.size());
        te_batch_query_index_by_sequence.reserve(shortlist.size());

        for (const auto& shortlisted : shortlist) {
            const EventReadEvidence& event_evidence =
                shortlisted.validator.event_evidence;
            const auto consensus_started = std::chrono::steady_clock::now();
            EventConsensus event_consensus = build_event_consensus(
                component,
                projection.records,
                fragments,
                event_evidence);
            bin_stats.timing.event_consensus_seconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - consensus_started).count();
            result.event_consensus_calls += 1;
            bin_stats.event_consensus_calls += 1;
            if (!event_consensus.qc_pass) {
                bin_stats.event_consensus_rejected += 1;
            }

            const GenotypeCall genotype = genotype_call(component, event_evidence);
            result.genotype_calls += 1;
            bin_stats.genotype_calls += 1;

            EventGenotypeInput existence_input;
            existence_input.alt_struct_reads = event_evidence.alt_struct_reads;
            existence_input.alt_split_reads = event_evidence.alt_split_reads;
            existence_input.alt_indel_reads = event_evidence.alt_indel_reads;
            existence_input.alt_left_clip_reads = event_evidence.alt_left_clip_reads;
            existence_input.alt_right_clip_reads = event_evidence.alt_right_clip_reads;
            existence_input.ref_span_reads = event_evidence.ref_span_reads;
            existence_input.min_depth = config_.genotype_min_depth;
            existence_input.min_gq = 20;
            existence_input.error_rate = config_.genotype_error_rate;
            existence_input.alt_observed_lengths =
                collect_alt_observed_lengths(component, event_evidence);
            existence_input.event_length =
                infer_event_length_from_alt_support(existence_input.alt_observed_lengths);
            const EventExistenceEvidence existence = build_event_existence_evidence(
                existence_input);

            EventSegmentation event_segmentation;
            const auto run_segmentation_for_consensus = [&](
                const EventConsensus& consensus_to_segment) {
                EventSegmentation segmentation;
                const std::string preseg_qc = pre_segmentation_gate_reason(
                    event_evidence,
                    consensus_to_segment);
                if (!preseg_qc.empty()) {
                    segmentation.qc_reason = preseg_qc;
                    ++bin_stats.pre_segmentation_gate_rejected;
                    if (preseg_qc == "PRESEG_NO_BILATERAL_ANCHOR") {
                        ++bin_stats.pre_segmentation_gate_no_bilateral_anchor;
                    } else if (preseg_qc == "PRESEG_NO_PRECISE_OR_FULL_CONTEXT") {
                        ++bin_stats.pre_segmentation_gate_no_precise_or_full_context;
                    }
                    return segmentation;
                }
                const std::string segmentation_key = event_segmentation_cache_key(
                    component,
                    event_evidence,
                    consensus_to_segment);
                const bool segmentation_cache_hit =
                    event_segmentation_cache.find(segmentation_key) != event_segmentation_cache.end();
                if (segmentation_cache_hit) {
                    segmentation = event_segmentation_cache.at(segmentation_key);
                    result.performance.segmentation_cache_hits += 1;
                } else {
                    const auto segmentation_started = std::chrono::steady_clock::now();
                    segmentation = segment_event_consensus_cached(
                        component,
                        event_evidence,
                        consensus_to_segment,
                        event_segmentation_cache);
                    bin_stats.timing.event_segmentation_seconds += std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - segmentation_started).count();
                    bin_stats.event_segmentation_calls += 1;
                    if (!segmentation.pass) {
                        bin_stats.event_segmentation_rejected += 1;
                    }
                    const EventSegmentationSearchStats segmentation_stats =
                        last_segmentation_stats();
                    result.performance.segmentation_edit_distance_calls +=
                        segmentation_stats.edit_distance_calls;
                    result.performance.segmentation_edit_distance_cache_hits +=
                        segmentation_stats.edit_distance_cache_hits;
                    result.performance.segmentation_edit_distance_cache_misses +=
                        segmentation_stats.edit_distance_cache_misses;
                    result.performance.segmentation_seed_bins_total +=
                        segmentation_stats.seed_bins_total;
                    result.performance.segmentation_paired_searches +=
                        segmentation_stats.paired_searches;
                    result.performance.segmentation_one_sided_searches +=
                        segmentation_stats.one_sided_searches;
                    result.performance.segmentation_endpoint_slack_searches +=
                        segmentation_stats.endpoint_slack_searches;
                    result.performance.segmentation_reverse_complement_retries +=
                        segmentation_stats.reverse_complement_retries;
                    result.performance.segmentation_cache_misses +=
                        segmentation_stats.segmentation_cache_misses;
                }
                return segmentation;
            };
            event_segmentation = run_segmentation_for_consensus(event_consensus);
            const bool can_retry_partial_context =
                !event_segmentation.pass &&
                event_consensus.used_full_context &&
                event_consensus.partial_context_input_reads >=
                    std::max(4, config_.event_consensus_poa_min_reads * 2) &&
                event_consensus.left_anchor_input_reads > 0 &&
                event_consensus.right_anchor_input_reads > 0 &&
                (event_segmentation.qc_reason == "NO_LEFT_FLANK_MATCH" ||
                 event_segmentation.qc_reason == "NO_RIGHT_FLANK_MATCH" ||
                 event_segmentation.qc_reason == "EMPTY_EVENT_INSERT_SEGMENT");
            if (can_retry_partial_context) {
                const auto partial_consensus_started = std::chrono::steady_clock::now();
                EventConsensus partial_consensus = build_event_consensus(
                    component,
                    projection.records,
                    fragments,
                    event_evidence,
                    ConsensusContextMode::kPartialOnly);
                bin_stats.timing.event_consensus_seconds += std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - partial_consensus_started).count();
                result.event_consensus_calls += 1;
                bin_stats.event_consensus_calls += 1;
                if (!partial_consensus.qc_pass) {
                    bin_stats.event_consensus_rejected += 1;
                } else {
                    EventSegmentation partial_segmentation =
                        run_segmentation_for_consensus(partial_consensus);
                    if (partial_segmentation.pass) {
                        event_consensus = std::move(partial_consensus);
                        event_segmentation = std::move(partial_segmentation);
                    }
                }
            }
            const EventSegmentationEvidence seg_evidence = analyze_event_segmentation(
                event_consensus,
                event_segmentation);

            EvaluatedShortlistCandidate evaluated;
            evaluated.shortlisted = shortlisted;
            evaluated.event_evidence = event_evidence;
            evaluated.event_consensus = std::move(event_consensus);
            evaluated.genotype = genotype;
            evaluated.existence = existence;
            evaluated.event_segmentation = std::move(event_segmentation);
            evaluated.seg_evidence = seg_evidence;
            if (seg_evidence.has_insert_seq) {
                const bool cacheable_te_alignment =
                    evaluated.event_segmentation.pass &&
                    !evaluated.event_segmentation.insert_seq.empty();
                if (cacheable_te_alignment) {
                    const auto cached_te_alignment =
                        te_alignment_cache.find(evaluated.event_segmentation.insert_seq);
                    if (cached_te_alignment != te_alignment_cache.end()) {
                        evaluated.te_alignment = cached_te_alignment->second;
                    } else {
                        auto [query_it, inserted] =
                            te_batch_query_index_by_sequence.emplace(
                                evaluated.event_segmentation.insert_seq,
                                te_batch_insert_sequences.size());
                        if (inserted) {
                            te_batch_insert_sequences.push_back(
                                evaluated.event_segmentation.insert_seq);
                            te_batch_candidate_index_groups.push_back({});
                        }
                        te_batch_candidate_index_groups[query_it->second].push_back(
                            evaluated_candidates.size());
                    }
                } else {
                    const auto te_started = std::chrono::steady_clock::now();
                    evaluated.te_alignment =
                        align_insert_seq_to_te(evaluated.event_segmentation);
                    bin_stats.timing.te_alignment_seconds += std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - te_started).count();
                    bin_stats.te_alignment_calls += 1;
                    if (!evaluated.te_alignment.pass) {
                        bin_stats.te_alignment_rejected += 1;
                    }
                }
            }
            evaluated_candidates.push_back(std::move(evaluated));
        }

        if (!te_batch_insert_sequences.empty()) {
            const auto te_started = std::chrono::steady_clock::now();
            const std::vector<TEAlignmentEvidence> batch_alignments =
                te_classifier_module_.align_insert_sequences(te_batch_insert_sequences);
            const TEAlignmentBatchStats te_batch_stats =
                te_classifier_module_.last_alignment_batch_stats();
            result.performance.te_blast_batches += te_batch_stats.blast_batches;
            result.performance.te_blast_queries += te_batch_stats.blast_queries;
            result.performance.te_blast_cache_hits += te_batch_stats.cache_hits;
            result.performance.te_blast_cache_misses += te_batch_stats.cache_misses;
            result.performance.te_blast_deduplicated_queries +=
                te_batch_stats.deduplicated_queries;
            bin_stats.timing.te_alignment_seconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - te_started).count();
            for (size_t i = 0; i < batch_alignments.size() &&
                               i < te_batch_candidate_index_groups.size();
                 ++i) {
                const TEAlignmentEvidence& alignment = batch_alignments[i];
                const std::string& insert_seq = te_batch_insert_sequences[i];
                te_alignment_cache.emplace(insert_seq, alignment);
                for (const size_t candidate_index : te_batch_candidate_index_groups[i]) {
                    evaluated_candidates[candidate_index].te_alignment = alignment;
                }
                bin_stats.te_alignment_calls += 1;
                if (!alignment.pass) {
                    bin_stats.te_alignment_rejected += 1;
                }
            }
        }

        for (const auto& evaluated : evaluated_candidates) {
            const auto& shortlisted = evaluated.shortlisted;
            const EventReadEvidence& event_evidence = evaluated.event_evidence;
            const EventConsensus& event_consensus = evaluated.event_consensus;
            const GenotypeCall& genotype = evaluated.genotype;
            const EventExistenceEvidence& existence = evaluated.existence;
            const EventSegmentation& event_segmentation = evaluated.event_segmentation;
            const EventSegmentationEvidence& seg_evidence = evaluated.seg_evidence;
            const TEAlignmentEvidence& te_alignment = evaluated.te_alignment;

            const int32_t breakpoint_envelope_width = std::max(
                1,
                event_evidence.bp_right - event_evidence.bp_left);
            BoundaryEvidence boundary_evidence;
            if (event_segmentation.pass) {
                boundary_evidence = analyze_boundary(
                    event_segmentation,
                    breakpoint_envelope_width);
            }

            const ClipInsertConcordanceEvidence clip_insert_concordance =
                analyze_clip_insert_concordance(
                    event_evidence,
                    event_segmentation,
                    fragments);

            const JointDecisionResult joint = evaluate_joint_hypotheses(
                existence,
                seg_evidence,
                te_alignment,
                boundary_evidence,
                &clip_insert_concordance);
            const double candidate_total = joint.best.hard_veto ? -1e9 : joint.best.total;
            result.evidence_ledger.push_back(make_joint_evidence_ledger_row(
                component,
                event_evidence,
                event_consensus,
                event_segmentation,
                te_alignment,
                joint,
                owner_bin_start,
                owner_bin_end,
                shortlisted.is_primary
                    ? "EXPENSIVE_STAGE_PRIMARY"
                    : "EXPENSIVE_STAGE_CHALLENGER"));

            if (config_.log_stage_components) {
                std::ostringstream candidate_joint_oss;
                const int32_t precise_reads =
                    std::max(0, event_evidence.alt_split_reads) +
                    std::max(0, event_evidence.alt_indel_reads);
                const int32_t bilateral_clip_reads = std::min(
                    std::max(0, event_evidence.alt_left_clip_reads),
                    std::max(0, event_evidence.alt_right_clip_reads));
                candidate_joint_oss << "[Pipeline][shortlist_result]"
                                    << " chrom=" << component.chrom
                                    << " anchor=" << component.anchor_pos
                                    << " bp=" << event_evidence.bp_left
                                    << "-" << event_evidence.bp_right
                                    << " primary=" << bool_as_int(shortlisted.is_primary)
                                    << " alt_reads=" << event_evidence.alt_struct_reads
                                    << " alt_split=" << event_evidence.alt_split_reads
                                    << " alt_indel=" << event_evidence.alt_indel_reads
                                    << " alt_left_clip=" << event_evidence.alt_left_clip_reads
                                    << " alt_right_clip=" << event_evidence.alt_right_clip_reads
                                    << " precise_reads=" << precise_reads
                                    << " bilateral_clip_reads=" << bilateral_clip_reads
                                    << " ref_reads=" << event_evidence.ref_span_reads
                                    << " consensus_qc=" << display_or_na(event_consensus.qc_reason)
                                    << " consensus_len=" << event_consensus.consensus_len
                                    << " seg_qc=" << display_or_na(event_segmentation.qc_reason)
                                    << " seg_insert_len=" << event_segmentation.insert_seq.size()
                                    << " seg_score=" << seg_evidence.score
                                    << " te_qc=" << display_or_na(te_alignment.qc_reason)
                                    << " te_pass=" << bool_as_int(te_alignment.pass)
                                    << " te_family=" << display_or_na(te_alignment.best_family)
                                    << " te_subfamily=" << display_or_na(te_alignment.best_subfamily)
                                    << " te_identity=" << te_alignment.best_identity
                                    << " te_cov=" << te_alignment.best_query_coverage
                                    << " te_margin=" << te_alignment.cross_family_margin
                                    << " te_model_label="
                                    << display_or_na(te_alignment.sequence_model_label)
                                    << " te_model_score="
                                    << te_alignment.sequence_model_score
                                    << " te_model_gc="
                                    << te_alignment.sequence_model_gc
                                    << " te_model_entropy="
                                    << te_alignment.sequence_model_entropy
                                    << " te_model_tandem="
                                    << te_alignment.sequence_model_tandem_fraction
                                    << " te_model_low_complexity="
                                    << te_alignment.sequence_model_low_complexity_fraction
                                    << " te_model_jsd_k5="
                                    << te_alignment.sequence_model_jsd_k5
                                    << " te_model_jsd_k6="
                                    << te_alignment.sequence_model_jsd_k6
                                    << " te_model_k9_containment="
                                    << te_alignment.sequence_model_k9_containment
                                    << " te_annotation_confidence="
                                    << display_or_na(te_alignment.annotation_confidence)
                                    << " te_annotation_residual="
                                    << te_alignment.annotation_residual_fraction
                                    << " te_annotation_masked="
                                    << te_alignment.annotation_masked_fraction
                                    << " te_explain_status="
                                    << te_annotation_status_name(
                                        te_alignment.te_sequence_explanation.status)
                                    << " te_explain_residual="
                                    << serialize_te_path_residual(
                                        te_alignment.te_sequence_explanation.residual)
                                    << " te_explain_path="
                                    << serialize_te_sequence_path(
                                        te_alignment.te_sequence_explanation)
                                    << " boundary_qc=" << display_or_na(boundary_evidence.qc)
                                    << " boundary_defined=" << bool_as_int(boundary_evidence.geometry_defined)
                                    << " boundary_canonical=" << bool_as_int(boundary_evidence.canonical_pass)
                                    << " boundary_consistent=" << bool_as_int(boundary_evidence.evidence_consistent)
                                    << " boundary_type=" << display_or_na(boundary_evidence.boundary_type)
                                    << " boundary_len=" << boundary_evidence.boundary_len
                                    << " clip_insert_qc=" << display_or_na(clip_insert_concordance.qc)
                                    << " clip_insert_full=" << clip_insert_concordance.full_insert_reads
                                    << " clip_insert_left=" << clip_insert_concordance.left_clip_reads
                                    << " clip_insert_right=" << clip_insert_concordance.right_clip_reads
                                    << " clip_insert_max_left_identity="
                                    << clip_insert_concordance.max_left_identity
                                    << " clip_insert_max_right_identity="
                                    << clip_insert_concordance.max_right_identity
                                    << " best_kind=" << final_hypothesis_kind_name(joint.best.kind)
                                    << " best_total=" << joint.best.total
                                    << " best_hard_veto=" << bool_as_int(joint.best.hard_veto)
                                    << " runner_kind="
                                    << final_hypothesis_kind_name(joint.runner_up.kind)
                                    << " runner_total=" << joint.runner_up.total
                                    << " runner_hard_veto=" << bool_as_int(joint.runner_up.hard_veto)
                                    << " final_qc=" << display_or_na(joint.final_qc)
                                    << " emit_te=" << bool_as_int(joint.emit_te_call)
                                    << " emit_structural="
                                    << bool_as_int(joint.emit_structural_event_call)
                                    << " emit_unknown=" << bool_as_int(joint.emit_unknown_te);
                emit_pipeline_log_line(candidate_joint_oss.str());
            }

            if (!have_best || candidate_total > best_total) {
                have_best = true;
                best_total = candidate_total;
                best_joint = joint;
                best_event_evidence = event_evidence;
                best_event_consensus = event_consensus;
                best_event_segmentation = event_segmentation;
                best_te_alignment = te_alignment;
                best_boundary_evidence = boundary_evidence;
                best_clip_insert_concordance = clip_insert_concordance;
                best_genotype = genotype;
            }

            if (joint.emit_te_call || joint.emit_structural_event_call) {
                FinalBoundaryDecision boundary;
                boundary.pass =
                    boundary_evidence.canonical_pass ||
                    boundary_evidence.evidence_consistent;
                boundary.boundary_type = boundary_evidence.boundary_type;
                boundary.boundary_len = boundary_evidence.boundary_len;
                boundary.qc = boundary_evidence.qc;

                ComponentFinalCallEvaluation evaluation;
                evaluation.call = emit_final_te_call(
                    component,
                    event_evidence,
                    event_consensus,
                    event_segmentation,
                    te_alignment,
                    genotype,
                    boundary,
                    joint);
                evaluation.genotype = genotype;
                evaluation.score = candidate_total;
                evaluation.anchor_pos = shortlisted.validator.summary.bp_left >= 0 &&
                        shortlisted.validator.summary.bp_right >= 0
                    ? (shortlisted.validator.summary.bp_left +
                       ((shortlisted.validator.summary.bp_right -
                         shortlisted.validator.summary.bp_left) / 2))
                    : -1;
                evaluation.anchor_support =
                    shortlisted.validator.summary.alt_split_reads +
                    shortlisted.validator.summary.alt_indel_reads +
                    std::max(
                        shortlisted.validator.summary.alt_left_clip_reads,
                        shortlisted.validator.summary.alt_right_clip_reads);
                evaluation.anchor_ref_span_reads =
                    shortlisted.validator.summary.ref_span_reads;
                evaluation.anchor_priority =
                    shortlisted.validator.summary.hypothesis_priority;
                evaluation.anchor_hypothesis_score =
                    shortlisted.validator.summary.hypothesis_score;
                if (evaluation.call.pos >= 0 && event_evidence.bp_left >= 0) {
                    evaluation.retether_direction =
                        event_segmentation.left_flank_align_len > 0 &&
                                event_segmentation.right_flank_align_len == 0
                            ? 1
                            : (event_segmentation.right_flank_align_len > 0 &&
                                       event_segmentation.left_flank_align_len == 0
                                   ? -1
                                   : 0);
                }
                evaluation.evidence_te =
                    joint.emit_evidence_te_call || joint.emit_structural_event_call;
                evaluation.one_sided_segmentation =
                    event_segmentation.pass &&
                    !event_segmentation.insert_seq.empty() &&
                    ((event_segmentation.left_flank_align_len > 0) !=
                     (event_segmentation.right_flank_align_len > 0));
                final_call_candidates.push_back(std::move(evaluation));
            }
        }

        if (config_.log_stage_components && have_best) {
            std::ostringstream oss;
            const int32_t best_precise_reads =
                std::max(0, best_event_evidence.alt_split_reads) +
                std::max(0, best_event_evidence.alt_indel_reads);
            const int32_t best_bilateral_clip_reads = std::min(
                std::max(0, best_event_evidence.alt_left_clip_reads),
                std::max(0, best_event_evidence.alt_right_clip_reads));
            oss << "[Pipeline][joint]"
                << " chrom=" << component.chrom
                << " anchor=" << component.anchor_pos
                << " bp=" << best_event_evidence.bp_left << "-" << best_event_evidence.bp_right
                << " gt=" << display_or_na(best_genotype.genotype)
                << " af=" << best_genotype.af
                << " gq=" << best_genotype.gq
                << " alt_reads=" << best_event_evidence.alt_struct_reads
                << " alt_split=" << best_event_evidence.alt_split_reads
                << " alt_indel=" << best_event_evidence.alt_indel_reads
                << " alt_left_clip=" << best_event_evidence.alt_left_clip_reads
                << " alt_right_clip=" << best_event_evidence.alt_right_clip_reads
                << " precise_reads=" << best_precise_reads
                << " bilateral_clip_reads=" << best_bilateral_clip_reads
                << " ref_reads=" << best_event_evidence.ref_span_reads
                << " consensus_qc=" << display_or_na(best_event_consensus.qc_reason)
                << " consensus_len=" << best_event_consensus.consensus_len
                << " seg_qc=" << display_or_na(best_event_segmentation.qc_reason)
                << " seg_insert_len=" << best_event_segmentation.insert_seq.size()
                << " seg_left=" << best_event_segmentation.left_flank_align_len
                << " seg_right=" << best_event_segmentation.right_flank_align_len
                << " te_qc=" << display_or_na(best_te_alignment.qc_reason)
                << " te_pass=" << bool_as_int(best_te_alignment.pass)
                << " te_family=" << display_or_na(best_te_alignment.best_family)
                << " te_subfamily=" << display_or_na(best_te_alignment.best_subfamily)
                << " te_identity=" << best_te_alignment.best_identity
                << " te_cov=" << best_te_alignment.best_query_coverage
                << " te_coarse_prefilter=" << best_te_alignment.coarse_prefilter_score
                << " te_coarse_cov=" << best_te_alignment.coarse_chain_coverage
                << " te_margin=" << best_te_alignment.cross_family_margin
                << " te_model_label="
                << display_or_na(best_te_alignment.sequence_model_label)
                << " te_model_score="
                << best_te_alignment.sequence_model_score
                << " te_model_gc="
                << best_te_alignment.sequence_model_gc
                << " te_model_entropy="
                << best_te_alignment.sequence_model_entropy
                << " te_model_tandem="
                << best_te_alignment.sequence_model_tandem_fraction
                << " te_model_low_complexity="
                << best_te_alignment.sequence_model_low_complexity_fraction
                << " te_model_jsd_k5="
                << best_te_alignment.sequence_model_jsd_k5
                << " te_model_jsd_k6="
                << best_te_alignment.sequence_model_jsd_k6
                << " te_model_k9_containment="
                << best_te_alignment.sequence_model_k9_containment
                << " te_annotation_confidence="
                << display_or_na(best_te_alignment.annotation_confidence)
                << " te_annotation_residual="
                << best_te_alignment.annotation_residual_fraction
                << " te_annotation_masked="
                << best_te_alignment.annotation_masked_fraction
                << " te_explain_status="
                << te_annotation_status_name(best_te_alignment.te_sequence_explanation.status)
                << " te_explain_residual="
                << serialize_te_path_residual(best_te_alignment.te_sequence_explanation.residual)
                << " te_explain_path="
                << serialize_te_sequence_path(best_te_alignment.te_sequence_explanation)
                << " boundary_qc=" << display_or_na(best_boundary_evidence.qc)
                << " boundary_defined=" << bool_as_int(best_boundary_evidence.geometry_defined)
                << " boundary_canonical=" << bool_as_int(best_boundary_evidence.canonical_pass)
                << " boundary_consistent=" << bool_as_int(best_boundary_evidence.evidence_consistent)
                << " boundary_type=" << display_or_na(best_boundary_evidence.boundary_type)
                << " boundary_len=" << best_boundary_evidence.boundary_len
                << " clip_insert_qc=" << display_or_na(best_clip_insert_concordance.qc)
                << " clip_insert_full=" << best_clip_insert_concordance.full_insert_reads
                << " clip_insert_left=" << best_clip_insert_concordance.left_clip_reads
                << " clip_insert_right=" << best_clip_insert_concordance.right_clip_reads
                << " clip_insert_max_left_identity="
                << best_clip_insert_concordance.max_left_identity
                << " clip_insert_max_right_identity="
                << best_clip_insert_concordance.max_right_identity
                << " best_kind=" << final_hypothesis_kind_name(best_joint.best.kind)
                << " best_total=" << best_joint.best.total
                << " best_hard_veto=" << bool_as_int(best_joint.best.hard_veto)
                << " best_exist=" << best_joint.best.existence
                << " best_seg=" << best_joint.best.segmentation
                << " best_te=" << best_joint.best.te
                << " best_boundary=" << best_joint.best.boundary
                << " runner_kind=" << final_hypothesis_kind_name(best_joint.runner_up.kind)
                << " runner_total=" << best_joint.runner_up.total
                << " runner_hard_veto=" << bool_as_int(best_joint.runner_up.hard_veto)
                << " final_qc=" << display_or_na(best_joint.final_qc)
                << " emit_te=" << bool_as_int(best_joint.emit_te_call)
                << " emit_structural="
                << bool_as_int(best_joint.emit_structural_event_call)
                << " emit_unknown=" << bool_as_int(best_joint.emit_unknown_te);
            emit_pipeline_log_line(oss.str());
        }

        if (shortlist.empty()) {
            log_component(
                "NO_FEASIBLE_EVENT_HYPOTHESIS",
                component,
                fragments,
                nullptr,
                nullptr);
            continue;
        }

        if (!have_best || final_call_candidates.empty()) {
            log_component(
                have_best ? best_joint.final_qc.c_str() : "NO_BREAKPOINT_HYPOTHESIS",
                component,
                fragments,
                have_best ? &best_genotype : nullptr,
                nullptr);
            continue;
        }

        std::vector<ComponentFinalCallCandidate> selection_candidates;
        selection_candidates.reserve(
            final_call_candidates.size() + component_anchor_candidates.size());
        for (const auto& candidate : final_call_candidates) {
            ComponentFinalCallCandidate selection_candidate;
            selection_candidate.pos = candidate.call.pos;
            selection_candidate.anchor_pos = candidate.anchor_pos;
            selection_candidate.score = candidate.score;
            selection_candidate.emit_te = true;
            selection_candidate.evidence_te = candidate.evidence_te;
            selection_candidate.resolved_te =
                candidate.call.family_committed &&
                candidate.call.final_qc.rfind("PASS_STRUCTURAL_INSERTION", 0) != 0;
            selection_candidate.one_sided_segmentation = candidate.one_sided_segmentation;
            selection_candidate.anchor_support = candidate.anchor_support;
            selection_candidate.anchor_ref_span_reads = candidate.anchor_ref_span_reads;
            selection_candidate.anchor_priority = candidate.anchor_priority;
            selection_candidate.anchor_hypothesis_score = candidate.anchor_hypothesis_score;
            selection_candidate.component_anchor_pos = component.anchor_pos;
            selection_candidate.retether_direction = candidate.retether_direction;
            selection_candidates.push_back(selection_candidate);
        }
        selection_candidates.insert(
            selection_candidates.end(),
            component_anchor_candidates.begin(),
            component_anchor_candidates.end());
        retether_evidence_supported_final_call_positions(selection_candidates);
        for (size_t i = 0; i < selection_candidates.size(); ++i) {
            if (i >= final_call_candidates.size()) {
                break;
            }
            auto& call = final_call_candidates[i].call;
            const int32_t retethered_pos = selection_candidates[i].pos;
            if (retethered_pos >= 0 && retethered_pos != call.pos) {
                call.pos = retethered_pos;
                call.bp_left = retethered_pos;
                call.bp_right = retethered_pos;
                append_qc_token(call.final_qc, "PETA_RETETHERED_BREAKPOINT");
            }
        }
        const auto selected_indices = select_component_final_call_indices(selection_candidates);
        if (selected_indices.empty()) {
            log_component(
                best_joint.final_qc.c_str(),
                component,
                fragments,
                &best_genotype,
                nullptr);
            continue;
        }

        for (size_t selected_index : selected_indices) {
            const auto& selected = final_call_candidates[selected_index];
            result.final_calls.push_back(selected.call);
            bin_stats.final_calls += 1;

            log_component(
                selected.call.final_qc.c_str(),
                component,
                fragments,
                &selected.genotype,
                &result.final_calls.back());
        }
    }

    result.performance.component_build_seconds += bin_stats.component_build_seconds;
    result.performance.local_fetch_seconds += bin_stats.local_fetch_seconds;
    result.performance.projection_seconds += bin_stats.projection_seconds;
    result.performance.local_component_refresh_seconds +=
        bin_stats.local_component_refresh_seconds;
    result.performance.fragment_extract_seconds += bin_stats.fragment_extract_seconds;
    result.performance.signal_cache_build_seconds += bin_stats.signal_cache_build_seconds;
    result.performance.hypothesis_summary_seconds += bin_stats.hypothesis_summary_seconds;
    result.performance.validator_input_count_seconds +=
        bin_stats.validator_input_count_seconds;
    result.performance.event_consensus_seconds += bin_stats.timing.event_consensus_seconds;
    result.performance.event_segmentation_seconds += bin_stats.timing.event_segmentation_seconds;
    result.performance.te_alignment_seconds += bin_stats.timing.te_alignment_seconds;
    result.performance.bin_processing_wall_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - bin_wall_started).count();

    if (config_.log_stage_bins) {
        std::ostringstream oss;
            oss << "[Pipeline][bin][done]"
                << " chrom=" << chrom
                << " tid=" << tid
                << " start=" << bin_start
                << " end=" << bin_end
                << " records=" << bin_records.size()
                << " components=" << components.size()
                << " breakpoint_hypotheses_total=" << bin_stats.breakpoint_hypotheses_total
                << " breakpoint_hypotheses_after_collapse="
                << bin_stats.breakpoint_hypotheses_after_collapse
                << " breakpoint_hypotheses_after_gate=" << bin_stats.breakpoint_hypotheses_after_gate
                << " validator_candidates_total=" << bin_stats.validator_candidates_total
                << " validator_feasible_total=" << bin_stats.validator_feasible_total
                << " shortlist_empty=" << bin_stats.shortlist_empty
                << " shortlist_primary_only=" << bin_stats.shortlist_primary_only
                << " shortlist_with_challenger=" << bin_stats.shortlist_with_challenger
                << " expensive_stage_hypotheses_total=" << bin_stats.expensive_stage_hypotheses_total
                << " event_consensus_calls=" << bin_stats.event_consensus_calls
                << " event_consensus_rejected=" << bin_stats.event_consensus_rejected
                << " pre_segmentation_gate_rejected=" << bin_stats.pre_segmentation_gate_rejected
                << " pre_segmentation_gate_no_bilateral_anchor="
                << bin_stats.pre_segmentation_gate_no_bilateral_anchor
                << " pre_segmentation_gate_no_precise_or_full_context="
                << bin_stats.pre_segmentation_gate_no_precise_or_full_context
                << " event_segmentation_calls=" << bin_stats.event_segmentation_calls
                << " event_segmentation_rejected=" << bin_stats.event_segmentation_rejected
                << " te_alignment_calls=" << bin_stats.te_alignment_calls
                << " te_alignment_rejected=" << bin_stats.te_alignment_rejected
                << " component_build_seconds=" << bin_stats.component_build_seconds
                << " local_fetch_seconds=" << bin_stats.local_fetch_seconds
                << " projection_seconds=" << bin_stats.projection_seconds
                << " local_component_refresh_seconds=" << bin_stats.local_component_refresh_seconds
                << " fragment_extract_seconds=" << bin_stats.fragment_extract_seconds
                << " signal_cache_build_seconds=" << bin_stats.signal_cache_build_seconds
                << " hypothesis_summary_seconds=" << bin_stats.hypothesis_summary_seconds
                << " validator_input_count_seconds=" << bin_stats.validator_input_count_seconds
                << " event_evidence_seconds=" << bin_stats.event_evidence_seconds
                << " event_consensus_seconds=" << bin_stats.timing.event_consensus_seconds
                << " event_segmentation_seconds=" << bin_stats.timing.event_segmentation_seconds
                << " te_alignment_seconds=" << bin_stats.timing.te_alignment_seconds
                << " genotype_calls=" << bin_stats.genotype_calls
                << " final_calls=" << bin_stats.final_calls
                << " local_reuse_ratio=" << local_interval_reuse_ratio(reuse_stats);
        emit_pipeline_log_line(oss.str());
    }
}
