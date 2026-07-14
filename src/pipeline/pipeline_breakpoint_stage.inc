std::vector<Pipeline::BreakpointHypothesis> Pipeline::enumerate_breakpoint_hypotheses(
    const ComponentCall& component,
    const std::vector<const bam1_t*>& local_records,
    const std::vector<InsertionFragment>& fragments,
    int32_t seed_left,
    int32_t seed_right,
    size_t top_k) const {
    const auto local_hypotheses = collect_breakpoint_hypotheses(
        component,
        local_records,
        fragments,
        seed_left,
        seed_right,
        0);
    std::vector<BreakpointHypothesis> hypotheses;
    hypotheses.reserve(local_hypotheses.size());
    for (const auto& local : local_hypotheses) {
        BreakpointHypothesis hypothesis;
        hypothesis.valid = local.valid;
        hypothesis.left = local.left;
        hypothesis.right = local.right;
        hypothesis.center = local.center;
        hypothesis.support = local.support;
        hypothesis.priority = local.priority;
        hypotheses.push_back(hypothesis);
    }
    return select_diverse_breakpoint_hypotheses(hypotheses, top_k, component.anchor_pos);
}

std::vector<Pipeline::BreakpointHypothesis> Pipeline::select_diverse_breakpoint_hypotheses(
    const std::vector<BreakpointHypothesis>& hypotheses,
    size_t top_k,
    int32_t anchor_pos) const {
    constexpr int32_t kBreakpointDiversitySlackBp = 30;
    constexpr int32_t kAnchorProximalSlackBp = 100;
    constexpr int32_t kAnchorProximalMinScore = 8;
    constexpr int32_t kRescuePreciseAnchorMinDistanceBp = 500;
    if (top_k == 0 || hypotheses.size() <= top_k) {
        return hypotheses;
    }

    const auto same_locus = [](const BreakpointHypothesis& lhs, const BreakpointHypothesis& rhs) {
        return std::abs(lhs.left - rhs.left) <= kBreakpointDiversitySlackBp &&
            std::abs(lhs.right - rhs.right) <= kBreakpointDiversitySlackBp;
    };
    const auto hypothesis_score = [](const BreakpointHypothesis& hypothesis) {
        return hypothesis.support * breakpoint_hypothesis_support_weight(hypothesis.priority);
    };
    const auto selected_contains = [&](const std::vector<BreakpointHypothesis>& selected,
                                       const BreakpointHypothesis& hypothesis) {
        return std::any_of(
            selected.begin(),
            selected.end(),
            [&](const auto& incumbent) {
                return incumbent.left == hypothesis.left && incumbent.right == hypothesis.right;
            });
    };

    std::vector<BreakpointHypothesis> selected;
    selected.reserve(std::min(top_k, hypotheses.size()));
    for (const auto& hypothesis : hypotheses) {
        bool redundant = false;
        for (const auto& incumbent : selected) {
            if (same_locus(hypothesis, incumbent)) {
                redundant = true;
                break;
            }
        }
        if (redundant) {
            continue;
        }
        selected.push_back(hypothesis);
        if (selected.size() >= top_k) {
            break;
        }
    }

    const bool has_anchor_proximal = std::any_of(
        selected.begin(),
        selected.end(),
        [&](const auto& hypothesis) {
            return std::abs(hypothesis.center - anchor_pos) <= kAnchorProximalSlackBp;
        });
    const auto append_rescue_precise_anchor = [&]() {
        const BreakpointHypothesis* rescue_precise_anchor = nullptr;
        for (const auto& hypothesis : hypotheses) {
            if (selected_contains(selected, hypothesis)) {
                continue;
            }
            if ((hypothesis.priority != 0 && hypothesis.priority != 1) ||
                hypothesis.left != hypothesis.right ||
                hypothesis.support > 1) {
                continue;
            }
            const bool spatially_distinct = std::all_of(
                selected.begin(),
                selected.end(),
                [&](const auto& incumbent) {
                    return std::abs(hypothesis.center - incumbent.center) >=
                        kRescuePreciseAnchorMinDistanceBp;
                });
            if (!spatially_distinct) {
                continue;
            }
            if (rescue_precise_anchor == nullptr) {
                rescue_precise_anchor = &hypothesis;
                continue;
            }
            const int32_t dist = std::abs(hypothesis.center - anchor_pos);
            const int32_t best_dist =
                std::abs(rescue_precise_anchor->center - anchor_pos);
            if (dist < best_dist ||
                (dist == best_dist &&
                 hypothesis.priority < rescue_precise_anchor->priority) ||
                (dist == best_dist &&
                 hypothesis.priority == rescue_precise_anchor->priority &&
                 hypothesis.center < rescue_precise_anchor->center)) {
                rescue_precise_anchor = &hypothesis;
            }
        }
        if (rescue_precise_anchor) {
            selected.push_back(*rescue_precise_anchor);
        }
    };
    if (has_anchor_proximal) {
        append_rescue_precise_anchor();
        return selected;
    }

    const BreakpointHypothesis* best_anchor_candidate = nullptr;
    for (const auto& hypothesis : hypotheses) {
        if (selected_contains(selected, hypothesis)) {
            continue;
        }
        if (hypothesis_score(hypothesis) < kAnchorProximalMinScore) {
            continue;
        }
        if (best_anchor_candidate == nullptr) {
            best_anchor_candidate = &hypothesis;
            continue;
        }
        const int32_t dist = std::abs(hypothesis.center - anchor_pos);
        const int32_t best_dist = std::abs(best_anchor_candidate->center - anchor_pos);
        if (dist != best_dist) {
            if (dist < best_dist) {
                best_anchor_candidate = &hypothesis;
            }
            continue;
        }
        const int32_t score = hypothesis_score(hypothesis);
        const int32_t best_score = hypothesis_score(*best_anchor_candidate);
        if (score != best_score) {
            if (score > best_score) {
                best_anchor_candidate = &hypothesis;
            }
            continue;
        }
        if (hypothesis.support != best_anchor_candidate->support) {
            if (hypothesis.support > best_anchor_candidate->support) {
                best_anchor_candidate = &hypothesis;
            }
            continue;
        }
        if (hypothesis.priority != best_anchor_candidate->priority) {
            if (hypothesis.priority < best_anchor_candidate->priority) {
                best_anchor_candidate = &hypothesis;
            }
            continue;
        }
        if (hypothesis.left < best_anchor_candidate->left ||
            (hypothesis.left == best_anchor_candidate->left &&
             hypothesis.right < best_anchor_candidate->right)) {
            best_anchor_candidate = &hypothesis;
        }
    }
    if (best_anchor_candidate == nullptr) {
        return selected;
    }

    if (selected.size() < top_k) {
        selected.push_back(*best_anchor_candidate);
    } else if (!selected.empty()) {
        selected.back() = *best_anchor_candidate;
    }
    append_rescue_precise_anchor();
    return selected;
}

std::pair<int32_t, int32_t> Pipeline::resolve_event_breakpoint_bounds(
    const ComponentCall& component,
    const std::vector<const bam1_t*>& local_records,
    const std::vector<InsertionFragment>& fragments,
    int32_t seed_left,
    int32_t seed_right) const {
    const auto hypotheses = enumerate_breakpoint_hypotheses(
        component,
        local_records,
        fragments,
        seed_left,
        seed_right,
        1);
    if (!hypotheses.empty()) {
        return {hypotheses.front().left, hypotheses.front().right};
    }
    return {std::max(0, component.anchor_pos), std::max(0, component.anchor_pos)};
}
