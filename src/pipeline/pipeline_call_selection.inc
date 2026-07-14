std::vector<size_t> select_component_final_call_indices(
    const std::vector<ComponentFinalCallCandidate>& candidates) {
    std::vector<size_t> order;
    order.reserve(candidates.size());
    for (size_t i = 0; i < candidates.size(); ++i) {
        if (candidates[i].emit_te) {
            order.push_back(i);
        }
    }
    if (order.empty()) {
        return {};
    }

    std::sort(order.begin(), order.end(), [&](size_t lhs, size_t rhs) {
        if (candidates[lhs].pos != candidates[rhs].pos) {
            return candidates[lhs].pos < candidates[rhs].pos;
        }
        return lhs < rhs;
    });

    std::vector<size_t> selected;
    selected.reserve(order.size());
    size_t cluster_start = 0;
    while (cluster_start < order.size()) {
        const int32_t anchor_pos = candidates[order[cluster_start]].pos;
        size_t best = order[cluster_start];
        size_t cluster_end = cluster_start;
        while (cluster_end < order.size() &&
               std::abs(candidates[order[cluster_end]].pos - anchor_pos) <=
                   kFinalCallDedupDistanceBp) {
            const size_t idx = order[cluster_end];
            if (candidates[idx].score > candidates[best].score ||
                (candidates[idx].score == candidates[best].score && idx < best)) {
                best = idx;
            }
            ++cluster_end;
        }
        selected.push_back(best);
        cluster_start = cluster_end;
    }

    std::sort(selected.begin(), selected.end(), [&](size_t lhs, size_t rhs) {
        if (candidates[lhs].pos != candidates[rhs].pos) {
            return candidates[lhs].pos < candidates[rhs].pos;
        }
        return lhs < rhs;
    });
    return selected;
}

void retether_evidence_supported_final_call_positions(
    std::vector<ComponentFinalCallCandidate>& candidates) {
    constexpr int32_t kMinRetetherAnchorSupport = 1;
    constexpr int32_t kMinRetetherDistanceBp = 100;
    constexpr int32_t kMaxRetetherDistanceBp = 3000;
    constexpr int32_t kCoreAnchorRadiusBp = 100;
    constexpr int32_t kLocalAnchorRadiusBp = 450;
    constexpr int32_t kMaxAnchorComponentDistanceBp = 1000;

    const auto is_precise_anchor = [](int32_t priority) {
        return priority == 0 || priority == 1 || priority == 2 || priority == 4;
    };
    const auto retether_anchor_score = [&](const ComponentFinalCallCandidate& candidate,
                                           const ComponentFinalCallCandidate& anchor) {
        const int32_t distance = std::abs(candidate.pos - anchor.pos);
        const int32_t component_distance =
            candidate.component_anchor_pos >= 0
                ? std::abs(anchor.pos - candidate.component_anchor_pos)
                : std::numeric_limits<int32_t>::max();
        double score = anchor.anchor_hypothesis_score;
        if (is_precise_anchor(anchor.anchor_priority)) {
            score += 1000.0;
        }
        if (component_distance <= kCoreAnchorRadiusBp) {
            score += 2000.0;
        }
        if (component_distance <= kLocalAnchorRadiusBp) {
            score += 50.0;
        }
        score -= static_cast<double>(distance) * 0.01;
        score -= static_cast<double>(component_distance) * 0.001;
        return score;
    };

    for (auto& candidate : candidates) {
        if (!candidate.emit_te ||
            !candidate.evidence_te ||
            candidate.pos < 0) {
            continue;
        }
        // A one-sided event is already positioned by its consensus-to-reference
        // flank. Moving it to an unrelated component hypothesis has no sequence
        // certificate and can turn a correctly anchored call into a distant
        // duplicate. Only two-sided calls remain eligible for candidate-level
        // retethering here.
        if (candidate.one_sided_segmentation) {
            continue;
        }
        if (!candidate.one_sided_segmentation) {
            if (candidate.anchor_pos >= 0 &&
                is_precise_anchor(candidate.anchor_priority) &&
                candidate.anchor_support >= kMinRetetherAnchorSupport &&
                (candidate.anchor_hypothesis_score >= 8.0 ||
                 candidate.resolved_te)) {
                const int32_t self_anchor_distance =
                    std::abs(candidate.pos - candidate.anchor_pos);
                if (self_anchor_distance <= kMaxRetetherDistanceBp &&
                    (self_anchor_distance >= kMinRetetherDistanceBp ||
                     candidate.resolved_te)) {
                    candidate.pos = candidate.anchor_pos;
                }
            }
        }
        if (candidate.resolved_te &&
            is_precise_anchor(candidate.anchor_priority) &&
            candidate.anchor_pos >= 0 &&
            candidate.anchor_support >= kMinRetetherAnchorSupport) {
            const int32_t self_anchor_distance =
                std::abs(candidate.pos - candidate.anchor_pos);
            if (self_anchor_distance <= kMaxRetetherDistanceBp) {
                candidate.pos = candidate.anchor_pos;
                continue;
            }
        }
        if (!candidate.resolved_te &&
            candidate.anchor_pos >= 0 &&
            candidate.anchor_support >= 2 &&
            candidate.anchor_hypothesis_score >= 10.0 &&
            candidate.anchor_ref_span_reads > 0) {
            const int32_t self_anchor_distance =
                std::abs(candidate.pos - candidate.anchor_pos);
            const int32_t self_component_distance =
                candidate.component_anchor_pos >= 0
                    ? std::abs(candidate.anchor_pos - candidate.component_anchor_pos)
                    : std::numeric_limits<int32_t>::max();
            if (self_component_distance <= kLocalAnchorRadiusBp &&
                self_anchor_distance <= kMaxRetetherDistanceBp) {
                if (self_anchor_distance >= kMinRetetherDistanceBp) {
                    candidate.pos = candidate.anchor_pos;
                }
                continue;
            }
        }
        int32_t tether_pos = candidate.anchor_pos;
        int32_t tether_support = candidate.anchor_support;
        int32_t tether_ref_span = candidate.anchor_ref_span_reads;
        int32_t tether_priority = candidate.anchor_priority;
        double tether_hypothesis_score = candidate.anchor_hypothesis_score;
        double best_tether_score = -std::numeric_limits<double>::infinity();

        for (const auto& anchor : candidates) {
            if (anchor.emit_te || anchor.pos < 0) {
                continue;
            }
            if (anchor.anchor_support < kMinRetetherAnchorSupport ||
                anchor.anchor_hypothesis_score <= 0.0) {
                continue;
            }
            if (!candidate.one_sided_segmentation &&
                !is_precise_anchor(anchor.anchor_priority)) {
                continue;
            }
            const int32_t distance = std::abs(candidate.pos - anchor.pos);
            if (distance < kMinRetetherDistanceBp ||
                distance > kMaxRetetherDistanceBp) {
                continue;
            }
            if (candidate.component_anchor_pos >= 0 &&
                std::abs(anchor.pos - candidate.component_anchor_pos) >
                    kMaxAnchorComponentDistanceBp) {
                continue;
            }
            if (!is_precise_anchor(anchor.anchor_priority) &&
                candidate.anchor_hypothesis_score > 0.0 &&
                !(candidate.component_anchor_pos >= 0 &&
                  std::abs(anchor.pos - candidate.component_anchor_pos) <=
                      kCoreAnchorRadiusBp) &&
                anchor.anchor_hypothesis_score + 1e-9 <
                    candidate.anchor_hypothesis_score) {
                continue;
            }
            if (candidate.retether_direction < 0 && anchor.pos > candidate.pos) {
                continue;
            }
            if (candidate.retether_direction > 0 && anchor.pos < candidate.pos) {
                continue;
            }
            const double score = retether_anchor_score(candidate, anchor);
            if (score > best_tether_score ||
                (score == best_tether_score &&
                 anchor.anchor_hypothesis_score > tether_hypothesis_score) ||
                (score == best_tether_score &&
                 anchor.anchor_hypothesis_score == tether_hypothesis_score &&
                 anchor.anchor_support > tether_support) ||
                (score == best_tether_score &&
                 anchor.anchor_hypothesis_score == tether_hypothesis_score &&
                 anchor.anchor_support == tether_support &&
                 anchor.anchor_ref_span_reads < tether_ref_span) ||
                (score == best_tether_score &&
                 anchor.anchor_hypothesis_score == tether_hypothesis_score &&
                 anchor.anchor_support == tether_support &&
                 anchor.anchor_ref_span_reads == tether_ref_span &&
                 anchor.pos < tether_pos)) {
                tether_pos = anchor.pos;
                tether_support = anchor.anchor_support;
                tether_ref_span = anchor.anchor_ref_span_reads;
                tether_priority = anchor.anchor_priority;
                tether_hypothesis_score = anchor.anchor_hypothesis_score;
                best_tether_score = score;
            }
        }

        if (tether_support < kMinRetetherAnchorSupport || tether_pos < 0) {
            continue;
        }
        const int32_t distance = std::abs(candidate.pos - tether_pos);
        if (distance < kMinRetetherDistanceBp ||
            distance > kMaxRetetherDistanceBp) {
            continue;
        }
        candidate.pos = tether_pos;
    }
}
