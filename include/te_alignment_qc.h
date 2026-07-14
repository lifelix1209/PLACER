#ifndef PLACER_TE_ALIGNMENT_QC_H
#define PLACER_TE_ALIGNMENT_QC_H

#include <string>

namespace placer {

// Strongly-typed classification of TEAlignmentEvidence::qc_reason for use in
// decision control flow.
//
// The qc_reason string field stays the source of truth -- it is serialised into
// the output QC token, so its exact spelling is part of the output contract.
// This enum lets the gating logic in decision_policy.cpp branch on a typed value
// instead of comparing raw string literals, so a mistyped reason is a compile
// error rather than a silently-false comparison. Reasons not relevant to gating
// map to kOther.
enum class TeAlignmentQc {
    kOther,
    kNoTeAlignment,
    kPassInsert,
    kPassInsertFamilyOnly,
    kPassInsertUnknown,
    kLowIdentity,
};

inline TeAlignmentQc classify_te_alignment_qc(const std::string& qc_reason) {
    if (qc_reason == "NO_TE_ALIGNMENT") {
        return TeAlignmentQc::kNoTeAlignment;
    }
    if (qc_reason == "PASS_INSERT_TE_ALIGNMENT") {
        return TeAlignmentQc::kPassInsert;
    }
    if (qc_reason == "PASS_INSERT_TE_ALIGNMENT_FAMILY_ONLY") {
        return TeAlignmentQc::kPassInsertFamilyOnly;
    }
    if (qc_reason == "PASS_INSERT_TE_ALIGNMENT_UNKNOWN") {
        return TeAlignmentQc::kPassInsertUnknown;
    }
    if (qc_reason == "TE_ALIGNMENT_LOW_IDENTITY") {
        return TeAlignmentQc::kLowIdentity;
    }
    return TeAlignmentQc::kOther;
}

// True for any placed PASS_INSERT_TE_ALIGNMENT* reason (resolved or unknown).
inline bool is_pass_insert_te_alignment(TeAlignmentQc qc) {
    return qc == TeAlignmentQc::kPassInsert ||
           qc == TeAlignmentQc::kPassInsertFamilyOnly ||
           qc == TeAlignmentQc::kPassInsertUnknown;
}

// True for the resolved (non-unknown) PASS_INSERT_TE_ALIGNMENT* reasons.
inline bool is_resolved_pass_insert_te_alignment(TeAlignmentQc qc) {
    return qc == TeAlignmentQc::kPassInsert ||
           qc == TeAlignmentQc::kPassInsertFamilyOnly;
}

}  // namespace placer

#endif  // PLACER_TE_ALIGNMENT_QC_H
