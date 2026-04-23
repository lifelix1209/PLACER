#ifndef PLACER_TE_FAMILY_ALIGNMENT_H
#define PLACER_TE_FAMILY_ALIGNMENT_H

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace placer {

using KmerPositionMap = std::unordered_map<uint64_t, std::vector<int32_t>>;

struct TeFamilyRepresentative {
    int32_t te_id = -1;
    std::string family;
    std::string family_key;
    std::string exact_name;
    std::string sequence;
    std::string reverse_complement_sequence;
    KmerPositionMap forward_kmer_positions;
    KmerPositionMap reverse_kmer_positions;
};

struct TeFamilyGroup {
    std::string family;
    std::string family_key;
    std::vector<int32_t> member_te_ids;
    std::vector<TeFamilyRepresentative> representatives;
};

struct TeFamilyGroupCache {
    std::vector<TeFamilyGroup> groups;
};

struct TeFamilyAlignmentIndex {
    int32_t kmer_size = 9;
    std::vector<TeFamilyGroup> groups;
    std::unordered_map<std::string, size_t> family_to_group;
    std::unordered_map<uint64_t, int32_t> kmer_family_counts;
    std::unordered_map<uint64_t, double> kmer_family_weights;
    std::vector<uint64_t> representative_seed_keys;
};

struct TeFamilyChainScore {
    int32_t group_index = -1;
    int32_t representative_index = -1;
    bool reverse = false;
    double best_rep_chain_norm = 0.0;
    double second_rep_chain_norm = 0.0;
    double chained_query_coverage = 0.0;
    double family_prefilter_score = 0.0;
};

struct TemplateBandEstimate {
    double seed_score = 0.0;
    double query_coverage = 0.0;
    int32_t band_center = 0;
    int32_t band_width = 32;
};

struct ExactAlignmentSummary {
    double identity = 0.0;
    double query_coverage = 0.0;
    double score = 0.0;
};

struct TeFamilyCacheBundle {
    std::shared_ptr<const TeFamilyAlignmentIndex> alignment_index;
    std::shared_ptr<const TeFamilyGroupCache> rep_groups;
};

TeFamilyCacheBundle build_te_family_cache(
    const std::vector<std::string>& te_names,
    const std::vector<std::string>& te_sequences,
    const std::vector<std::string>& te_reverse_complement_sequences,
    int32_t family_representatives,
    int32_t kmer_size = 9);

std::vector<TeFamilyChainScore> score_family_prefilter(
    const TeFamilyAlignmentIndex& index,
    const std::string& query,
    int32_t kmer_size,
    int32_t family_topn);

TemplateBandEstimate estimate_template_band(
    const std::string& query,
    const std::string& target,
    int32_t kmer_size);

TemplateBandEstimate estimate_template_band_from_positions(
    const KmerPositionMap& query_positions,
    int32_t query_len,
    const KmerPositionMap& target_positions,
    int32_t target_len,
    int32_t kmer_size);

ExactAlignmentSummary banded_semiglobal_affine_align(
    const std::string& query,
    const std::string& target,
    int32_t band_center,
    int32_t band_width);

}  // namespace placer

#endif  // PLACER_TE_FAMILY_ALIGNMENT_H
