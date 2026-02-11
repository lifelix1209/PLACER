#ifndef PLACER_GENOTYPING_H
#define PLACER_GENOTYPING_H

#include "placeability.h"
#include "assembly.h"
#include "local_realign.h"
#include <vector>
#include <string>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <random>
#include <unordered_map>
#include <array>

namespace placer {

// ============================================================================
// Genotyping Configuration
// ============================================================================

struct GenotypeConfig {
    // EM 迭代参数
    int max_em_iterations = 20;           // 最大迭代次数
    double em_convergence_threshold = 1e-6; // 收敛阈值 (对数似然增益)
    double em_init_alt_fraction = 0.5;   // ALT 初始比例

    // 空间先验参数
    double spatial_lambda = 100.0;       // 距离衰减参数 λ (bp)
    double null_base_prior = 0.1;        // NULL 基础先验
    double prior_min_bound = 1e-6;       // prior_{ik} 下限（防止某类被完全打死）

    // 结构先验参数
    double family_match_bonus = 2.0;      // TE 家族匹配时的先验加成
    double family_mismatch_penalty = 0.1;  // TE 家族不匹配时的惩罚

    // 似然权重
    double geom_ok_weight = 2.0;         // 几何一致性权重
    double contig_support_weight = 1.5;    // Contig 支持度权重
    double align_score_weight = 1.0;      // 对齐分数权重

    // 过滤阈值
    double min_likelihood_ratio = 1e-10; // 最小似然比（避免数值下溢）
    int min_reads_for_genotyping = 2;   // 最小支持 reads

    // Genotype-level 参数
    double genotype_error_rate = 0.01;   // 基因型错误率 (eps)
    double population_af = 0.001;        // 群体等位基因频率 (用于 HWE prior)

    // LRT 检验参数
    double lrt_threshold = 5.99;        // LRT 阈值 (chi2 df=2, p=0.05)

    // 输出参数
    bool output_confidence_interval = true; // 输出置信区间
    double confidence_level = 0.95;       // 置信水平
};

// ============================================================================
// Genotyping Results
// ============================================================================

struct GenotypeResult {
    // ========== Global Mixture Weights π ==========
    // 注意：mix_* 是 EM 估计的全局混合权重 π，不是等位基因频率 AF！
    // π_alt + π_ref + π_null = 1
    double mix_alt = 0.0;      // π_alt: ALT 结构成分的全局权重
    double mix_ref = 0.0;       // π_ref: REF (无插入) 成分的全局权重
    double mix_null = 0.0;      // π_null: NULL (背景/不可解释) 成分的全局权重

    // ========== Expected Counts ==========
    // 每个成分对所有 reads 的期望责任度之和
    double e_alt_sum = 0.0;     // Σ_i γ_{i,ALT}
    double e_ref_sum = 0.0;     // Σ_i γ_{i,REF}
    double e_null_sum = 0.0;     // Σ_i γ_{i,NULL}

    // ========== Allele Frequency (Derived) ==========
    // AF 从 π 导出：AF = π_alt / (π_alt + π_ref)，NULL 不计入分母
    double af = 0.0;               // ALT 频率 = E[ALT] / (E[ALT] + E[REF])
    double af_ci_low = 0.0;
    double af_ci_high = 0.0;

    // ========== Genotype判定 ==========
    std::string genotype = "./.";  // 0/0, 0/1, 1/1, ./.
    int gq = 0;                    // 基因型质量 (Phred)
    std::array<double, 3> gt_log_lik = {0.0, 0.0, 0.0};  // log P(data | GT)

    // ========== EM 统计 ==========
    double log_likelihood = 0.0;    // 最终对数似然 log P(Data | π)
    int em_iterations = 0;          // 实际迭代次数
    bool converged = false;         // 是否收敛

    // ========== LRT 检验 ==========
    double lrt_statistic = 0.0;    // 似然比统计量
    double likelihood_ratio = 0.0;  // 似然比 (exp(LRT))
    bool lrt_significant = false;  // LRT 是否显著

    // ========== 标记 ==========
    bool high_background = false;   // HIGH_BACKGROUND 标记
    bool low_complexity = false;    // LOW_COMPLEXITY 标记
    bool hotspot = false;          // HOTSPOT 标记

    std::string to_string() const;
};

// ============================================================================
// Read Evidence for Genotyping
// ============================================================================

struct ReadEvidence {
    size_t read_idx = 0;

    // ========== 空间特征 ==========
    int32_t primary_pos = -1;       // Primary alignment 位置
    int32_t locus_pos = -1;         // 最佳候选 locus 位置
    double d_spatial = 0.0;         // 与候选位点的距离

    // ========== 几何一致性 ==========
    bool geom_ok = false;          // 断点几何解释一致
    double geom_score = 0.0;        // 几何一致性分数 [0, 1]

    // ========== REF 跨越特征 ==========
    bool ref_span_ok = false;      // Read 是否跨越断点区间
    double ref_span_score = 0.0;    // REF 跨越分数 [0, 1]

    // ========== 对齐特征 ==========
    double align_score = 0.0;       // 对齐分数
    double normalized_score = 0.0; // 归一化分数
    float alignment_identity = 0.0f;  // 对齐相似度

    // ========== Contig 支持 ==========
    bool contig_support = false;    // 是否与代表 contig 一致
    double contig_score = 0.0;     // Contig 支持度 [0, 1]

    // ========== TE 家族 ==========
    int te_family_id = -1;          // TE 家族 ID
    bool family_match = false;      // Read TE 家族是否匹配代表
    bool family_mismatch = false;   // Read TE 家族是否不匹配代表

    // ========== 结构评分（用于空间先验） ==========
    double structural_score = 0.5;  // 综合结构评分 [0, 1]，由 structural_prior 计算

    // ========== 旧别名（兼容） ==========
    double identity = 0.0f;         // 对齐相似度 (alignment_identity 的别名)
};

// ============================================================================
// Prior Calculators
// ============================================================================

class SpatialPriorCalculator {
public:
    explicit SpatialPriorCalculator(const GenotypeConfig& config);

    // 计算空间先验
    void calculate_prior(double distance, double structural_score,
                        double& pi_alt, double& pi_ref, double& pi_null) const;

    // 批量计算（返回 pi_alt, pi_ref, pi_null）
    std::vector<std::array<double, 3>> calculate_priors(
        const std::vector<ReadEvidence>& evidence) const;

private:
    GenotypeConfig config_;
};

class StructuralPriorCalculator {
public:
    explicit StructuralPriorCalculator(const GenotypeConfig& config);

    // 计算家族加成
    double calculate_family_bonus(int read_family, int rep_family) const;

    // 计算结构评分
    double calculate_structural_score(
        const ReadEvidence& e, int representative_family) const;

    // 批量计算结构评分
    std::vector<double> calculate_structural_scores(
        const std::vector<ReadEvidence>& evidence,
        int representative_family) const;

private:
    GenotypeConfig config_;
};

// ============================================================================
// EM Engine
// ============================================================================

class EMEngine {
public:
    explicit EMEngine(const GenotypeConfig& config);

    // 运行 EM 迭代
    // 返回: (log_likelihood, iterations)
    std::pair<double, int> run_em(
        const std::vector<ReadEvidence>& evidence,
        const std::vector<std::array<double, 3>>& priors,
        GenotypeResult& result) const;

    // E-step: 计算责任度
    void e_step(
        const std::vector<ReadEvidence>& evidence,
        const std::vector<std::array<double, 3>>& priors,
        double pi_alt, double pi_ref, double pi_null,
        std::vector<std::array<double, 3>>& responsibilities) const;

    // M-step: 更新参数
    void m_step(
        const std::vector<std::array<double, 3>>& responsibilities,
        double& pi_alt, double& pi_ref, double& pi_null) const;

    // 计算对数似然
    double calculate_log_likelihood(
        const std::vector<ReadEvidence>& evidence,
        const std::vector<std::array<double, 3>>& priors,
        double pi_alt, double pi_ref, double pi_null) const;

    // 计算 genotype likelihoods
    void calculate_genotype_likelihoods(
        double e_alt_sum, double e_ref_sum, int n_reads,
        std::array<double, 3>& gt_log_lik) const;

    // 从 genotype likelihoods 计算 GQ
    int calculate_gq_from_gt(
        const std::array<double, 3>& gt_log_lik,
        std::string& best_gt) const;

private:
    GenotypeConfig config_;

    // 计算单个似然
    double likelihood_alt(const ReadEvidence& e) const;
    double likelihood_ref(const ReadEvidence& e) const;
    double likelihood_null(const ReadEvidence& e) const;

    // Beta-Binomial 后验区间
    void calculate_confidence_interval(
        double alpha, double beta, double level,
        double& ci_low, double& ci_high) const;
};

// ============================================================================
// Genotyper (Main Class)
// ============================================================================

class Genotyper {
public:
    explicit Genotyper(const GenotypeConfig& config = GenotypeConfig());

    // 主分型接口
    GenotypeResult genotype(
        const StructuralRepresentative& representative,
        const std::vector<LocusEvidence>& evidence,
        const std::vector<ReadSketch>& reads,
        const GenomeAccessor& genome) const;

    // 批量分型
    std::vector<GenotypeResult> genotype_batch(
        const std::vector<StructuralRepresentative>& representatives,
        const std::vector<std::vector<LocusEvidence>>& all_evidence,
        const std::vector<ReadSketch>& all_reads,
        const GenomeAccessor& genome) const;

    // 从 LocusEvidence 构建 ReadEvidence
    std::vector<ReadEvidence> build_evidence(
        const StructuralRepresentative& representative,
        const std::vector<LocusEvidence>& locus_evidence,
        const std::vector<ReadSketch>& reads,
        const GenomeAccessor& genome) const;

private:
    GenotypeConfig config_;
    SpatialPriorCalculator spatial_prior_;
    StructuralPriorCalculator structural_prior_;
    EMEngine em_engine_;

    // 构建空间特征
    void compute_spatial_features(
        const StructuralRepresentative& rep,
        const std::vector<LocusEvidence>& evidence,
        std::vector<ReadEvidence>& out_evidence) const;

    // 构建几何特征
    void compute_geom_features(
        const StructuralRepresentative& rep,
        const std::vector<LocusEvidence>& evidence,
        std::vector<ReadEvidence>& out_evidence) const;

    // 构建对齐特征
    void compute_align_features(
        const std::vector<LocusEvidence>& evidence,
        std::vector<ReadEvidence>& out_evidence) const;

    // 构建 REF 跨越特征
    void compute_ref_span_features(
        const StructuralRepresentative& rep,
        const std::vector<LocusEvidence>& evidence,
        const std::vector<ReadSketch>& reads,
        const GenomeAccessor& genome,
        std::vector<ReadEvidence>& out_evidence) const;

    // 构建 contig 支持特征
    void compute_contig_features(
        const StructuralRepresentative& rep,
        const std::vector<LocusEvidence>& evidence,
        const std::vector<ReadSketch>& reads,
        std::vector<ReadEvidence>& out_evidence) const;

    // 构建 TE 家族特征
    void compute_family_features(
        const StructuralRepresentative& rep,
        const std::vector<LocusEvidence>& evidence,
        const std::vector<ReadSketch>& reads,
        std::vector<ReadEvidence>& out_evidence) const;

    // 判定基因型
    std::string determine_genotype(double e_alt, double e_ref, double e_null) const;

    // 计算 GQ (Phred-scaled)
    int calculate_gq(double best_prob, double second_best_prob) const;

    // 计算 AF CI
    void calculate_af_ci(double e_alt, double e_ref, double e_null, GenotypeResult& result) const;
};

// ============================================================================
// Utility Functions
// ============================================================================

// Beta-Binomial 计算
double beta_pdf(double x, double alpha, double beta);
double beta_cdf(double x, double alpha, double beta);
std::array<double, 2> beta_binomial_ci(int successes, int trials, double level);

// Phred 转换
double phred_to_prob(int phred);
int prob_to_phred(double prob);

// 似然比检验
double likelihood_ratio_test(double log_lik_null, double log_lik_alt);

}  // namespace placer

#endif  // PLACER_GENOTYPING_H
