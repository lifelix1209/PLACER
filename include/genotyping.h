#ifndef PLACER_GENOTYPING_H
#define PLACER_GENOTYPING_H

#include "placeability.h"
#include "assembly.h"
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

    // 输出参数
    bool output_confidence_interval = true; // 输出置信区间
    double confidence_level = 0.95;       // 置信水平
};

// ============================================================================
// Genotyping Results
// ============================================================================

struct GenotypeResult {
    // 混合模型比例（期望值）
    double e_alt = 0.0;      // ALT 结构比例
    double e_ref = 0.0;       // REF (无插入) 比例
    double e_null = 0.0;       // NULL (背景/不可解释) 比例

    // 后验置信区间
    double alt_ci_low = 0.0;
    double alt_ci_high = 0.0;

    // 基因型判定
    std::string genotype = "./.";  // 0/0, 0/1, 1/1, ./.
    int gq = 0;                    // 基因型质量 (Phred)

    // 等位基因频率
    double af = 0.0;               // ALT 频率 = E[ALT] / (E[ALT] + E[REF])
    double af_ci_low = 0.0;
    double af_ci_high = 0.0;

    // 详细信息
    double log_likelihood = 0.0;    // 最终对数似然
    int em_iterations = 0;          // 实际迭代次数
    bool converged = false;         // 是否收敛
    double likelihood_ratio = 0.0;  // ALT vs NULL 似然比

    // 标记
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

    // 空间特征
    int32_t primary_pos = -1;       // Primary alignment 位置
    int32_t locus_pos = -1;         // 最佳候选 locus 位置
    double d_spatial = 0.0;         // 与候选位点的距离

    // 几何一致性
    bool geom_ok = false;          // 断点几何解释一致
    double geom_score = 0.0;        // 几何一致性分数 [0, 1]

    // 对齐特征
    double align_score = 0.0;       // 对齐分数
    double normalized_score = 0.0; // 归一化分数
    float identity = 0.0f;         // 相似度

    // Contig 支持
    bool contig_support = false;    // 是否与代表 contig 一致
    double contig_score = 0.0;     // Contig 支持度 [0, 1]

    // TE 家族
    int te_family_id = -1;          // TE 家族 ID
    int representative_family = -1; // 代表结构 TE 家族

    // 证据强度
    double total_evidence = 0.0;   // 总证据强度
    uint32_t evidence_mask = 0;     // bit0=SA, bit1=CLIP, bit2=INS, bit3=ALIGN
};

// ============================================================================
// Prior Calculators
// ============================================================================

class SpatialPriorCalculator {
public:
    explicit SpatialPriorCalculator(const GenotypeConfig& config);

    // 计算空间先验
    void calculate_prior(double distance, double& pi_alt, double& pi_ref, double& pi_null) const;

    // 批量计算
    std::vector<double> calculate_priors(const std::vector<ReadEvidence>& evidence) const;

private:
    GenotypeConfig config_;
};

class StructuralPriorCalculator {
public:
    explicit StructuralPriorCalculator(const GenotypeConfig& config);

    // 计算结构先验加成
    double calculate_family_bonus(int read_family, int rep_family) const;

    // 批量计算
    std::vector<double> calculate_bonuses(
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
        const std::vector<double>& priors,
        GenotypeResult& result) const;

    // E-step: 计算责任度
    void e_step(
        const std::vector<ReadEvidence>& evidence,
        const std::vector<double>& priors,
        std::vector<std::array<double, 3>>& responsibilities) const;

    // M-step: 更新参数
    void m_step(
        const std::vector<std::array<double, 3>>& responsibilities,
        double& pi_alt, double& pi_ref, double& pi_null) const;

    // 计算对数似然
    double calculate_log_likelihood(
        const std::vector<ReadEvidence>& evidence,
        double pi_alt, double pi_ref, double pi_null) const;

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
        const std::vector<ReadSketch>& reads) const;

    // 批量分型
    std::vector<GenotypeResult> genotype_batch(
        const std::vector<StructuralRepresentative>& representatives,
        const std::vector<std::vector<LocusEvidence>>& all_evidence,
        const std::vector<ReadSketch>& all_reads) const;

    // 从 LocusEvidence 构建 ReadEvidence
    std::vector<ReadEvidence> build_evidence(
        const StructuralRepresentative& representative,
        const std::vector<LocusEvidence>& locus_evidence,
        const std::vector<ReadSketch>& reads) const;

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

    // 构建 contig 支持特征
    void compute_contig_features(
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
