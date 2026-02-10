#include "genotyping.h"
#include "assembly.h"
#include "placeability.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <sstream>

namespace placer {

// ============================================================================
// GenotypeConfig
// ============================================================================

// ============================================================================
// GenotypeResult
// ============================================================================

std::string GenotypeResult::to_string() const {
    std::ostringstream oss;
    oss << "GT=" << genotype
        << ";AF=" << std::fixed << std::setprecision(4) << af
        << ";E[ALT]=" << e_alt
        << ";E[REF]=" << e_ref
        << ";E[NULL]=" << e_null
        << ";GQ=" << gq
        << ";AF_CI=[" << alt_ci_low << "," << alt_ci_high << "]"
        << ";ITER=" << em_iterations
        << ";CONVERGED=" << (converged ? "1" : "0");
    return oss.str();
}

// ============================================================================
// SpatialPriorCalculator
// ============================================================================

SpatialPriorCalculator::SpatialPriorCalculator(const GenotypeConfig& config)
    : config_(config) {}

void SpatialPriorCalculator::calculate_prior(
    double distance,
    double& pi_alt, double& pi_ref, double& pi_null) const {

    // 空间衰减
    double decay = std::exp(-distance / config_.spatial_lambda);
    double noise = config_.null_base_prior;

    // ALT 和 REF 有相似的空间先验（都依赖于接近位点）
    pi_alt = decay;
    pi_ref = decay;
    pi_null = 1.0 - decay + noise;

    // 归一化
    double total = pi_alt + pi_ref + pi_null;
    if (total > 0) {
        pi_alt /= total;
        pi_ref /= total;
        pi_null /= total;
    }
}

std::vector<double> SpatialPriorCalculator::calculate_priors(
    const std::vector<ReadEvidence>& evidence) const {

    std::vector<double> priors;
    priors.reserve(evidence.size());

    for (const auto& e : evidence) {
        double pi_alt, pi_ref, pi_null;
        calculate_prior(e.d_spatial, pi_alt, pi_ref, pi_null);
        // 返回 ALT 先验（用于 EM）
        priors.push_back(pi_alt);
    }

    return priors;
}

// ============================================================================
// StructuralPriorCalculator
// ============================================================================

StructuralPriorCalculator::StructuralPriorCalculator(const GenotypeConfig& config)
    : config_(config) {}

double StructuralPriorCalculator::calculate_family_bonus(
    int read_family, int rep_family) const {

    if (read_family < 0 || rep_family < 0) {
        return 1.0;  // 未知家族，无加成
    }

    if (read_family == rep_family) {
        return config_.family_match_bonus;
    } else {
        return config_.family_mismatch_penalty;
    }
}

std::vector<double> StructuralPriorCalculator::calculate_bonuses(
    const std::vector<ReadEvidence>& evidence,
    int representative_family) const {

    std::vector<double> bonuses;
    bonuses.reserve(evidence.size());

    for (const auto& e : evidence) {
        bonuses.push_back(calculate_family_bonus(e.te_family_id, representative_family));
    }

    return bonuses;
}

// ============================================================================
// EM Engine
// ============================================================================

EMEngine::EMEngine(const GenotypeConfig& config)
    : config_(config) {}

double EMEngine::likelihood_alt(const ReadEvidence& e) const {
    if (!e.geom_ok) {
        return config_.min_likelihood_ratio;
    }

    double geom_component = e.geom_score * config_.geom_ok_weight;
    double contig_component = e.contig_support ? e.contig_score * config_.contig_support_weight : 1.0;
    double align_component = e.normalized_score * config_.align_score_weight;

    double likelihood = geom_component * contig_component * align_component;
    return std::max(likelihood, config_.min_likelihood_ratio);
}

double EMEngine::likelihood_ref(const ReadEvidence& e) const {
    // REF: 不需要几何一致，但需要对齐分数好
    double align_component = e.normalized_score;
    return std::max(align_component, config_.min_likelihood_ratio);
}

double EMEngine::likelihood_null(const ReadEvidence& e) const {
    // NULL: 距离远、几何不一致、或不支持当前结构
    double distance_component = 1.0 - std::exp(-e.d_spatial / config_.spatial_lambda);
    double geom_component = e.geom_ok ? (1.0 - e.geom_score) : 1.0;

    double likelihood = distance_component * geom_component;
    return std::max(likelihood, config_.min_likelihood_ratio);
}

void EMEngine::e_step(
    const std::vector<ReadEvidence>& evidence,
    const std::vector<double>& priors,
    std::vector<std::array<double, 3>>& responsibilities) const {

    responsibilities.resize(evidence.size());

    for (size_t i = 0; i < evidence.size(); ++i) {
        const auto& e = evidence[i];
        double prior = (i < priors.size()) ? priors[i] : 0.5;

        double l_alt = likelihood_alt(e) * prior;
        double l_ref = likelihood_ref(e) * prior;
        double l_null = likelihood_null(e) * (1.0 - prior);

        double total = l_alt + l_ref + l_null;
        if (total > 0) {
            responsibilities[i][0] = l_alt / total;
            responsibilities[i][1] = l_ref / total;
            responsibilities[i][2] = l_null / total;
        } else {
            responsibilities[i] = {1.0/3.0, 1.0/3.0, 1.0/3.0};
        }
    }
}

void EMEngine::m_step(
    const std::vector<std::array<double, 3>>& responsibilities,
    double& pi_alt, double& pi_ref, double& pi_null) const {

    double sum_alt = 0.0, sum_ref = 0.0, sum_null = 0.0;
    double total_weight = 0.0;

    for (const auto& r : responsibilities) {
        sum_alt += r[0];
        sum_ref += r[1];
        sum_null += r[2];
        total_weight += 1.0;
    }

    if (total_weight > 0) {
        pi_alt = sum_alt / total_weight;
        pi_ref = sum_ref / total_weight;
        pi_null = sum_null / total_weight;
    }
}

double EMEngine::calculate_log_likelihood(
    const std::vector<ReadEvidence>& evidence,
    double pi_alt, double pi_ref, double pi_null) const {

    double log_lik = 0.0;

    for (const auto& e : evidence) {
        double l_alt = likelihood_alt(e);
        double l_ref = likelihood_ref(e);
        double l_null = likelihood_null(e);

        double likelihood = pi_alt * l_alt + pi_ref * l_ref + pi_null * l_null;
        if (likelihood > 0) {
            log_lik += std::log(likelihood);
        }
    }

    return log_lik;
}

void EMEngine::calculate_confidence_interval(
    double alpha, double beta, double level,
    double& ci_low, double& ci_high) const {

    // 简化的 Beta 分布 CI 计算
    double mode = (alpha > 1) ? (alpha - 1) / (alpha + beta - 2) : 0.5;
    double variance = (alpha * beta) / ((alpha + beta) * (alpha + beta) * (alpha + beta + 1));

    // 近似正态分布的 z 值
    double z = 1.96;  // 95% CI
    double std_dev = std::sqrt(variance);

    ci_low = std::max(0.0, mode - z * std_dev);
    ci_high = std::min(1.0, mode + z * std_dev);
}

std::pair<double, int> EMEngine::run_em(
    const std::vector<ReadEvidence>& evidence,
    const std::vector<double>& priors,
    GenotypeResult& result) const {

    if (evidence.empty()) {
        result.e_alt = 0.0;
        result.e_ref = 0.0;
        result.e_null = 1.0;
        return {0.0, 0};
    }

    // 初始化
    double pi_alt = config_.em_init_alt_fraction;
    double pi_ref = (1.0 - config_.em_init_alt_fraction) / 2.0;
    double pi_null = (1.0 - config_.em_init_alt_fraction) / 2.0;

    double prev_log_lik = std::numeric_limits<double>::lowest();
    int iterations = 0;

    std::vector<std::array<double, 3>> responsibilities;

    // EM 迭代
    for (int iter = 0; iter < config_.max_em_iterations; ++iter) {
        iterations = iter + 1;

        // E-step
        e_step(evidence, priors, responsibilities);

        // M-step
        m_step(responsibilities, pi_alt, pi_ref, pi_null);

        // 归一化
        double total = pi_alt + pi_ref + pi_null;
        if (total > 0) {
            pi_alt /= total;
            pi_ref /= total;
            pi_null /= total;
        }

        // 计算对数似然
        double log_lik = calculate_log_likelihood(evidence, pi_alt, pi_ref, pi_null);

        // 检查收敛
        if (iter > 0 && std::abs(log_lik - prev_log_lik) < config_.em_convergence_threshold) {
            result.converged = true;
            break;
        }

        prev_log_lik = log_lik;
    }

    result.e_alt = pi_alt;
    result.e_ref = pi_ref;
    result.e_null = pi_null;
    result.log_likelihood = prev_log_lik;
    result.em_iterations = iterations;

    // 计算置信区间（使用 Beta-Binomial 近似）
    int n = static_cast<int>(evidence.size());
    int k = static_cast<int>(pi_alt * n + 0.5);
    calculate_confidence_interval(k + 1, n - k + 1, config_.confidence_level,
                                  result.alt_ci_low, result.alt_ci_high);

    // 似然比 (ALT vs NULL)
    double null_log_lik = 0.0;
    for (const auto& e : evidence) {
        double l_null = likelihood_null(e);
        if (l_null > 0) {
            null_log_lik += std::log(l_null);
        }
    }
    result.likelihood_ratio = std::exp(result.log_likelihood - null_log_lik);

    return {prev_log_lik, iterations};
}

// ============================================================================
// Genotyper
// ============================================================================

Genotyper::Genotyper(const GenotypeConfig& config)
    : config_(config),
      spatial_prior_(config),
      structural_prior_(config),
      em_engine_(config) {}

void Genotyper::compute_spatial_features(
    const StructuralRepresentative& rep,
    const std::vector<LocusEvidence>& evidence,
    std::vector<ReadEvidence>& out_evidence) const {

    // 从 fingerprint 获取候选位点
    int32_t candidate_pos = rep.fingerprint.breakpoint_l;

    for (size_t i = 0; i < evidence.size() && i < out_evidence.size(); ++i) {
        const auto& e = evidence[i];
        auto& re = out_evidence[i];

        re.primary_pos = e.locus_pos;
        re.locus_pos = (candidate_pos > 0) ? candidate_pos : e.locus_pos;

        re.d_spatial = std::abs(static_cast<double>(e.locus_pos - re.locus_pos));
    }
}

void Genotyper::compute_geom_features(
    const StructuralRepresentative& rep,
    const std::vector<LocusEvidence>& evidence,
    std::vector<ReadEvidence>& out_evidence) const {

    for (size_t i = 0; i < evidence.size() && i < out_evidence.size(); ++i) {
        const auto& e = evidence[i];
        auto& re = out_evidence[i];

        // 几何一致性：取决于证据类型
        bool has_breakpoint = (e.evidence_bits & 0x01) || (e.evidence_bits & 0x02);
        bool has_high_identity = (e.evidence_bits & 0x08);

        if (has_breakpoint && has_high_identity) {
            re.geom_ok = true;
            re.geom_score = 0.9;
        } else if (has_breakpoint) {
            re.geom_ok = true;
            re.geom_score = 0.7;
        } else if (e.normalized_score > 0.8) {
            re.geom_ok = true;
            re.geom_score = 0.6;
        } else {
            re.geom_ok = false;
            re.geom_score = 0.3;
        }
    }
}

void Genotyper::compute_align_features(
    const std::vector<LocusEvidence>& evidence,
    std::vector<ReadEvidence>& out_evidence) const {

    for (size_t i = 0; i < evidence.size() && i < out_evidence.size(); ++i) {
        const auto& e = evidence[i];
        auto& re = out_evidence[i];

        re.align_score = e.total_score;
        re.normalized_score = e.normalized_score;
        re.identity = e.weighted_identity;
    }
}

void Genotyper::compute_contig_features(
    const StructuralRepresentative& rep,
    const std::vector<LocusEvidence>& evidence,
    const std::vector<ReadSketch>& reads,
    std::vector<ReadEvidence>& out_evidence) const {

    // 简化的 contig 支持计算
    for (size_t i = 0; i < evidence.size() && i < out_evidence.size(); ++i) {
        const auto& e = evidence[i];
        auto& re = out_evidence[i];

        // 如果 read 支持的结构与代表结构一致
        re.contig_support = (e.normalized_score > 0.6);
        re.contig_score = e.normalized_score;
    }
}

std::vector<ReadEvidence> Genotyper::build_evidence(
    const StructuralRepresentative& representative,
    const std::vector<LocusEvidence>& locus_evidence,
    const std::vector<ReadSketch>& reads) const {

    std::vector<ReadEvidence> evidence;
    evidence.reserve(locus_evidence.size());

    for (size_t i = 0; i < locus_evidence.size(); ++i) {
        ReadEvidence re;
        re.read_idx = i;
        evidence.push_back(re);
    }

    // 计算特征
    compute_spatial_features(representative, locus_evidence, evidence);
    compute_geom_features(representative, locus_evidence, evidence);
    compute_align_features(locus_evidence, evidence);
    compute_contig_features(representative, locus_evidence, reads, evidence);

    return evidence;
}

std::string Genotyper::determine_genotype(
    double e_alt, double e_ref, double e_null) const {

    // 判定逻辑
    if (e_null > 0.5 && e_alt < 0.2) {
        return "./.";  // 无法判定
    }

    double alt_rate = e_alt / (e_alt + e_ref + 1e-10);

    if (alt_rate < 0.1) {
        return "0/0";
    } else if (alt_rate < 0.5) {
        return "0/1";
    } else {
        return "1/1";
    }
}

int Genotyper::calculate_gq(double best_prob, double second_best_prob) const {
    double ratio = best_prob / (second_best_prob + 1e-10);
    if (ratio <= 1.0) return 0;

    int phred = static_cast<int>(-10.0 * std::log10(1.0 - 1.0/ratio));
    return std::min(phred, 99);
}

void Genotyper::calculate_af_ci(
    double e_alt, double e_ref, double e_null,
    GenotypeResult& result) const {

    double total = e_alt + e_ref + e_null;
    if (total <= 0) {
        result.af = 0.0;
        return;
    }

    result.af = e_alt / total;

    // Beta-Binomial CI
    int n = 100;  // 假设 100x 覆盖
    int k = static_cast<int>(result.af * n + 0.5);

    double alpha = k + 1;
    double beta = n - k + 1;

    double z = 1.96;  // 95% CI
    double mean = alpha / (alpha + beta);
    double variance = (alpha * beta) / ((alpha + beta) * (alpha + beta) * (alpha + beta + 1));

    result.af_ci_low = std::max(0.0, mean - z * std::sqrt(variance));
    result.af_ci_high = std::min(1.0, mean + z * std::sqrt(variance));
}

GenotypeResult Genotyper::genotype(
    const StructuralRepresentative& representative,
    const std::vector<LocusEvidence>& evidence,
    const std::vector<ReadSketch>& reads) const {

    GenotypeResult result;

    // 检查最小 reads 数
    if (evidence.size() < static_cast<size_t>(config_.min_reads_for_genotyping)) {
        result.genotype = "./.";
        result.gq = 0;
        result.e_null = 1.0;
        return result;
    }

    // 构建 ReadEvidence
    auto read_evidence = build_evidence(representative, evidence, reads);

    // 计算空间先验
    auto priors = spatial_prior_.calculate_priors(read_evidence);

    // 运行 EM
    auto [log_lik, iterations] = em_engine_.run_em(read_evidence, priors, result);

    // 判定基因型
    result.genotype = determine_genotype(result.e_alt, result.e_ref, result.e_null);

    // 计算 GQ
    double probs[] = {result.e_alt, result.e_ref, result.e_null};
    std::sort(probs, probs + 3, std::greater<double>());
    result.gq = calculate_gq(probs[0], probs[1]);

    // 计算 AF CI
    calculate_af_ci(result.e_alt, result.e_ref, result.e_null, result);

    // 设置标记
    result.high_background = (result.e_null > 0.5);

    // 结构简单检测
    if (representative.poly_summary.size() < 2) {
        result.low_complexity = true;
    }

    return result;
}

std::vector<GenotypeResult> Genotyper::genotype_batch(
    const std::vector<StructuralRepresentative>& representatives,
    const std::vector<std::vector<LocusEvidence>>& all_evidence,
    const std::vector<ReadSketch>& all_reads) const {

    std::vector<GenotypeResult> results;
    results.reserve(representatives.size());

    for (size_t i = 0; i < representatives.size(); ++i) {
        const auto& rep = representatives[i];
        const auto& evidence = (i < all_evidence.size()) ? all_evidence[i] : std::vector<LocusEvidence>{};

        results.push_back(genotype(rep, evidence, all_reads));
    }

    return results;
}

// ============================================================================
// Utility Functions
// ============================================================================

double beta_pdf(double x, double alpha, double beta) {
    if (x <= 0 || x >= 1) return 0.0;
    // 简化的 PDF（实际应该用 Gamma 函数）
    double log_pdf = (alpha - 1) * std::log(x) + (beta - 1) * std::log(1 - x);
    return std::exp(log_pdf);
}

double beta_cdf(double x, double alpha, double beta) {
    // 简化的 CDF（实际应该用不完全 Beta 函数）
    if (x <= 0) return 0.0;
    if (x >= 1) return 1.0;
    return x;  // 简化：假设均匀先验
}

std::array<double, 2> beta_binomial_ci(int successes, int trials, double level) {
    double alpha = successes + 1;
    double beta_param = trials - successes + 1;
    // 近似正态分布的 z 值 (95% CI = 1.96)
    double z = 1.96;

    double mean = alpha / (alpha + beta_param);
    double variance = (alpha * beta_param) / ((alpha + beta_param) * (alpha + beta_param) * (alpha + beta_param + 1));
    double std_dev = std::sqrt(variance);

    return {
        std::max(0.0, mean - z * std_dev),
        std::min(1.0, mean + z * std_dev)
    };
}

double phred_to_prob(int phred) {
    return std::pow(10.0, -phred / 10.0);
}

int prob_to_phred(double prob) {
    if (prob >= 1.0) return 0;
    if (prob <= 0) return 99;
    return static_cast<int>(-10.0 * std::log10(prob));
}

double likelihood_ratio_test(double log_lik_null, double log_lik_alt) {
    return 2.0 * (log_lik_alt - log_lik_null);
}

}  // namespace placer
