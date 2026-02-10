#include "genotyping.h"
#include "assembly.h"
#include "placeability.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <sstream>
#include <array>
#include <cassert>
#include <unordered_map>
#include <iomanip>    // [修正] std::setprecision, std::fixed
#include <limits>     // [修正] std::numeric_limits

namespace placer {

// ============================================================================
// GenotypeResult::to_string()
// ============================================================================

std::string GenotypeResult::to_string() const {
    std::ostringstream oss;
    oss << "GT=" << genotype;

    oss << ";AF=" << std::fixed << std::setprecision(4) << af;

    oss << ";MIX_ALT=" << std::setprecision(4) << mix_alt;
    oss << ";MIX_REF=" << std::setprecision(4) << mix_ref;
    oss << ";MIX_NULL=" << std::setprecision(4) << mix_null;

    oss << ";GQ=" << gq;

    oss << ";GL_00=" << std::setprecision(2) << gt_log_lik[0];
    oss << ";GL_01=" << std::setprecision(2) << gt_log_lik[1];
    oss << ";GL_11=" << std::setprecision(2) << gt_log_lik[2];

    oss << ";AF_CI=[" << std::setprecision(4)
        << af_ci_low << "," << af_ci_high << "]";

    oss << ";ITER=" << em_iterations;
    oss << ";CONVERGED=" << (converged ? "1" : "0");

    if (high_background) oss << ";HIGH_BG=1";
    if (low_complexity) oss << ";LOW_COMPLEXITY=1";

    return oss.str();
}

// ============================================================================
// SpatialPriorCalculator
// ============================================================================

SpatialPriorCalculator::SpatialPriorCalculator(
    const GenotypeConfig& config)
    : config_(config) {}

void SpatialPriorCalculator::calculate_prior(
    double distance,
    double structural_score,
    double& pi_alt, double& pi_ref, double& pi_null) const {

    double pi_inlocus = std::exp(-distance / config_.spatial_lambda);

    double s_i = std::clamp(structural_score, 0.01, 0.99);

    pi_alt = pi_inlocus * s_i;
    pi_ref = pi_inlocus * (1.0 - s_i);
    pi_null = (1.0 - pi_inlocus) + config_.null_base_prior;

    double total = pi_alt + pi_ref + pi_null;
    if (total > 0) {
        pi_alt /= total;
        pi_ref /= total;
        pi_null /= total;
    }
}

std::vector<std::array<double, 3>>
SpatialPriorCalculator::calculate_priors(
    const std::vector<ReadEvidence>& evidence) const {

    std::vector<std::array<double, 3>> priors;
    priors.reserve(evidence.size());

    for (const auto& e : evidence) {
        double pi_alt, pi_ref, pi_null;
        calculate_prior(e.d_spatial, e.structural_score,
                        pi_alt, pi_ref, pi_null);
        priors.push_back({pi_alt, pi_ref, pi_null});
    }

    return priors;
}

// ============================================================================
// StructuralPriorCalculator
// ============================================================================

StructuralPriorCalculator::StructuralPriorCalculator(
    const GenotypeConfig& config)
    : config_(config) {}

double StructuralPriorCalculator::calculate_family_bonus(
    int read_family, int rep_family) const {

    if (read_family < 0 || rep_family < 0) {
        return 0.5;
    }

    if (read_family == rep_family) {
        return config_.family_match_bonus;
    } else {
        return config_.family_mismatch_penalty;
    }
}

double StructuralPriorCalculator::calculate_structural_score(
    const ReadEvidence& e,
    int representative_family) const {

    double geom_component = e.geom_ok ?
        (0.3 + 0.4 * e.geom_score) : 0.1;

    double family_component = calculate_family_bonus(
        e.te_family_id, representative_family);
    family_component = std::clamp(family_component, 0.0, 1.0);

    double contig_component = e.contig_support ?
        (0.3 + 0.5 * e.contig_score) : 0.1;

    double score = geom_component * 0.4 +
                   family_component * 0.3 +
                   contig_component * 0.3;

    return std::clamp(score, 0.01, 0.99);
}

std::vector<double> StructuralPriorCalculator::calculate_structural_scores(
    const std::vector<ReadEvidence>& evidence,
    int representative_family) const {

    std::vector<double> scores;
    scores.reserve(evidence.size());

    for (const auto& e : evidence) {
        scores.push_back(
            calculate_structural_score(e, representative_family));
    }

    return scores;
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

    double geom_component = e.geom_score;

    double contig_component = e.contig_support ?
        e.contig_score : 0.1;

    double align_component = e.normalized_score;

    double family_factor = 1.0;
    if (e.family_mismatch) {
        family_factor = config_.family_mismatch_penalty;
    } else if (e.family_match) {
        family_factor = config_.family_match_bonus;
    }

    double likelihood = geom_component *
                        contig_component *
                        align_component *
                        family_factor;

    return std::max(likelihood, config_.min_likelihood_ratio);
}

double EMEngine::likelihood_ref(const ReadEvidence& e) const {
    if (e.ref_span_ok) {
        double span_component = std::max(e.ref_span_score, 0.1);
        double align_component = e.normalized_score;

        double likelihood = span_component * 0.7 +
                            align_component * 0.3;

        return std::max(likelihood, config_.min_likelihood_ratio);
    }

    if (!e.geom_ok && e.normalized_score > 0.8) {
        return e.normalized_score * 0.3;
    }

    return config_.min_likelihood_ratio;
}

double EMEngine::likelihood_null(const ReadEvidence& e) const {
    double distance_component = 1.0 -
        std::exp(-e.d_spatial / config_.spatial_lambda);

    double geom_component = e.geom_ok ?
        (1.0 - e.geom_score) * 0.5 : 0.8;

    double align_component = 1.0 - e.normalized_score;

    double span_component = e.ref_span_ok ? 0.2 : 0.6;

    double likelihood = distance_component * 0.3 +
                        geom_component * 0.3 +
                        align_component * 0.2 +
                        span_component * 0.2;

    return std::max(likelihood, config_.min_likelihood_ratio);
}

// ============================================================================
// [关键修正] E-step：加入 mixture weights pi_k
//
// 标准 gated mixture EM:
//   r_{ik} ∝ pi_k * p(x_i | k) * prior_{ik}
//
// 之前缺少 pi_k，导致 responsibilities 与 mixture weights 无关，
// EM 实质上退化为单次软分类，不是在最大化声明的 log-likelihood。
// ============================================================================

void EMEngine::e_step(
    const std::vector<ReadEvidence>& evidence,
    const std::vector<std::array<double, 3>>& priors,
    double pi_alt, double pi_ref, double pi_null,
    std::vector<std::array<double, 3>>& responsibilities) const {

    responsibilities.resize(evidence.size());

    // prior_{ik} 下限，防止某类被完全打死导致 EM 不可逆
    double prior_min = config_.prior_min_bound;

    for (size_t i = 0; i < evidence.size(); ++i) {
        const auto& e = evidence[i];

        // [修正] prior_{ik} 加下限
        double prior_alt  = std::max((i < priors.size()) ? priors[i][0] : 1.0 / 3.0, prior_min);
        double prior_ref  = std::max((i < priors.size()) ? priors[i][1] : 1.0 / 3.0, prior_min);
        double prior_null = std::max((i < priors.size()) ? priors[i][2] : 1.0 / 3.0, prior_min);

        // [关键修正] r_{ik} ∝ pi_k * likelihood_k(x_i) * prior_{ik}
        double l_alt  = pi_alt  * likelihood_alt(e)  * prior_alt;
        double l_ref  = pi_ref  * likelihood_ref(e)  * prior_ref;
        double l_null = pi_null * likelihood_null(e)  * prior_null;

        double total = l_alt + l_ref + l_null;
        if (total > 0) {
            responsibilities[i][0] = l_alt / total;
            responsibilities[i][1] = l_ref / total;
            responsibilities[i][2] = l_null / total;
        } else {
            responsibilities[i] = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
        }
    }
}

void EMEngine::m_step(
    const std::vector<std::array<double, 3>>& responsibilities,
    double& pi_alt, double& pi_ref, double& pi_null) const {

    double sum_alt = 0.0, sum_ref = 0.0, sum_null = 0.0;

    for (const auto& r : responsibilities) {
        sum_alt  += r[0];
        sum_ref  += r[1];
        sum_null += r[2];
    }

    double total = sum_alt + sum_ref + sum_null;
    if (total > 0) {
        pi_alt  = sum_alt  / total;
        pi_ref  = sum_ref  / total;
        pi_null = sum_null / total;
    }
}

// ============================================================================
// [修正] log-likelihood 与 E-step 使用完全一致的公式
//
// log L = Σ_i log( pi_alt * L_alt(x_i) * prior_{i,alt}
//                + pi_ref * L_ref(x_i) * prior_{i,ref}
//                + pi_null * L_null(x_i) * prior_{i,null} )
// ============================================================================

double EMEngine::calculate_log_likelihood(
    const std::vector<ReadEvidence>& evidence,
    const std::vector<std::array<double, 3>>& priors,
    double pi_alt, double pi_ref, double pi_null) const {

    double log_lik = 0.0;

    // prior_{ik} 下限
    double prior_min = config_.prior_min_bound;

    for (size_t i = 0; i < evidence.size(); ++i) {
        const auto& e = evidence[i];

        // [修正] prior_{ik} 加下限（与 E-step 一致）
        double prior_alt  = std::max((i < priors.size()) ? priors[i][0] : 1.0 / 3.0, prior_min);
        double prior_ref  = std::max((i < priors.size()) ? priors[i][1] : 1.0 / 3.0, prior_min);
        double prior_null = std::max((i < priors.size()) ? priors[i][2] : 1.0 / 3.0, prior_min);

        // 与 E-step 完全一致
        double l_alt  = pi_alt  * likelihood_alt(e)  * prior_alt;
        double l_ref  = pi_ref  * likelihood_ref(e)  * prior_ref;
        double l_null = pi_null * likelihood_null(e)  * prior_null;

        double mixture = l_alt + l_ref + l_null;

        // [修正] 更稳健的 mixture 保护
        double min_mix = config_.min_likelihood_ratio;
        if (mixture < min_mix) {
            mixture = min_mix;
        }
        log_lik += std::log(mixture);
    }

    return log_lik;
}

// ============================================================================
// [修正] CI 计算：使用 Wilson interval 替代纯正态近似
//
// Wilson interval 对小样本和极端 AF 更稳健：
//   p̃ = (X + z²/2) / (n + z²)
//   CI = p̃ ± z * sqrt(p̃(1-p̃) / (n + z²))
//
// 当 e_alt + e_ref 很小时，回退到 Jeffreys Beta(0.5, 0.5) 先验的
// 正态近似。
// ============================================================================

void EMEngine::calculate_confidence_interval(
    double e_alt, double e_ref, double level,
    double& ci_low, double& ci_high) const {

    double n = e_alt + e_ref;

    if (n <= 0.0) {
        ci_low = 0.0;
        ci_high = 1.0;
        return;
    }

    // z 值查表
    double tail = (1.0 - level) / 2.0;
    double z;
    if      (tail <= 0.001)  z = 3.090;
    else if (tail <= 0.005)  z = 2.576;
    else if (tail <= 0.01)   z = 2.326;
    else if (tail <= 0.025)  z = 1.960;
    else if (tail <= 0.05)   z = 1.645;
    else                     z = 1.282;

    double z2 = z * z;

    if (n >= 5.0) {
        // Wilson interval
        double p_hat = e_alt / n;
        double denom = n + z2;
        double center = (e_alt + z2 / 2.0) / denom;
        double half_width = z * std::sqrt(
            (p_hat * (1.0 - p_hat) + z2 / (4.0 * n)) / denom);

        ci_low  = std::max(0.0, center - half_width);
        ci_high = std::min(1.0, center + half_width);
    } else {
        // 小样本：Jeffreys Beta posterior 正态近似
        // Beta(e_alt + 0.5, e_ref + 0.5)
        double alpha = e_alt + 0.5;
        double beta_param = e_ref + 0.5;
        double mean = alpha / (alpha + beta_param);
        double variance = (alpha * beta_param) /
            ((alpha + beta_param) * (alpha + beta_param) *
             (alpha + beta_param + 1.0));
        double std_dev = std::sqrt(variance);

        ci_low  = std::max(0.0, mean - z * std_dev);
        ci_high = std::min(1.0, mean + z * std_dev);
    }
}

// ============================================================================
// [修正] EM 主循环
//
// 修正点：
// 1. E-step 传入 pi_alt/pi_ref/pi_null
// 2. log-likelihood 与 E-step 公式一致
// 3. 收敛 break 前更新 prev_log_lik
// 4. e_null_sum 直接累加而非间接计算
// 5. LRT null 模型使用 calculate_log_likelihood(pi_alt=0,pi_ref=0,pi_null=1)
//    保证嵌套关系
// 6. CI 使用 Wilson interval
// ============================================================================

std::pair<double, int> EMEngine::run_em(
    const std::vector<ReadEvidence>& evidence,
    const std::vector<std::array<double, 3>>& priors,
    GenotypeResult& result) const {

    if (evidence.empty()) {
        result.mix_alt = 0.0;
        result.mix_ref = 0.0;
        result.mix_null = 1.0;
        return {0.0, 0};
    }

    // 初始化 mixture weights
    double pi_alt  = config_.em_init_alt_fraction;
    double pi_ref  = (1.0 - config_.em_init_alt_fraction) * 0.6;
    double pi_null = (1.0 - config_.em_init_alt_fraction) * 0.4;

    // 归一化
    {
        double init_total = pi_alt + pi_ref + pi_null;
        pi_alt  /= init_total;
        pi_ref  /= init_total;
        pi_null /= init_total;
    }

    double prev_log_lik = std::numeric_limits<double>::lowest();
    int iterations = 0;
    result.converged = false;

    std::vector<std::array<double, 3>> responsibilities;

    for (int iter = 0; iter < config_.max_em_iterations; ++iter) {
        iterations = iter + 1;

        // [修正] E-step：传入当前 pi_k
        e_step(evidence, priors, pi_alt, pi_ref, pi_null,
               responsibilities);

        // M-step
        m_step(responsibilities, pi_alt, pi_ref, pi_null);

        // 归一化（加微小值防止退化到 0）
        {
            pi_alt  = std::max(pi_alt,  1e-6);
            pi_ref  = std::max(pi_ref,  1e-6);
            pi_null = std::max(pi_null, 1e-6);
            double total = pi_alt + pi_ref + pi_null;
            pi_alt  /= total;
            pi_ref  /= total;
            pi_null /= total;
        }

        // log-likelihood（与 E-step 公式一致）
        double log_lik = calculate_log_likelihood(
            evidence, priors, pi_alt, pi_ref, pi_null);

        // 收敛检查：per-read 平均增益
        if (iter > 0) {
            double avg_gain = std::abs(log_lik - prev_log_lik)
                              / static_cast<double>(evidence.size());

            if (avg_gain < config_.em_convergence_threshold) {
                prev_log_lik = log_lik;
                result.converged = true;
                break;
            }
        }

        prev_log_lik = log_lik;
    }

    // 存储 mixture weights
    result.mix_alt  = pi_alt;
    result.mix_ref  = pi_ref;
    result.mix_null = pi_null;
    result.log_likelihood = prev_log_lik;
    result.em_iterations  = iterations;

    // [修正] 直接累加三个分量，避免间接计算的浮点误差
    double e_alt_sum = 0.0, e_ref_sum = 0.0, e_null_sum = 0.0;
    for (const auto& r : responsibilities) {
        e_alt_sum  += r[0];
        e_ref_sum  += r[1];
        e_null_sum += r[2];
    }
    result.e_alt_sum  = e_alt_sum;
    result.e_ref_sum  = e_ref_sum;
    result.e_null_sum = e_null_sum;

    // AF = E_alt / (E_alt + E_ref)，NULL 不进 AF
    if (e_alt_sum + e_ref_sum > 0) {
        result.af = e_alt_sum / (e_alt_sum + e_ref_sum);
    } else {
        result.af = 0.0;
    }

    // [修正] CI 使用 Wilson interval
    calculate_confidence_interval(
        e_alt_sum, e_ref_sum, config_.confidence_level,
        result.af_ci_low, result.af_ci_high);

    // [修正] LRT：null 模型 = full model 在 pi_null=1 下的特例
    // 这保证了嵌套关系，使 chi2(df=2) 近似有理论基础
    double null_log_lik = calculate_log_likelihood(
        evidence, priors,
        0.0,   // pi_alt  = 0
        0.0,   // pi_ref  = 0
        1.0    // pi_null = 1
    );

    result.lrt_statistic = 2.0 * (prev_log_lik - null_log_lik);
    result.likelihood_ratio = std::exp(
        std::min(prev_log_lik - null_log_lik, 700.0));

    return {prev_log_lik, iterations};
}

// ============================================================================
// Genotype-level Layer
// [修正] 将 eps 和 pop_af 从 config 读取，不再硬编码
// ============================================================================

void EMEngine::calculate_genotype_likelihoods(
    double e_alt_sum, double e_ref_sum, int n_reads,
    std::array<double, 3>& gt_log_lik) const {

    if (n_reads == 0) {
        gt_log_lik = {0.0, 0.0, 0.0};
        return;
    }

    // [修正] 从 config 读取，并 clamp 防止 log(0)
    double eps = std::clamp(config_.genotype_error_rate, 1e-10, 0.5);
    double pop_af = std::clamp(config_.population_af, 1e-10, 1.0 - 1e-10);

    // 三种 genotype 的预期 AF
    const double expected_af[3] = {eps, 0.5, 1.0 - eps};

    // log P(data | GT=g) ∝ E_alt * log(af_g) + E_ref * log(1 - af_g)
    for (int g = 0; g < 3; ++g) {
        double af_g = expected_af[g];
        gt_log_lik[g] = e_alt_sum * std::log(af_g) +
                        e_ref_sum * std::log(1.0 - af_g);
    }

    // Hardy-Weinberg genotype prior
    double prior_00 = (1.0 - pop_af) * (1.0 - pop_af);
    double prior_01 = 2.0 * pop_af * (1.0 - pop_af);
    double prior_11 = pop_af * pop_af;

    gt_log_lik[0] += std::log(prior_00);
    gt_log_lik[1] += std::log(prior_01);
    gt_log_lik[2] += std::log(prior_11);
}

int EMEngine::calculate_gq_from_gt(
    const std::array<double, 3>& gt_log_lik,
    std::string& best_gt) const {

    // 找最大 log-likelihood
    int best_idx = 0;
    double best_ll = gt_log_lik[0];
    for (int i = 1; i < 3; ++i) {
        if (gt_log_lik[i] > best_ll) {
            best_ll = gt_log_lik[i];
            best_idx = i;
        }
    }

    static const char* gt_strings[] = {"0/0", "0/1", "1/1"};
    best_gt = gt_strings[best_idx];

    // log-sum-exp 转 posterior
    double max_ll = best_ll;
    double sum_exp = 0.0;
    for (int i = 0; i < 3; ++i) {
        sum_exp += std::exp(gt_log_lik[i] - max_ll);
    }
    double log_total = max_ll + std::log(sum_exp);

    double p_best = std::exp(best_ll - log_total);

    // GQ = -10 * log10(1 - p_best)
    double p_error = std::max(1.0 - p_best, 1e-10);

    int gq = static_cast<int>(-10.0 * std::log10(p_error));
    return std::clamp(gq, 0, 99);
}

// ============================================================================
// Genotyper
// ============================================================================

Genotyper::Genotyper(const GenotypeConfig& config)
    : config_(config),
      spatial_prior_(config),
      structural_prior_(config),
      em_engine_(config) {}

// ============================================================================
// Feature Computation
// ============================================================================

void Genotyper::compute_spatial_features(
    const StructuralRepresentative& rep,
    const std::vector<LocusEvidence>& evidence,
    std::vector<ReadEvidence>& out_evidence) const {

    int32_t candidate_pos = rep.fingerprint.breakpoint_l;

    for (size_t i = 0; i < evidence.size() &&
                        i < out_evidence.size(); ++i) {
        const auto& e = evidence[i];
        auto& re = out_evidence[i];

        re.primary_pos = e.locus_pos;
        re.locus_pos = (candidate_pos > 0) ?
            candidate_pos : e.locus_pos;

        re.d_spatial = std::abs(
            static_cast<double>(e.locus_pos - re.locus_pos));
    }
}

void Genotyper::compute_geom_features(
    const StructuralRepresentative& rep,
    const std::vector<LocusEvidence>& evidence,
    std::vector<ReadEvidence>& out_evidence) const {

    for (size_t i = 0; i < evidence.size() &&
                        i < out_evidence.size(); ++i) {
        const auto& e = evidence[i];
        auto& re = out_evidence[i];

        bool has_seq_match = (e.evidence_bits & 0x01) != 0;
        bool has_distance  = (e.evidence_bits & 0x02) != 0;
        bool has_sa        = (e.evidence_bits & 0x04) != 0;
        bool has_indel     = (e.evidence_bits & 0x08) != 0;

        if ((has_seq_match || has_sa) && has_indel) {
            re.geom_ok = true;
            re.geom_score = 0.95;
        } else if (has_seq_match || has_sa) {
            re.geom_ok = true;
            re.geom_score = 0.8;
        } else if (has_distance && e.normalized_score > 0.7) {
            re.geom_ok = true;
            re.geom_score = 0.6;
        } else if (e.normalized_score > 0.85) {
            re.geom_ok = true;
            re.geom_score = 0.5;
        } else {
            re.geom_ok = false;
            re.geom_score = 0.2;
        }
    }
}

void Genotyper::compute_align_features(
    const std::vector<LocusEvidence>& evidence,
    std::vector<ReadEvidence>& out_evidence) const {

    for (size_t i = 0; i < evidence.size() &&
                        i < out_evidence.size(); ++i) {
        const auto& e = evidence[i];
        auto& re = out_evidence[i];

        re.align_score = e.total_score;
        re.normalized_score = e.normalized_score;
        re.alignment_identity = e.weighted_identity;
        re.identity = e.weighted_identity;  // 别名兼容
    }
}

void Genotyper::compute_contig_features(
    const StructuralRepresentative& rep,
    const std::vector<LocusEvidence>& evidence,
    const std::vector<ReadSketch>& reads,
    std::vector<ReadEvidence>& out_evidence) const {

    for (size_t i = 0; i < evidence.size() &&
                        i < out_evidence.size(); ++i) {
        const auto& e = evidence[i];
        auto& re = out_evidence[i];

        bool score_ok = (e.normalized_score > 0.5);
        bool identity_ok = (e.weighted_identity > 0.7f);

        re.contig_support = score_ok && identity_ok;
        re.contig_score = std::max(
            0.0, e.normalized_score * 0.6 +
                 e.weighted_identity * 0.4);
    }
}

// ============================================================================
// [修正] compute_ref_span_features
// SA tag 时直接将 ref_span_score 归零，避免残留非零值造成歧义
// ============================================================================

void Genotyper::compute_ref_span_features(
    const StructuralRepresentative& rep,
    const std::vector<LocusEvidence>& evidence,
    const std::vector<ReadSketch>& reads,
    const GenomeAccessor& genome,
    std::vector<ReadEvidence>& out_evidence) const {

    int32_t bp_left = rep.fingerprint.breakpoint_l;
    int32_t bp_right = rep.fingerprint.breakpoint_r;

    if (bp_left <= 0 || bp_right <= 0 || bp_left >= bp_right) {
        for (auto& re : out_evidence) {
            re.ref_span_ok = false;
            re.ref_span_score = 0.0;
        }
        return;
    }

    for (size_t i = 0; i < evidence.size() &&
                        i < out_evidence.size(); ++i) {
        const auto& e = evidence[i];
        auto& re = out_evidence[i];

        re.ref_span_ok = false;
        re.ref_span_score = 0.0;

        if (e.read_idx >= reads.size()) continue;
        const auto& read = reads[e.read_idx];

        int32_t read_start = read.pos;
        int32_t read_end = read.end_pos;

        // read 是否跨越断点
        bool spans_left  = (read_start < bp_left  && read_end > bp_left);
        bool spans_right = (read_start < bp_right && read_end > bp_right);
        bool spans_both  = (read_start < bp_left  && read_end > bp_right);

        if (!spans_left && !spans_right) {
            re.ref_span_ok = false;
            re.ref_span_score = 0.0;
            continue;
        }

        // 检查 CIGAR 中是否有大 clip / 大 indel
        bool has_large_clip = false;
        bool has_large_indel = false;

        for (const auto& [op, len] : read.cigar_ops) {
            if ((op == 'S' || op == 'H') && len >= 20) {
                has_large_clip = true;
            }
            if ((op == 'I' || op == 'D') && len >= 50) {
                has_large_indel = true;
            }
        }

        // 跨越两个断点且无大 clip/indel → 强 REF 证据
        if (spans_both && !has_large_clip && !has_large_indel) {
            re.ref_span_ok = true;

            double span_coverage = 1.0;
            double align_quality = e.normalized_score;
            double mapq_factor = std::min(
                static_cast<double>(read.mapq) / 60.0, 1.0);

            re.ref_span_score = span_coverage * 0.4 +
                                align_quality * 0.4 +
                                mapq_factor * 0.2;
        }
        else if ((spans_left || spans_right) &&
                 !has_large_clip && !has_large_indel) {
            // 只跨越一个断点
            re.ref_span_ok = true;
            double align_quality = e.normalized_score;
            double mapq_factor = std::min(
                static_cast<double>(read.mapq) / 60.0, 1.0);

            re.ref_span_score = (0.3 +
                                 align_quality * 0.4 +
                                 mapq_factor * 0.2) * 0.7;
        }
        else {
            re.ref_span_ok = false;
            re.ref_span_score = 0.0;
        }

        // [修正] SA tag → 参考上不连续，强烈反对 REF
        // 直接归零，不保留残留值
        if (read.has_sa) {
            re.ref_span_ok = false;
            re.ref_span_score = 0.0;
        }
    }
}

void Genotyper::compute_family_features(
    const StructuralRepresentative& rep,
    const std::vector<LocusEvidence>& evidence,
    const std::vector<ReadSketch>& reads,
    std::vector<ReadEvidence>& out_evidence) const {

    int rep_family = rep.fingerprint.te_family_id;

    for (size_t i = 0; i < evidence.size() &&
                        i < out_evidence.size(); ++i) {
        auto& re = out_evidence[i];

        // ReadSketch 没有 te_family_id，设为 -1 (未知)
        re.te_family_id = -1;

        // 如果代表有家族信息，设为不匹配（因为无法确认 read 的家族）
        if (rep_family >= 0) {
            re.family_match = false;
            re.family_mismatch = false;  // 未知情况不设为 mismatch
        } else {
            re.family_match = false;
            re.family_mismatch = false;
        }
    }
}

std::vector<ReadEvidence> Genotyper::build_evidence(
    const StructuralRepresentative& representative,
    const std::vector<LocusEvidence>& locus_evidence,
    const std::vector<ReadSketch>& reads,
    const GenomeAccessor& genome) const {

    std::vector<ReadEvidence> evidence;
    evidence.reserve(locus_evidence.size());

    for (size_t i = 0; i < locus_evidence.size(); ++i) {
        ReadEvidence re;
        re.read_idx = locus_evidence[i].read_idx;
        evidence.push_back(re);
    }

    compute_spatial_features(representative, locus_evidence, evidence);
    compute_geom_features(representative, locus_evidence, evidence);
    compute_align_features(locus_evidence, evidence);
    compute_contig_features(representative, locus_evidence,
                            reads, evidence);
    compute_ref_span_features(representative, locus_evidence,
                              reads, genome, evidence);
    compute_family_features(representative, locus_evidence,
                            reads, evidence);

    int rep_family = representative.fingerprint.te_family_id;
    for (auto& re : evidence) {
        re.structural_score =
            structural_prior_.calculate_structural_score(
                re, rep_family);
    }

    return evidence;
}

// ============================================================================
// Main Genotyping Entry Point
// ============================================================================

GenotypeResult Genotyper::genotype(
    const StructuralRepresentative& representative,
    const std::vector<LocusEvidence>& evidence,
    const std::vector<ReadSketch>& reads,
    const GenomeAccessor& genome) const {

    GenotypeResult result;

    if (evidence.size() < static_cast<size_t>(
            config_.min_reads_for_genotyping)) {
        result.genotype = "./.";
        result.gq = 0;
        result.mix_null = 1.0;
        return result;
    }

    // ===== Layer 1: Read-level 分类 =====

    auto read_evidence = build_evidence(
        representative, evidence, reads, genome);

    auto priors = spatial_prior_.calculate_priors(read_evidence);

    auto [log_lik, iterations] = em_engine_.run_em(
        read_evidence, priors, result);

    // ===== Layer 2: Genotype-level =====

    std::array<double, 3> gt_log_lik;
    em_engine_.calculate_genotype_likelihoods(
        result.e_alt_sum, result.e_ref_sum,
        static_cast<int>(evidence.size()),
        gt_log_lik);

    result.gt_log_lik = gt_log_lik;

    std::string best_gt;
    result.gq = em_engine_.calculate_gq_from_gt(
        gt_log_lik, best_gt);
    result.genotype = best_gt;

    // NULL 占比过高
    if (result.mix_null > 0.7 && result.mix_alt < 0.15) {
        result.high_background = true;
        if (result.gq < 5) {
            result.genotype = "./.";
        }
    }

    // 低复杂度检测
    if (representative.fingerprint.trunc_level > 0 ||
        representative.rep_sequence.size() < 50) {
        result.low_complexity = true;
    }

    // LRT 显著性检查 (df=2, 默认 chi2 临界值 5.99 at p=0.05)
    // [修正] 从 config 读取阈值
    if (result.lrt_statistic < config_.lrt_threshold) {
        result.lrt_significant = false;
        result.gq = std::min(result.gq, 10);
    } else {
        result.lrt_significant = true;
    }

    return result;
}

// ============================================================================
// Batch Genotyping
// ============================================================================

std::vector<GenotypeResult> Genotyper::genotype_batch(
    const std::vector<StructuralRepresentative>& representatives,
    const std::vector<std::vector<LocusEvidence>>& all_evidence,
    const std::vector<ReadSketch>& all_reads,
    const GenomeAccessor& genome) const {

    std::vector<GenotypeResult> results;
    results.reserve(representatives.size());

    for (size_t i = 0; i < representatives.size(); ++i) {
        const auto& rep = representatives[i];
        const auto& ev = (i < all_evidence.size()) ?
            all_evidence[i] : std::vector<LocusEvidence>{};

        results.push_back(genotype(rep, ev, all_reads, genome));
    }

    return results;
}

// ============================================================================
// Utility Functions
// ============================================================================

double beta_pdf(double x, double alpha, double beta_param) {
    if (x <= 0.0 || x >= 1.0 || alpha <= 0.0 || beta_param <= 0.0) {
        return 0.0;
    }
    double log_pdf = (alpha - 1.0) * std::log(x) +
                     (beta_param - 1.0) * std::log(1.0 - x);
    return std::exp(log_pdf);
}

double beta_cdf_approx(double x, double alpha, double beta_param) {
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;
    if (alpha <= 0.0 || beta_param <= 0.0) return x;

    double mean = alpha / (alpha + beta_param);
    double variance = (alpha * beta_param) /
        ((alpha + beta_param) * (alpha + beta_param) *
         (alpha + beta_param + 1.0));
    double std_dev = std::sqrt(variance);

    if (std_dev <= 0.0) return (x >= mean) ? 1.0 : 0.0;

    double z = (x - mean) / std_dev;
    double t = 1.0 / (1.0 + 0.2316419 * std::abs(z));
    double d = 0.3989422804014327;
    double p = d * std::exp(-z * z / 2.0) *
        (t * (0.319381530 +
         t * (-0.356563782 +
         t * (1.781477937 +
         t * (-1.821255978 +
         t * 1.330274429)))));

    return (z >= 0) ? (1.0 - p) : p;
}

std::array<double, 2> beta_binomial_ci(
    int successes, int trials, double level) {

    double alpha = successes + 1.0;
    double beta_param = trials - successes + 1.0;

    double tail = (1.0 - level) / 2.0;

    auto find_quantile = [&](double target_cdf) -> double {
        double lo = 0.0, hi = 1.0;
        for (int iter = 0; iter < 50; ++iter) {
            double mid = (lo + hi) / 2.0;
            double cdf = beta_cdf_approx(mid, alpha, beta_param);
            if (cdf < target_cdf) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        return (lo + hi) / 2.0;
    };

    return {
        find_quantile(tail),
        find_quantile(1.0 - tail)
    };
}

double phred_to_prob(int phred) {
    return std::pow(10.0, -phred / 10.0);
}

int prob_to_phred(double prob) {
    if (prob >= 1.0) return 0;
    if (prob <= 0.0) return 99;
    return std::clamp(
        static_cast<int>(-10.0 * std::log10(prob)), 0, 99);
}

double likelihood_ratio_test(
    double log_lik_null, double log_lik_alt) {
    return 2.0 * (log_lik_alt - log_lik_null);
}

}  // namespace placer

