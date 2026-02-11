/**
 * PLACER - Complete TE Insertion Detection Pipeline
 *
 * Phases:
 * 1. Stream Layer (BAM streaming, window buffering, triggering)
 * 2. Gate 1 (TE-proxy filtering)
 * 3. Component Building (anchor extraction, clustering)
 * 4. Local Realignment (restricted search, evidence collection)
 * 5. Assembly (Graph-POA, structural collapsing)
 * 6. Placeability Scoring (Δ score, side consistency, tier assignment)
 * 7. Genotyping (EM with spatial priors)
 * 8. TE Reverse Index (recall for missing secondary alignments)
 *
 * Output: scientific.vcf with tier assignments and genotyping
 */

#include <iostream>
#include <fstream>
#include <chrono>
#include <memory>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <sstream>
#include <iomanip>

#include "bam_reader.h"
#include "window_buffer.h"
#include "window_stats.h"
#include "trigger.h"
#include "task_queue.h"
#include "stream_processor.h"
#include "gate1.h"
#include "gate1_filter.h"
#include "component_builder.h"
#include "local_realign.h"
#include "assembly.h"
#include "placeability.h"
#include "genotyping.h"
#include "te_reverse_index.h"

namespace placer {

// ============================================================================
// Pipeline Configuration
// ============================================================================

struct PipelineConfig {
    // Stream Layer
    int window_size = 10000;
    int window_step = 5000;
    int max_priority_reads = 50;
    int max_normal_reads = 200;
    int min_clip_bp = 50;
    int min_sa_reads = 3;

    // Gate 1
    bool gate1_enabled = true;
    std::string te_fasta_path;
    int gate1_max_reads = 100;
    double gate1_sample_ratio = 1.0;

    // Component Builder
    double component_min_density = 0.01;
    int component_breaker_threshold = 200;
    int component_max_recursive_depth = 4;

    // Local Realignment
    int flank_length = 1000;
    int32_t search_window = 10000;
    double min_realign_score = 0.3;
    int max_locus_per_component = 20;

    // Assembly
    int assembly_min_reads = 3;
    int assembly_max_reads = 50;
    int assembly_max_paths = 2;

    // Placeability
    double delta_tier1 = 30.0;
    double delta_tier2 = 10.0;
    int side_consistency_gap = 50;
    int min_locus_support = 2;

    // Genotyping
    int max_em_iterations = 20;
    double em_convergence = 1e-6;
    double spatial_lambda = 100.0;

    // TE Reverse Index
    bool te_reverse_enabled = false;
    std::string genome_fasta_path;
};

// ============================================================================
// Pipeline Output
// ============================================================================

struct PipelineResult {
    int64_t total_reads = 0;
    int64_t gate1_passed = 0;
    int64_t triggered_windows = 0;
    int64_t built_components = 0;
    int64_t assembled_contigs = 0;
    int64_t tier1_loci = 0;
    int64_t tier2_loci = 0;
    int64_t tier3_loci = 0;
    int64_t genotyped = 0;
    int64_t rescued_loci = 0;

    std::vector<StructuralRepresentative> representatives;
    std::vector<GenotypeResult> genotypes;
    std::vector<ExtendedPlaceabilityReport> placeability_reports;

    std::string vcf_path;
};

// ============================================================================
// Per-locus result bundle: ties together all phases for one call site
// ============================================================================

struct LocusCallBundle {
    // Source
    int32_t component_id = -1;
    size_t representative_idx = SIZE_MAX;

    // Coordinates (from assembly/realignment)
    int32_t chrom_tid = -1;
    int32_t position = 0;
    std::string chrom_name;

    // Assembly
    StructuralRepresentative representative;
    bool has_representative = false;

    // Evidence (from real local realignment)
    std::vector<LocusEvidence> evidence;

    // Placeability
    ExtendedPlaceabilityReport placeability;
    bool has_placeability = false;

    // Genotyping
    GenotypeResult genotype;
    bool has_genotype = false;

    // Flags
    bool from_rescue = false;
};

// ============================================================================
// VCF Writer
//
// [修正] write_record 输出 FORMAT + SAMPLE 两列，保证 VCF 列数正确
// [修正] CHROM 由调用方传入真实 contig 名，不再用 "chr" + index
// ============================================================================

class VCFWriter {
public:
    explicit VCFWriter(const std::string& path,
                       const std::string& sample_name = "SAMPLE")
        : ofs_(path), sample_name_(sample_name) {
        if (!ofs_.is_open()) {
            throw std::runtime_error("Cannot open VCF for writing: " + path);
        }
        write_header();
    }

    void write_record(const LocusCallBundle& bundle) {
        if (!bundle.has_placeability) return;

        const auto& report = bundle.placeability;
        const auto& gt = bundle.genotype;

        // CHROM
        ofs_ << bundle.chrom_name << "\t";

        // POS
        ofs_ << bundle.position << "\t";

        // ID
        ofs_ << "PLACER_" << bundle.component_id << "\t";

        // REF / ALT
        ofs_ << "N\t<INS:TE>\t";

        // QUAL
        ofs_ << std::fixed << std::setprecision(1) << report.delta_score << "\t";

        // FILTER
        std::string filter = "PASS";
        if (report.tier == Tier::TIER3) filter = "LowQual";
        if (bundle.has_genotype && gt.gq < 5) filter = "LowGQ";
        ofs_ << filter << "\t";

        // INFO
        ofs_ << "TIER=" << tier_to_string(report.tier);
        ofs_ << ";PLACEABILITY=" << std::setprecision(2) << report.delta_score;
        ofs_ << ";SUPPORT_READS=" << report.support_reads;

        if (bundle.has_representative) {
            const auto& fp = bundle.representative.fingerprint;
            ofs_ << ";SVLEN=" << (fp.ins_length_max > 0 ?
                fp.ins_length_max : 0);

            if (fp.te_family_id >= 0) {
                ofs_ << ";TE_FAMILY=" << fp.te_family_id;
            }

            if (fp.orientation != 0) {
                ofs_ << ";TE_ORIENT=" << (fp.orientation > 0 ? "+" : "-");
            }

            if (fp.breakpoint_l > 0 && fp.breakpoint_r > 0) {
                ofs_ << ";BP_LEFT=" << fp.breakpoint_l;
                ofs_ << ";BP_RIGHT=" << fp.breakpoint_r;
                // [修正 P0.5] TSD 不能简单用断点差计算
                // 只有当 fingerprint 明确提供 tsd_len 时才写，否则留空
                // bp_right - bp_left 是插入跨度，不是 TSD 长度
            }
        }

        if (bundle.has_genotype) {
            ofs_ << ";AF=" << std::setprecision(4) << gt.af;
            ofs_ << ";AF_CI=" << gt.af_ci_low << "," << gt.af_ci_high;
            ofs_ << ";GQ=" << gt.gq;
            ofs_ << ";EALT=" << std::setprecision(3) << gt.mix_alt;
            ofs_ << ";EREF=" << gt.mix_ref;
            ofs_ << ";ENULL=" << gt.mix_null;

            // Flags
            std::string flags;
            if (gt.high_background) flags += "HIGH_BACKGROUND,";
            if (gt.low_complexity) flags += "LOW_COMPLEXITY,";
            if (!gt.lrt_significant) flags += "LRT_NS,";
            if (!flags.empty()) {
                flags.pop_back();  // remove trailing comma
                ofs_ << ";FLAGS=" << flags;
            }
        }

        if (bundle.from_rescue) {
            ofs_ << ";RESCUED=1";
        }

        ofs_ << "\t";

        // FORMAT
        ofs_ << "GT:GQ:AF\t";

        // SAMPLE
        if (bundle.has_genotype) {
            ofs_ << gt.genotype << ":"
                 << gt.gq << ":"
                 << std::setprecision(4) << gt.af;
        } else {
            ofs_ << "./.:.:.";
        }

        ofs_ << "\n";
    }

private:
    std::ofstream ofs_;
    std::string sample_name_;

    static std::string tier_to_string(Tier t) {
        switch (t) {
            case Tier::TIER1: return "1";
            case Tier::TIER2: return "2";
            default: return "3";
        }
    }

    void write_header() {
        ofs_ << "##fileformat=VCFv4.2\n";
        ofs_ << "##source=PLACER_v0.2\n";
        ofs_ << "##INFO=<ID=TIER,Number=1,Type=String,"
                "Description=\"Tier classification (1/2/3)\">\n";
        ofs_ << "##INFO=<ID=PLACEABILITY,Number=1,Type=Float,"
                "Description=\"Placeability delta score\">\n";
        ofs_ << "##INFO=<ID=SUPPORT_READS,Number=1,Type=Integer,"
                "Description=\"Number of supporting reads\">\n";
        ofs_ << "##INFO=<ID=SVLEN,Number=1,Type=Integer,"
                "Description=\"Insertion length estimate\">\n";
        ofs_ << "##INFO=<ID=TE_FAMILY,Number=1,Type=Integer,"
                "Description=\"TE family ID\">\n";
        ofs_ << "##INFO=<ID=TE_ORIENT,Number=1,Type=String,"
                "Description=\"TE orientation (+/-)\">\n";
        ofs_ << "##INFO=<ID=BP_LEFT,Number=1,Type=Integer,"
                "Description=\"Left breakpoint\">\n";
        ofs_ << "##INFO=<ID=BP_RIGHT,Number=1,Type=Integer,"
                "Description=\"Right breakpoint\">\n";
        ofs_ << "##INFO=<ID=TSD_LEN,Number=1,Type=Integer,"
                "Description=\"Target site duplication length\">\n";
        ofs_ << "##INFO=<ID=AF,Number=1,Type=Float,"
                "Description=\"Allele frequency\">\n";
        ofs_ << "##INFO=<ID=AF_CI,Number=2,Type=Float,"
                "Description=\"AF confidence interval\">\n";
        ofs_ << "##INFO=<ID=GQ,Number=1,Type=Integer,"
                "Description=\"Genotype quality (Phred)\">\n";
        ofs_ << "##INFO=<ID=EALT,Number=1,Type=Float,"
                "Description=\"Expected ALT fraction\">\n";
        ofs_ << "##INFO=<ID=EREF,Number=1,Type=Float,"
                "Description=\"Expected REF fraction\">\n";
        ofs_ << "##INFO=<ID=ENULL,Number=1,Type=Float,"
                "Description=\"Expected NULL fraction\">\n";
        ofs_ << "##INFO=<ID=FLAGS,Number=.,Type=String,"
                "Description=\"Quality flags\">\n";
        ofs_ << "##INFO=<ID=RESCUED,Number=0,Type=Flag,"
                "Description=\"Locus recovered by TE reverse index\">\n";
        ofs_ << "##FORMAT=<ID=GT,Number=1,Type=String,"
                "Description=\"Genotype\">\n";
        ofs_ << "##FORMAT=<ID=GQ,Number=1,Type=Integer,"
                "Description=\"Genotype quality\">\n";
        ofs_ << "##FORMAT=<ID=AF,Number=1,Type=Float,"
                "Description=\"Allele frequency\">\n";
        ofs_ << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
             << sample_name_ << "\n";
    }
};

// ============================================================================
// Main Pipeline
// ============================================================================

class PLACERPipeline {
public:
    explicit PLACERPipeline(const PipelineConfig& config)
        : config_(config),
          component_builder_(create_builder_config()),
          local_realigner_(create_realign_config()),
          assembly_engine_(create_assembly_config()),
          placeability_scorer_(create_placeability_config()),
          genotyper_(create_genotype_config()),
          genome_accessor_(config.genome_fasta_path) {

        // 初始化 Gate1（如果启用）
        if (config.gate1_enabled && !config.te_fasta_path.empty()) {
            Gate1Config g1cfg;
            g1cfg.probe_len = 200;
            g1cfg.min_hit_count = 3;
            g1cfg.min_density = 0.1;

            HashTEIndexConfig htcfg;
            htcfg.kmer_size = 15;
            htcfg.min_hit_count = 3;
            htcfg.min_hit_density = 0.05;

            auto te_index = HashTEIndex::build_from_fasta(
                config.te_fasta_path, htcfg);
            if (!te_index) {
                throw std::runtime_error(
                    "Failed to build TE index from: " + config.te_fasta_path);
            }

            gate1_ = std::make_unique<Gate1>(
                std::move(te_index), g1cfg);
        }

        // 初始化 TE Reverse Index（如果启用）
        if (config.te_reverse_enabled && !config.genome_fasta_path.empty()) {
            TEReverseIndexConfig teri_cfg;
            if (gate1_) {
                teri_cfg.gate1_config = gate1_->get_config();  // [修正] 使用 get_config()
            }
            te_reverse_index_ = std::make_unique<TEReverseIndex>(teri_cfg);
            te_reverse_index_->initialize(config.genome_fasta_path);
        }

        // 从 genome_accessor 获取 chrom 名称映射
        build_chrom_name_map();
    }

    PipelineResult run(const std::string& bam_path,
                       const std::vector<ReadSketch>& raw_reads) {

        PipelineResult result;

        std::cout << "=== PLACER Pipeline ===" << std::endl;
        auto start_time = std::chrono::high_resolution_clock::now();

        result.total_reads = raw_reads.size();
        std::cout << "Input reads: " << result.total_reads << std::endl;

        // ================================================================
        // Phase 2: Gate1 filtering
        // ================================================================
        std::cout << "\n[Phase 2] Gate1 TE-proxy filtering..." << std::endl;
        std::vector<ReadSketch> filtered_reads;
        std::vector<std::vector<ProbeFragment>> filtered_probes;  // [修正 P0.1] 与 filtered_reads 对齐

        if (gate1_) {
            gate1_filter(raw_reads, filtered_reads, filtered_probes);
            result.gate1_passed = filtered_reads.size();
            std::cout << "  Passed Gate1: " << result.gate1_passed
                      << " / " << raw_reads.size() << std::endl;
        } else {
            // Gate1 未启用：所有 reads 直接通过
            filtered_reads = raw_reads;
            result.gate1_passed = filtered_reads.size();
            std::cout << "  Gate1 disabled, using all reads" << std::endl;
        }

        // ================================================================
        // Phase 3: Component Building
        // ================================================================
        std::cout << "\n[Phase 3] Building components..." << std::endl;
        auto components = component_builder_.build(filtered_reads);
        result.built_components = components.size();
        std::cout << "  Built " << result.built_components
                  << " components" << std::endl;

        // ================================================================
        // Phase 4: Local Realignment → 真实 LocusEvidence
        // [修正] 使用 realign_and_collect 真正的序列比对
        // ================================================================
        std::cout << "\n[Phase 4] Local realignment..." << std::endl;

        std::vector<std::vector<LocusEvidence>> all_evidence;
        std::vector<PlaceabilityReport> all_reports;

        all_evidence.reserve(components.size());
        all_reports.reserve(components.size());

        for (auto& comp : components) {
            auto result = local_realigner_.realign_and_collect(
                comp, filtered_reads, genome_accessor_);
            all_evidence.push_back(std::move(result.evidence));
            all_reports.push_back(std::move(result.report));
        }

        std::cout << "  Collected evidence for "
                  << all_evidence.size() << " components" << std::endl;

        // ================================================================
        // Phase 5: Assembly → StructuralRepresentative
        // ================================================================
        std::cout << "\n[Phase 5] Assembly..." << std::endl;
        std::vector<Contig> all_contigs;
        std::vector<StructuralRepresentative> representatives;

        for (const auto& comp : components) {
            // 只处理 5-50 reads 的组件
            if (comp.read_count < 5 || comp.read_count > 50) continue;

            auto contigs = assembly_engine_.assemble_component(
                comp, filtered_reads, genome_accessor_);

            for (auto& c : contigs) {
                c.source_component_id = comp.id;
                all_contigs.push_back(std::move(c));
            }
        }

        // 结构级合并
        representatives = assembly_engine_.collapse_structurally(all_contigs);
        result.assembled_contigs = representatives.size();
        std::cout << "  Assembled " << representatives.size()
                  << " structural representatives" << std::endl;

        // ================================================================
        // Phase 8: TE Reverse Index recall（在 placeability/genotyping 之前）
        // 把 rescued loci 作为新 component 加入流程
        // ================================================================
        if (te_reverse_index_) {
            std::cout << "\n[Phase 8] TE Reverse Index recall..."
                      << std::endl;
            auto rescued = te_reverse_recall(
                filtered_reads, components, filtered_probes);  // [修正 P0.1] 使用对齐后的 probes
            result.rescued_loci = rescued.size();
            std::cout << "  Rescued " << rescued.size() << " loci"
                      << std::endl;

            // 把 rescued loci 转换为额外的 component + evidence + rep
            integrate_rescued_loci(
                rescued, filtered_reads,
                components, all_evidence, all_reports, representatives);

            std::cout << "  After integration: "
                      << components.size() << " components, "
                      << representatives.size() << " representatives"
                      << std::endl;
        }

        // ================================================================
        // 绑定 representative ↔ component ↔ evidence
        // ================================================================
        auto bundles = build_call_bundles(
            components, all_evidence, representatives);
        std::cout << "\n  Call bundles: " << bundles.size() << std::endl;

        // ================================================================
        // Phase 6: Placeability scoring
        // [修正] 使用 realign_and_collect 产生的 all_reports
        // ================================================================
        std::cout << "\n[Phase 6] Placeability scoring..." << std::endl;
        int tier1 = 0, tier2 = 0, tier3 = 0;

        // 建立 component_id → report 映射
        std::unordered_map<int32_t, PlaceabilityReport> comp_id_to_report;
        for (size_t i = 0; i < components.size() && i < all_reports.size(); ++i) {
            comp_id_to_report[components[i].id] = all_reports[i];
        }

        for (auto& bundle : bundles) {
            if (bundle.evidence.empty()) continue;

            // 使用已有的 report
            auto it = comp_id_to_report.find(bundle.component_id);
            if (it != comp_id_to_report.end()) {
                // 从 PlaceabilityReport 复制到 ExtendedPlaceabilityReport
                const auto& src = it->second;
                bundle.placeability.best_locus = src.best_locus;
                bundle.placeability.second_best_locus = -1;
                bundle.placeability.delta_score = src.delta_score;
                bundle.placeability.tier = static_cast<Tier>(src.tier);
                bundle.placeability.candidate_count = 1;
                bundle.placeability.support_reads = src.best_support_count;
                bundle.placeability.forward_count = src.forward_count;
                bundle.placeability.reverse_count = src.reverse_count;
                bundle.placeability.strand_balanced = src.strand_balanced;
                bundle.placeability.strand_ratio = src.strand_ratio;
                bundle.placeability.is_ambiguous = src.is_ambiguous;
                bundle.has_placeability = true;

                if (bundle.placeability.tier == Tier::TIER1) tier1++;
                else if (bundle.placeability.tier == Tier::TIER2) tier2++;
                else tier3++;
            } else {
                // Fallback: 重新计算
                bundle.placeability =
                    placeability_scorer_.calculate_full_placeability(
                        bundle.evidence);
                bundle.has_placeability = true;
                tier3++;
            }
        }

        result.tier1_loci = tier1;
        result.tier2_loci = tier2;
        result.tier3_loci = tier3;
        std::cout << "  Tier1: " << tier1
                  << ", Tier2: " << tier2
                  << ", Tier3: " << tier3 << std::endl;

        // ================================================================
        // Phase 7: Genotyping
        // [修正] 传入真实 reads 和完整 representative
        // ================================================================
        std::cout << "\n[Phase 7] Genotyping..." << std::endl;
        int genotyped_count = 0;

        for (auto& bundle : bundles) {
            if (bundle.evidence.empty()) continue;
            if (!bundle.has_representative) continue;

            bundle.genotype = genotyper_.genotype(
                bundle.representative,
                bundle.evidence,
                filtered_reads,       // [修正] 传入真实 reads
                genome_accessor_);
            bundle.has_genotype = true;
            genotyped_count++;
        }

        result.genotyped = genotyped_count;
        std::cout << "  Genotyped " << genotyped_count
                  << " loci" << std::endl;

        // ================================================================
        // Output VCF
        // ================================================================
        std::cout << "\n[Output] Writing VCF..." << std::endl;
        result.vcf_path = "scientific.vcf";
        write_vcf(bundles, result.vcf_path);

        // 存储结果
        result.representatives = std::move(representatives);
        for (auto& b : bundles) {
            if (b.has_genotype) {
                result.genotypes.push_back(b.genotype);
            }
            if (b.has_placeability) {
                result.placeability_reports.push_back(b.placeability);
            }
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            end_time - start_time);
        std::cout << "\nTotal time: " << duration.count() << "ms"
                  << std::endl;

        return result;
    }

private:
    PipelineConfig config_;

    // Modules
    ComponentBuilder component_builder_;
    LocalRealigner local_realigner_;
    AssemblyEngine assembly_engine_;
    PlaceabilityScorer placeability_scorer_;
    Genotyper genotyper_;
    GenomeAccessor genome_accessor_;
    std::unique_ptr<Gate1> gate1_;
    std::unique_ptr<TEReverseIndex> te_reverse_index_;

    // Chrom name mapping: tid → name
    std::unordered_map<int32_t, std::string> chrom_names_;

    // ====================================================================
    // Config factory methods
    // ====================================================================

    ComponentBuilderConfig create_builder_config() {
        ComponentBuilderConfig cfg;
        cfg.min_density = config_.component_min_density;
        cfg.max_cluster_span = config_.component_breaker_threshold;
        cfg.max_recursive_depth = config_.component_max_recursive_depth;
        return cfg;
    }

    RealignConfig create_realign_config() {
        RealignConfig cfg;
        cfg.flank_length = config_.flank_length;
        cfg.search_window = config_.search_window;
        cfg.min_normalized_score = config_.min_realign_score;
        cfg.max_locus_per_component = config_.max_locus_per_component;
        return cfg;
    }

    AssemblyConfig create_assembly_config() {
        AssemblyConfig cfg;
        cfg.min_reads_for_poa = config_.assembly_min_reads;
        cfg.max_reads_for_poa = config_.assembly_max_reads;
        cfg.max_output_paths = config_.assembly_max_paths;
        return cfg;
    }

    PlaceabilityConfig create_placeability_config() {
        PlaceabilityConfig cfg;
        cfg.delta_score_threshold = config_.delta_tier1;
        cfg.delta_score_tier2 = config_.delta_tier2;
        cfg.side_consistency_gap = config_.side_consistency_gap;
        cfg.min_locus_support = config_.min_locus_support;
        return cfg;
    }

    GenotypeConfig create_genotype_config() {
        GenotypeConfig cfg;
        cfg.max_em_iterations = config_.max_em_iterations;
        cfg.em_convergence_threshold = config_.em_convergence;
        cfg.spatial_lambda = config_.spatial_lambda;
        return cfg;
    }

    // ====================================================================
    // Chrom name mapping
    // ====================================================================

    void build_chrom_name_map() {
        // 从 genome_accessor 获取 contig 名称
        auto names = genome_accessor_.get_contig_names();
        for (size_t i = 0; i < names.size(); ++i) {
            chrom_names_[static_cast<int32_t>(i)] = names[i];
        }
    }

    std::string get_chrom_name(int32_t tid) const {
        auto it = chrom_names_.find(tid);
        if (it != chrom_names_.end()) {
            return it->second;
        }
        // fallback
        return "chr" + std::to_string(tid);
    }

    // ====================================================================
    // Phase 2: Gate1 filtering
    // [修正 P0.1] 维护与 filtered_reads 同步的 probes 向量
    // 不再使用 raw_reads 索引，而是：
    // - filtered_reads[idx] 的 probes 在 filtered_probes[idx]
    // ====================================================================

    void gate1_filter(
        const std::vector<ReadSketch>& raw_reads,
        std::vector<ReadSketch>& filtered_reads,
        std::vector<std::vector<ProbeFragment>>& filtered_probes) {

        filtered_reads.clear();
        filtered_probes.clear();

        for (size_t i = 0; i < raw_reads.size(); ++i) {
            const auto& read = raw_reads[i];

            // SA reads 直接通过（它们已有定位信息）
            if (read.has_sa) {
                filtered_reads.push_back(read);
                // SA reads 没有从 Gate1 来的 probes，存入空 vector
                filtered_probes.emplace_back();
                continue;
            }

            // 通过 Gate1 评估
            auto result = gate1_->evaluate(read);

            if (result.passed) {
                filtered_reads.push_back(read);
                // [修正] probes 存入 filtered_probes 的末尾（与 filtered_reads 对齐）
                filtered_probes.push_back(std::move(result.probes));
            }
        }
    }

    // ====================================================================
    // Phase 5: Assembly
    //
    // [修正] assembly 输出的 representative 与 component 绑定：
    //   - representative.component_ids 记录来源 component
    //   - representative 包含 fingerprint/breakpoints/sequence/family
    // ====================================================================

    std::vector<StructuralRepresentative> assemble(
        const std::vector<Component>& components,
        const std::vector<ReadSketch>& reads) {

        std::vector<Contig> all_contigs;

        for (const auto& comp : components) {
            auto contigs = assembly_engine_.assemble_component(
                comp, reads, genome_accessor_);

            for (auto& c : contigs) {
                // 标记来源 component
                c.source_component_id = comp.id;
                all_contigs.push_back(std::move(c));
            }
        }

        // 结构级合并
        auto representatives =
            assembly_engine_.collapse_structurally(all_contigs);

        return representatives;
    }

    // ====================================================================
    // Phase 8: TE Reverse Index recall
    //
    // [修正] rescued loci 集成进 component/evidence/representative 流程
    // ====================================================================

    std::vector<RescuedLocus> te_reverse_recall(
        const std::vector<ReadSketch>& reads,
        const std::vector<Component>& existing_components,
        const std::vector<std::vector<ProbeFragment>>& filtered_probes) {

        if (!te_reverse_index_) return {};

        // 方式 1：如果 filtered_probes 不为空（Gate1 过滤后），
        // 检查是否有非空的 probes（跳过 SA reads 的空 probes）
        bool has_valid_probes = false;
        for (const auto& probes : filtered_probes) {
            if (!probes.empty()) {
                has_valid_probes = true;
                break;
            }
        }

        if (has_valid_probes) {
            // [修正] 传入 reads 和 filtered_probes
            return te_reverse_index_->rescue_with_evidence(
                reads, filtered_probes);
        }

        // 方式 2：否则用 rescue_reads（内部会调用 Gate1 提取 probes）
        std::vector<bool> has_sa(reads.size(), false);
        for (size_t i = 0; i < reads.size(); ++i) {
            has_sa[i] = reads[i].has_sa;
        }

        return te_reverse_index_->rescue_reads(
            reads, existing_components, has_sa);
    }

    // ====================================================================
    // 把 rescued loci 转换为 component + evidence + representative
    // 并追加到现有数据中
    // ====================================================================

    void integrate_rescued_loci(
        const std::vector<RescuedLocus>& rescued,
        const std::vector<ReadSketch>& reads,
        std::vector<Component>& components,
        std::vector<std::vector<LocusEvidence>>& all_evidence,
        std::vector<PlaceabilityReport>& all_reports,
        std::vector<StructuralRepresentative>& representatives) {

        int32_t next_comp_id = components.empty() ?
            0 : components.back().id + 1;

        for (const auto& locus : rescued) {
            if (!locus.passed_filter) continue;

            // [修正 P1.3] 使用更合理的 overlap 判定策略
            // 基于 breakpoint tolerance 和 locus_cluster_radius，而非粗粒的 search_window/2
            static constexpr int32_t kOverlapTolerance = 50;  // 断点容忍度 (bp)
            bool overlaps = false;
            for (const auto& comp : components) {
                // 检查染色体是否匹配
                if (comp.chrom_tid != locus.chrom_tid) continue;

                // 检查 component 范围是否重叠
                int32_t comp_center = (comp.start + comp.end) / 2;
                int32_t locus_center = locus.position;
                int32_t distance = std::abs(comp_center - locus_center);
                int32_t comp_span = comp.end - comp.start;

                // 如果距离小于两者的 span 之和，判定为重叠
                if (distance < (comp_span / 2 + kOverlapTolerance)) {
                    overlaps = true;
                    break;
                }
            }

            if (overlaps) continue;  // 已有覆盖，跳过

            // 创建新 component
            Component new_comp;
            new_comp.id = next_comp_id++;
            new_comp.read_indices = locus.probe_read_indices;  // [修正 P0.3] 语法错误
            new_comp.read_count = locus.support_reads;

            LocusCandidate lc;
            lc.chrom_tid = locus.chrom_tid;
            lc.pos = locus.position;
            lc.score = locus.placeability_score;
            lc.support_reads = locus.support_reads;
            new_comp.locus_set.push_back(lc);

            components.push_back(new_comp);

            // 对新 component 做 local realignment 获取真实 evidence
            // [修正] 使用 realign_and_collect
            auto result = local_realigner_.realign_and_collect(
                new_comp, reads, genome_accessor_);
            all_evidence.push_back(std::move(result.evidence));
            all_reports.push_back(std::move(result.report));

            // 对新 component 做 assembly
            auto contigs = assembly_engine_.assemble_component(
                new_comp, reads, genome_accessor_);

            if (!contigs.empty()) {
                for (auto& c : contigs) {
                    c.source_component_id = new_comp.id;
                }
                auto new_reps =
                    assembly_engine_.collapse_structurally(contigs);
                for (auto& rep : new_reps) {
                    rep.from_rescue = true;
                    representatives.push_back(std::move(rep));
                }
            } else {
                // 没有 assembly 结果：创建一个最小 representative
                StructuralRepresentative min_rep;
                min_rep.rep_id = new_comp.id;
                min_rep.component_ids.push_back(new_comp.id);
                min_rep.total_reads = new_comp.read_count;
                min_rep.fingerprint.tid = locus.chrom_tid;
                min_rep.fingerprint.breakpoint_l = locus.position;
                min_rep.fingerprint.breakpoint_r = locus.position;
                min_rep.from_rescue = true;
                representatives.push_back(std::move(min_rep));
            }
        }
    }

    // ====================================================================
    // 绑定 representative ↔ component ↔ evidence → LocusCallBundle
    //
    // 策略：
    //   1. 每个 representative 通过 component_ids 找到对应 component
    //   2. 每个 component 通过 index 找到对应 evidence
    //   3. 如果某个 component 没有 representative，仍然创建 bundle
    //      （只有 evidence，没有 assembly 结果）
    // ====================================================================

    std::vector<LocusCallBundle> build_call_bundles(
        const std::vector<Component>& components,
        const std::vector<std::vector<LocusEvidence>>& all_evidence,
        const std::vector<StructuralRepresentative>& representatives) {

        // component_id → index in components vector
        std::unordered_map<int32_t, size_t> comp_id_to_idx;
        for (size_t i = 0; i < components.size(); ++i) {
            comp_id_to_idx[components[i].id] = i;
        }

        // 标记哪些 component 已被 representative 覆盖
        std::unordered_set<int32_t> covered_comp_ids;

        std::vector<LocusCallBundle> bundles;

        // 1. 从 representative 出发创建 bundle
        for (size_t ri = 0; ri < representatives.size(); ++ri) {
            const auto& rep = representatives[ri];

            LocusCallBundle bundle;
            bundle.representative = rep;
            bundle.has_representative = true;
            bundle.representative_idx = ri;
            bundle.from_rescue = rep.from_rescue;

            // 找到对应的 component 和 evidence
            // representative 可能关联多个 component（collapse 后）
            // 取第一个有 evidence 的 component
            for (int32_t comp_id : rep.component_ids) {
                covered_comp_ids.insert(comp_id);

                auto it = comp_id_to_idx.find(comp_id);
                if (it == comp_id_to_idx.end()) continue;

                size_t ci = it->second;
                const auto& comp = components[ci];

                // 设置坐标
                if (bundle.chrom_tid < 0) {
                    bundle.component_id = comp.id;
                    if (!comp.locus_set.empty()) {
                        bundle.chrom_tid = comp.locus_set[0].chrom_tid;
                        bundle.position = comp.locus_set[0].pos;
                    } else if (rep.fingerprint.tid >= 0) {
                        bundle.chrom_tid = rep.fingerprint.tid;
                        bundle.position = rep.fingerprint.breakpoint_l;
                    }
                }

                // [修正 P1.4] 合并 evidence 时去重（按 read_idx + locus_pos + side）
                if (ci < all_evidence.size()) {
                    std::unordered_set<std::string> seen_evidence;
                    for (const auto& ev : all_evidence[ci]) {
                        // 构建唯一键：read_idx + locus_pos + evidence_bits
                        std::string key = std::to_string(ev.read_idx) + "_" +
                                          std::to_string(ev.locus_pos) + "_" +
                                          std::to_string(ev.evidence_bits);
                        if (!seen_evidence.count(key)) {
                            seen_evidence.insert(key);
                            bundle.evidence.push_back(ev);
                        }
                    }
                }
            }

            // 使用 fingerprint 的坐标作为 fallback
            if (bundle.chrom_tid < 0 && rep.fingerprint.tid >= 0) {
                bundle.chrom_tid = rep.fingerprint.tid;
                bundle.position = rep.fingerprint.breakpoint_l;
            }

            // 设置 chrom name
            bundle.chrom_name = get_chrom_name(bundle.chrom_tid);

            // 如果 representative 有更精确的 breakpoint，使用它
            if (rep.fingerprint.breakpoint_l > 0) {
                bundle.position = rep.fingerprint.breakpoint_l;
            }

            bundles.push_back(std::move(bundle));
        }

        // 2. 没有被 representative 覆盖的 component 也创建 bundle
        //    （assembly 失败但仍有 evidence 的情况）
        for (size_t ci = 0; ci < components.size(); ++ci) {
            const auto& comp = components[ci];

            if (covered_comp_ids.count(comp.id)) continue;
            if (ci >= all_evidence.size() || all_evidence[ci].empty()) continue;

            LocusCallBundle bundle;
            bundle.component_id = comp.id;
            bundle.has_representative = false;

            if (!comp.locus_set.empty()) {
                bundle.chrom_tid = comp.locus_set[0].chrom_tid;
                bundle.position = comp.locus_set[0].pos;
            }

            bundle.chrom_name = get_chrom_name(bundle.chrom_tid);
            bundle.evidence = all_evidence[ci];

            bundles.push_back(std::move(bundle));
        }

        // 按 (chrom_tid, position) 排序
        std::sort(bundles.begin(), bundles.end(),
            [](const LocusCallBundle& a, const LocusCallBundle& b) {
                if (a.chrom_tid != b.chrom_tid)
                    return a.chrom_tid < b.chrom_tid;
                return a.position < b.position;
            });

        return bundles;
    }

    // ====================================================================
    // Output VCF
    // ====================================================================

    void write_vcf(const std::vector<LocusCallBundle>& bundles,
                   const std::string& vcf_path) {
        try {
            VCFWriter writer(vcf_path, "SAMPLE");

            int written = 0;
            for (const auto& bundle : bundles) {
                if (!bundle.has_placeability) continue;
                if (bundle.evidence.empty()) continue;

                writer.write_record(bundle);
                written++;
            }

            std::cout << "  VCF written: " << written
                      << " records to " << vcf_path << std::endl;

        } catch (const std::exception& e) {
            std::cerr << "Warning: Could not write VCF: "
                      << e.what() << std::endl;
        }
    }
};

}  // namespace placer

// ============================================================================
// Main
// ============================================================================

using namespace placer;

void print_usage(const char* prog) {
    std::cout << "Usage: " << prog
              << " <bam_file> <ref_fasta> [te_fasta]" << std::endl;
    std::cout << std::endl;
    std::cout << "PLACER v0.2 - TE Insertion Detection Pipeline"
              << std::endl;
    std::cout << std::endl;
    std::cout << "Arguments:" << std::endl;
    std::cout << "  bam_file    Input BAM file (required)" << std::endl;
    std::cout << "  ref_fasta   Reference genome FASTA "
                 "(required, indexed with .fai)" << std::endl;
    std::cout << "  te_fasta    TE library FASTA "
                 "(optional, enables Gate1 + TE Reverse Index)"
              << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        print_usage(argv[0]);
        return 1;
    }

    std::string bam_path = argv[1];
    std::string ref_fasta_path = argv[2];
    std::string te_fasta_path = (argc >= 4) ? argv[3] : "";

    std::cout << "=== PLACER v0.2 - TE Insertion Detection ==="
              << std::endl;
    std::cout << "Input BAM: " << bam_path << std::endl;
    std::cout << "Reference: " << ref_fasta_path << std::endl;
    if (!te_fasta_path.empty()) {
        std::cout << "TE Library: " << te_fasta_path << std::endl;
    }

    // Configure pipeline
    PipelineConfig config;
    config.genome_fasta_path = ref_fasta_path;

    if (!te_fasta_path.empty()) {
        config.gate1_enabled = true;
        config.te_fasta_path = te_fasta_path;
        config.te_reverse_enabled = true;
    }

    // Phase 1: Streaming BAM with window buffering
    // Uses reference genome for chromosome count, not BAM header
    std::cout << "\n[Phase 1] Streaming BAM with window buffering..." << std::endl;

    // First, get chromosome count from reference genome
    GenomeAccessor genome_accessor(ref_fasta_path);
    int32_t ref_chromosomes = static_cast<int32_t>(genome_accessor.num_chroms());
    std::cout << "  Reference genome: " << ref_chromosomes << " chromosomes" << std::endl;

    // Open BAM for streaming
    BamReader reader(bam_path);
    if (!reader.is_valid()) {
        std::cerr << "Error: Cannot open BAM file: " << bam_path
                  << std::endl;
        return 1;
    }
    int32_t bam_chromosomes = reader.get_num_chromosomes();
    std::cout << "  BAM header: " << bam_chromosomes << " chromosomes" << std::endl;

    // Configure streaming processor
    AccumulatingStreamProcessor::Config stream_config;
    stream_config.window_size = config.window_size;
    stream_config.window_step = config.window_step;
    stream_config.min_clip_bp = config.min_clip_bp;
    stream_config.min_sa_reads = config.min_sa_reads;
    stream_config.min_density = static_cast<int>(config.component_min_density * 100);
    stream_config.max_cluster_span = config.component_breaker_threshold;
    stream_config.max_recursive_depth = config.component_max_recursive_depth;
    stream_config.verbose = true;
    stream_config.progress_interval = 500000;

    AccumulatingStreamProcessor stream_processor(stream_config);
    auto stream_result = stream_processor.process(reader);

    std::cout << "\n=== Pipeline Summary ===" << std::endl;
    std::cout << "Total reads:       " << stream_result.total_reads << std::endl;
    std::cout << "Triggered windows: " << stream_result.triggered_windows << std::endl;
    std::cout << "Components built:  " << stream_result.built_components << std::endl;
    std::cout << "Reference chroms:  " << ref_chromosomes << std::endl;

    // Continue with full pipeline if components were found
    if (stream_result.built_components > 0) {
        std::cout << "\n[Phase 2+] Processing " << stream_result.built_components
                  << " components..." << std::endl;
        // Continue with remaining pipeline phases...
    }

    return 0;
}

