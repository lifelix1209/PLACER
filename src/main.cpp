#include "denovo.h"
#include "null_control.h"
#include "pipeline.h"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace placer
{
    namespace
    {

        bool env_flag_enabled(const char *key)
        {
            const char *value = std::getenv(key);
            if (!value)
            {
                return false;
            }
            const std::string v(value);
            return v == "1" || v == "true" || v == "TRUE" || v == "on" || v == "ON";
        }

        bool env_try_bool(const char *key, bool &out)
        {
            const char *value = std::getenv(key);
            if (!value)
            {
                return false;
            }
            const std::string v(value);
            if (v == "1" || v == "true" || v == "TRUE" || v == "on" || v == "ON")
            {
                out = true;
                return true;
            }
            if (v == "0" || v == "false" || v == "FALSE" || v == "off" || v == "OFF")
            {
                out = false;
                return true;
            }
            return false;
        }

        bool env_try_double(const char *key, double &out)
        {
            const char *value = std::getenv(key);
            if (!value || !*value)
            {
                return false;
            }
            char *end = nullptr;
            const double parsed = std::strtod(value, &end);
            if (end == value || (end && *end != '\0'))
            {
                return false;
            }
            out = parsed;
            return true;
        }

        bool env_try_int32(const char *key, int32_t &out)
        {
            const char *value = std::getenv(key);
            if (!value || !*value)
            {
                return false;
            }
            char *end = nullptr;
            const long parsed = std::strtol(value, &end, 10);
            if (end == value || (end && *end != '\0'))
            {
                return false;
            }
            out = static_cast<int32_t>(parsed);
            return true;
        }

        bool env_try_string(const char *key, std::string &out)
        {
            const char *value = std::getenv(key);
            if (!value)
            {
                return false;
            }
            out = value;
            return true;
        }

        std::string safe_current_path()
        {
            try
            {
                return std::filesystem::current_path().string();
            }
            catch (const std::exception &)
            {
                return ".";
            }
        }

        std::string safe_absolute_path(const std::string &path)
        {
            try
            {
                return std::filesystem::absolute(path).string();
            }
            catch (const std::exception &)
            {
                return path;
            }
        }

        std::string serialize_support_qnames(const std::vector<std::string> &qnames)
        {
            std::ostringstream out;
            bool first = true;
            for (const auto &qname : qnames)
            {
                if (!first)
                {
                    out << ",";
                }
                first = false;
                for (char ch : qname)
                {
                    if (ch == '\t' || ch == '\n' || ch == '\r' || ch == ',')
                    {
                        out << "_";
                    }
                    else
                    {
                        out << ch;
                    }
                }
            }
            return out.str();
        }

        bool include_candidate_insert_seq_debug_output()
        {
            return env_flag_enabled("PLACER_DEBUG_OUTPUT_CANDIDATE_INSERT_SEQ") ||
                   env_flag_enabled("PLACER_DEBUG_OUTPUT_INSERT_SEQ");
        }

        void print_usage()
        {
            std::cerr << "Usage:\n"
                      << "  placer [--region <chrom:start-end>] [--final-fdr-q <q>] "
                         "[--final-report-mode <legacy|te-calibrated>] "
                         "[--min-final-raw-cigar-insert-len-bp <bp>] <input.bam> <ref.fa> <te.fa>"
                      << std::endl;
        }

        const char *final_report_mode_name(FinalReportMode mode)
        {
            switch (mode)
            {
            case FinalReportMode::Legacy:
                return "legacy";
            case FinalReportMode::TeCalibrated:
                return "te-calibrated";
            }
            return "unknown";
        }

        bool parse_final_report_mode(const std::string &text, FinalReportMode &out)
        {
            if (text == "legacy")
            {
                out = FinalReportMode::Legacy;
                return true;
            }
            if (text == "te-calibrated")
            {
                out = FinalReportMode::TeCalibrated;
                return true;
            }
            return false;
        }

        bool parse_positive_i64(const std::string &text, int64_t &out)
        {
            if (text.empty())
            {
                return false;
            }
            for (char ch : text)
            {
                if (!std::isdigit(static_cast<unsigned char>(ch)))
                {
                    return false;
                }
            }
            try
            {
                out = std::stoll(text);
            }
            catch (const std::exception &)
            {
                return false;
            }
            return out > 0;
        }

        BamRegionScope parse_region_scope_or_throw(const std::string &region)
        {
            if (region.empty())
            {
                throw std::runtime_error("empty region");
            }

            BamRegionScope scope;
            scope.enabled = true;

            const size_t colon = region.rfind(':');
            if (colon == std::string::npos)
            {
                scope.chrom = region;
                scope.start = 0;
                scope.end = -1;
                return scope;
            }

            scope.chrom = region.substr(0, colon);
            const std::string range = region.substr(colon + 1);
            const size_t dash = range.find('-');
            if (scope.chrom.empty() || dash == std::string::npos)
            {
                throw std::runtime_error("invalid region: " + region);
            }

            int64_t start_1based = 0;
            int64_t end_1based = 0;
            if (!parse_positive_i64(range.substr(0, dash), start_1based) ||
                !parse_positive_i64(range.substr(dash + 1), end_1based) ||
                end_1based < start_1based ||
                end_1based > std::numeric_limits<int32_t>::max())
            {
                throw std::runtime_error("invalid region: " + region);
            }

            scope.start = static_cast<int32_t>(start_1based - 1);
            scope.end = static_cast<int32_t>(end_1based);
            return scope;
        }

        void apply_environment_config(PipelineConfig &config)
        {
            config.enable_parallel = env_flag_enabled("PLACER_PARALLEL");

            double v = 0.0;
            bool b = false;
            int32_t i = 0;
            if (env_try_int32("PLACER_PROGRESS_INTERVAL", i))
            {
                config.progress_interval = std::max<int64_t>(0, static_cast<int64_t>(i));
            }
            if (env_try_bool("PLACER_LOG_STAGE_BINS", b))
            {
                config.log_stage_bins = b;
            }
            if (env_try_bool("PLACER_LOG_STAGE_COMPONENTS", b))
            {
                config.log_stage_components = b;
            }
            std::string s;
            if (env_try_string("PLACER_INS_FRAGMENTS_FASTA_PATH", s))
            {
                config.ins_fragments_fasta_path = s;
            }
            if (env_try_string("PLACER_INS_FRAGMENT_HITS_TSV_PATH", s))
            {
                config.ins_fragment_hits_tsv_path = s;
            }
            if (env_try_int32("PLACER_TE_KMER_SIZE", i))
            {
                config.te_kmer_size = std::max(7, i);
            }
            if (env_try_string("PLACER_TE_KMER_SIZES", s))
            {
                config.te_kmer_sizes_csv = s;
            }
            if (env_try_int32("PLACER_TE_FAMILY_TOPN", i))
            {
                config.te_family_topn = std::max(1, i);
            }
            if (env_try_int32("PLACER_TE_FAMILY_REPRESENTATIVES", i))
            {
                config.te_family_representatives = std::max(1, i);
            }
            if (env_try_int32("PLACER_TE_TEMPLATE_REFINE_TOPN", i))
            {
                config.te_template_refine_topn = std::max(1, i);
            }
            if (env_try_int32("PLACER_TE_EXACT_ALIGN_TOPN", i))
            {
                config.te_exact_align_topn = std::max(1, i);
            }
            if (env_try_double("PLACER_TE_FAMILY_MARGIN_MIN", v))
            {
                config.te_family_margin_min = std::clamp(v, 0.0, 1.0);
            }
            if (env_try_double("PLACER_TE_SUBFAMILY_MARGIN_MIN", v))
            {
                config.te_subfamily_margin_min = std::clamp(v, 0.0, 1.0);
            }
            if (env_try_int32("PLACER_BAM_THREADS", i))
            {
                config.bam_threads = std::max(1, i);
            }
            if (env_try_int32("PLACER_PARALLEL_WORKERS", i))
            {
                config.parallel_workers = std::max(1, i);
            }
            if (env_try_int32("PLACER_PARALLEL_QUEUE_MAX_TASKS", i))
            {
                config.parallel_queue_max_tasks = i;
            }
            if (env_try_int32("PLACER_PARALLEL_RESULT_BUFFER_MAX", i))
            {
                config.parallel_result_buffer_max = i;
            }
            if (env_try_bool("PLACER_LOG_PARALLEL_PROGRESS", b))
            {
                config.log_parallel_progress = b;
            }
            if (env_try_bool("PLACER_TSD_ENABLE", b))
            {
                config.tsd_enable = b;
            }
            if (env_try_int32("PLACER_TSD_MIN_LEN", i))
            {
                config.tsd_min_len = std::max(1, i);
            }
            if (env_try_int32("PLACER_TSD_MAX_LEN", i))
            {
                config.tsd_max_len = std::max(1, i);
            }
            if (env_try_int32("PLACER_TSD_FLANK_WINDOW", i))
            {
                config.tsd_flank_window = std::max(10, i);
            }
            if (env_try_double("PLACER_TSD_BG_P_MAX", v))
            {
                config.tsd_bg_p_max = std::clamp(v, 0.0, 1.0);
            }
            if (env_try_int32("PLACER_GENOTYPE_MIN_DEPTH", i))
            {
                config.genotype_min_depth = std::max(1, i);
            }
            if (env_try_double("PLACER_GENOTYPE_ERROR_RATE", v))
            {
                config.genotype_error_rate = std::clamp(v, 1e-4, 0.25);
            }
            if (env_try_double("PLACER_FINAL_FDR_Q", v))
            {
                config.final_fdr_q = std::clamp(v, 0.0, 1.0);
            }
            if (env_try_int32("PLACER_MIN_FINAL_RAW_CIGAR_INSERT_LEN_BP", i))
            {
                config.min_final_raw_cigar_insert_len_bp = std::max(0, i);
            }
            if (env_try_string("PLACER_FINAL_REPORT_MODE", s))
            {
                FinalReportMode mode = config.final_report_mode;
                if (parse_final_report_mode(s, mode))
                {
                    config.final_report_mode = mode;
                }
            }
            if (env_try_int32("PLACER_EVENT_CONSENSUS_POA_MIN_READS", i))
            {
                config.event_consensus_poa_min_reads = std::max(2, i);
            }
            if (env_try_int32("PLACER_EVENT_CONSENSUS_POA_MAX_READS", i))
            {
                config.event_consensus_poa_max_reads = std::max(2, i);
            }
            config.tsd_max_len = std::max(config.tsd_min_len, config.tsd_max_len);
            config.event_consensus_poa_min_reads = std::min(
                config.event_consensus_poa_min_reads,
                config.event_consensus_poa_max_reads);
        }

        void write_scientific_txt(const PipelineResult &result, const std::string &output_path)
        {
            std::ofstream out(output_path);
            if (!out.is_open())
            {
                throw std::runtime_error("Failed to open output file: " + output_path);
            }
            const bool include_insert_seq = env_flag_enabled("PLACER_DEBUG_OUTPUT_INSERT_SEQ");

            out << "#PLACER streaming pipeline summary\n";
            out << "total_reads\t" << result.total_reads << "\n";
            out << "gate1_passed\t" << result.gate1_passed << "\n";
            out << "processed_bins\t" << result.processed_bins << "\n";
            out << "components\t" << result.built_components << "\n";
            out << "event_consensus_calls\t" << result.event_consensus_calls << "\n";
            out << "genotype_calls\t" << result.genotype_calls << "\n";
            out << "final_pass_calls\t" << result.final_pass_calls << "\n";
            out << "performance_pipeline_wall_seconds\t"
                << result.performance.pipeline_wall_seconds << "\n";
            out << "performance_pipeline_construction_seconds\t"
                << result.performance.pipeline_construction_seconds << "\n";
            out << "performance_bam_stream_wall_seconds\t"
                << result.performance.bam_stream_wall_seconds << "\n";
            out << "performance_bin_processing_wall_seconds\t"
                << result.performance.bin_processing_wall_seconds << "\n";
            out << "performance_finalization_seconds\t"
                << result.performance.finalization_seconds << "\n";
            out << "performance_gate1_preliminary_seconds\t"
                << result.performance.gate1_preliminary_seconds << "\n";
            out << "performance_exact_bin_scan_seconds\t"
                << result.performance.exact_bin_scan_seconds << "\n";
            out << "performance_component_build_seconds\t"
                << result.performance.component_build_seconds << "\n";
            out << "performance_local_fetch_seconds\t"
                << result.performance.local_fetch_seconds << "\n";
            out << "performance_projection_seconds\t"
                << result.performance.projection_seconds << "\n";
            out << "performance_local_component_refresh_seconds\t"
                << result.performance.local_component_refresh_seconds << "\n";
            out << "performance_fragment_extract_seconds\t"
                << result.performance.fragment_extract_seconds << "\n";
            out << "performance_signal_cache_build_seconds\t"
                << result.performance.signal_cache_build_seconds << "\n";
            out << "performance_hypothesis_summary_seconds\t"
                << result.performance.hypothesis_summary_seconds << "\n";
            out << "performance_validator_input_count_seconds\t"
                << result.performance.validator_input_count_seconds << "\n";
            out << "performance_event_consensus_seconds\t"
                << result.performance.event_consensus_seconds << "\n";
            out << "performance_event_segmentation_seconds\t"
                << result.performance.event_segmentation_seconds << "\n";
            out << "performance_te_alignment_seconds\t"
                << result.performance.te_alignment_seconds << "\n";
            out << "performance_te_blast_batches\t"
                << result.performance.te_blast_batches << "\n";
            out << "performance_te_blast_queries\t"
                << result.performance.te_blast_queries << "\n";
            out << "performance_te_blast_cache_hits\t"
                << result.performance.te_blast_cache_hits << "\n";
            out << "performance_te_blast_cache_misses\t"
                << result.performance.te_blast_cache_misses << "\n";
            out << "performance_te_blast_deduplicated_queries\t"
                << result.performance.te_blast_deduplicated_queries << "\n";
            out << "performance_segmentation_edit_distance_calls\t"
                << result.performance.segmentation_edit_distance_calls << "\n";
            out << "performance_segmentation_edit_distance_cache_hits\t"
                << result.performance.segmentation_edit_distance_cache_hits << "\n";
            out << "performance_segmentation_edit_distance_cache_misses\t"
                << result.performance.segmentation_edit_distance_cache_misses << "\n";
            out << "performance_segmentation_seed_bins_total\t"
                << result.performance.segmentation_seed_bins_total << "\n";
            out << "performance_segmentation_paired_searches\t"
                << result.performance.segmentation_paired_searches << "\n";
            out << "performance_segmentation_one_sided_searches\t"
                << result.performance.segmentation_one_sided_searches << "\n";
            out << "performance_segmentation_endpoint_slack_searches\t"
                << result.performance.segmentation_endpoint_slack_searches << "\n";
            out << "performance_segmentation_reverse_complement_retries\t"
                << result.performance.segmentation_reverse_complement_retries << "\n";
            out << "performance_segmentation_cache_hits\t"
                << result.performance.segmentation_cache_hits << "\n";
            out << "performance_segmentation_cache_misses\t"
                << result.performance.segmentation_cache_misses << "\n";
            out << "schema_version\t1.4.0\n";

            out << "\n#chrom\tpos\tbp_left\tbp_right\tte\tfamily\tsubfamily\tfamily_status"
                << "\tstrand\tinsert_len"
                << "\tsupport_reads\talt_struct_reads\traw_cigar_insert_reads"
                << "\tmax_raw_cigar_insert_len\tref_span_reads\tlow_mapq_ref_span_reads"
                << "\tgt\taf\tgq"
                << "\tbest_te_identity\tbest_te_query_coverage\tcross_family_margin"
                << "\tte_sequence_model_label\tte_sequence_model_score"
                << "\tte_sequence_model_gc\tte_sequence_model_entropy"
                << "\tte_sequence_model_tandem_fraction"
                << "\tte_sequence_model_low_complexity_fraction"
                << "\tte_sequence_model_jsd_k5\tte_sequence_model_jsd_k6"
                << "\tte_sequence_model_k9_containment"
                << "\tte_annotation_confidence\tte_annotation_class\tte_annotation_order"
                << "\tte_annotation_intervals"
                << "\tte_annotation_residual_fraction\tte_annotation_masked_fraction"
                << "\ttsd_type\ttsd_len"
                << "\tleft_flank_align_len\tright_flank_align_len\tconsensus_len";
            if (include_insert_seq)
            {
                out << "\tinsert_seq";
            }
            out << "\tqc\tbest_explanation\texplanation_residual\texplanation_path"
                << "\tte_structure_path\tte_structure_log_evidence"
                << "\tnonte_structure_log_evidence"
                << "\tartifact_structure_log_evidence"
                << "\tte_structure_path_confidence"
                << "\tpolyA_posterior\ttransduction_posterior"
                << "\tte_core_coverage\tunexplained_high_complexity_bp"
                << "\tte_posterior\tnon_te_posterior\tartifact_posterior"
                << "\tte_vs_artifact_log_odds\tte_vs_non_te_log_odds\tposterior_qc"
                << "\tlatent_mechanism\tfamily_activity_prior\tlfdr"
                << "\tworst_case_lfdr\tlfdr_qc"
                << "\tmechanistic_lower_log_bf_te_vs_artifact"
                << "\tmechanistic_lower_log_bf_te_vs_non_te"
                << "\tmechanistic_ref_conflict_signal"
                << "\tmechanistic_ambiguity_width"
                << "\tmechanistic_blocks"
                << "\trobust_mechanistic_lfdr"
                << "\trobust_mechanistic_worst_case_lfdr"
                << "\trobust_mechanistic_qc"
                << "\tconformal_null_p"
                << "\tconformal_by_threshold"
                << "\tconformal_dominated_nulls"
                << "\tconformal_null_count"
                << "\tconformal_qc"
                << "\tbp_ci_width\tbp_posterior_entropy"
                << "\tte_consensus_start\tte_consensus_end\n";
            for (const auto &call : result.final_calls)
            {
                out << call.chrom << "\t"
                    << call.pos << "\t"
                    << call.bp_left << "\t"
                    << call.bp_right << "\t"
                    << (call.te_name.empty() ? "NA" : call.te_name) << "\t"
                    << call.family << "\t"
                    << call.subfamily << "\t"
                    << (call.family_committed ? "COMMITTED" : "ABSTAINED") << "\t"
                    << call.strand << "\t"
                    << call.insert_len << "\t"
                    << call.support_reads << "\t"
                    << call.alt_struct_reads << "\t"
                    << call.raw_cigar_insert_reads << "\t"
                    << call.max_raw_cigar_insert_len << "\t"
                    << call.ref_span_reads << "\t"
                    << call.low_mapq_ref_span_reads << "\t"
                    << call.genotype << "\t"
                    << call.af << "\t"
                    << call.gq << "\t"
                    << call.best_te_identity << "\t"
                    << call.best_te_query_coverage << "\t"
                    << call.cross_family_margin << "\t"
                    << call.te_sequence_model_label << "\t"
                    << call.te_sequence_model_score << "\t"
                    << call.te_sequence_model_gc << "\t"
                    << call.te_sequence_model_entropy << "\t"
                    << call.te_sequence_model_tandem_fraction << "\t"
                    << call.te_sequence_model_low_complexity_fraction << "\t"
                    << call.te_sequence_model_jsd_k5 << "\t"
                    << call.te_sequence_model_jsd_k6 << "\t"
                    << call.te_sequence_model_k9_containment << "\t"
                    << call.te_annotation_confidence << "\t"
                    << call.te_annotation_class << "\t"
                    << call.te_annotation_order << "\t"
                    << call.te_annotation_intervals << "\t"
                    << call.te_annotation_residual_fraction << "\t"
                    << call.te_annotation_masked_fraction << "\t"
                    << call.tsd_type << "\t"
                    << call.tsd_len << "\t"
                    << call.left_flank_align_len << "\t"
                    << call.right_flank_align_len << "\t"
                    << call.event_consensus_len << "\t";
                if (include_insert_seq)
                {
                    out << call.insert_seq << "\t";
                }
                out << call.final_qc << "\t"
                    << call.best_explanation << "\t"
                    << call.explanation_residual << "\t"
                    << call.explanation_path << "\t"
                    << call.te_structure_path << "\t"
                    << call.te_structure_log_evidence << "\t"
                    << call.nonte_structure_log_evidence << "\t"
                    << call.artifact_structure_log_evidence << "\t"
                    << call.te_structure_path_confidence << "\t"
                    << call.polyA_posterior << "\t"
                    << call.transduction_posterior << "\t"
                    << call.te_core_coverage << "\t"
                    << call.unexplained_high_complexity_bp << "\t"
                    << call.te_posterior << "\t"
                    << call.non_te_posterior << "\t"
                    << call.artifact_posterior << "\t"
                    << call.te_vs_artifact_log_odds << "\t"
                    << call.te_vs_non_te_log_odds << "\t"
                    << call.posterior_qc << "\t"
                    << call.latent_mechanism << "\t"
                    << call.family_activity_prior << "\t"
                    << call.lfdr << "\t"
                    << call.worst_case_lfdr << "\t"
                    << call.lfdr_qc << "\t"
                    << call.mechanistic_lower_log_bf_te_vs_artifact << "\t"
                    << call.mechanistic_lower_log_bf_te_vs_non_te << "\t"
                    << call.mechanistic_ref_conflict_signal << "\t"
                    << call.mechanistic_ambiguity_width << "\t"
                    << call.mechanistic_blocks << "\t"
                    << call.robust_mechanistic_lfdr << "\t"
                    << call.robust_mechanistic_worst_case_lfdr << "\t"
                    << call.robust_mechanistic_qc << "\t"
                    << call.conformal_null_p << "\t"
                    << call.conformal_by_threshold << "\t"
                    << call.conformal_dominated_nulls << "\t"
                    << call.conformal_null_count << "\t"
                    << call.conformal_qc << "\t"
                    << call.bp_ci_width << "\t"
                    << call.bp_posterior_entropy << "\t"
                    << call.te_consensus_start << "\t"
                    << call.te_consensus_end << "\n";
            }
        }

        void write_evidence_ledger_tsv(const PipelineResult &result, const std::string &output_path)
        {
            std::ofstream out(output_path);
            if (!out.is_open())
            {
                throw std::runtime_error("Failed to open output file: " + output_path);
            }
            const bool include_support_qnames =
                env_flag_enabled("PLACER_DEBUG_OUTPUT_SUPPORT_QNAMES");
            const bool include_insert_seq = include_candidate_insert_seq_debug_output();
            out << "chrom\tpos\tbp_left\tbp_right\tcoverage_left\tcoverage_right\tfamily\tsubfamily"
                << "\tfinal_qc\tposterior_qc\tlfdr_qc\tcandidate_retention_reason"
                << "\talt_struct_reads\talt_split_reads\talt_indel_reads"
                << "\talt_left_clip_reads\talt_right_clip_reads"
                << "\traw_cigar_insert_reads\tmax_raw_cigar_insert_len"
                << "\tfull_context_input_reads\tpartial_context_input_reads"
                << "\tleft_anchor_input_reads\tright_anchor_input_reads"
                << "\tinput_event_reads\tevent_consensus_len";
            if (include_insert_seq)
            {
                out << "\tinsert_seq";
            }
            out
                << "\tleft_flank_align_len\tright_flank_align_len"
                << "\tref_span_reads\tsupport_qname_count";
            if (include_support_qnames)
            {
                out << "\tsupport_qnames";
            }
            out << "\tbest_te_identity\tbest_te_query_coverage\tcross_family_margin"
                << "\tte_structure_path\tte_structure_log_evidence"
                << "\tnonte_structure_log_evidence"
                << "\tartifact_structure_log_evidence"
                << "\tte_structure_path_confidence"
                << "\tpolyA_posterior\ttransduction_posterior"
                << "\tte_core_coverage\tunexplained_high_complexity_bp"
                << "\tte_posterior\tnon_te_posterior\tartifact_posterior"
                << "\tlfdr\tworst_case_lfdr"
                << "\tmechanistic_lower_log_bf_te_vs_artifact"
                << "\tmechanistic_lower_log_bf_te_vs_non_te"
                << "\tmechanistic_ref_conflict_signal"
                << "\tmechanistic_ambiguity_width"
                << "\tmechanistic_blocks"
                << "\trobust_mechanistic_lfdr"
                << "\trobust_mechanistic_worst_case_lfdr"
                << "\trobust_mechanistic_qc"
                << "\tconformal_null_p"
                << "\tconformal_by_threshold"
                << "\tconformal_dominated_nulls"
                << "\tconformal_null_count"
                << "\tconformal_qc\n";
            for (const auto &row : result.evidence_ledger)
            {
                out << row.chrom << "\t"
                    << row.pos << "\t"
                    << row.bp_left << "\t"
                    << row.bp_right << "\t"
                    << row.coverage_left << "\t"
                    << row.coverage_right << "\t"
                    << row.family << "\t"
                    << row.subfamily << "\t"
                    << row.final_qc << "\t"
                    << row.posterior_qc << "\t"
                    << row.lfdr_qc << "\t"
                    << row.candidate_retention_reason << "\t"
                    << row.alt_struct_reads << "\t"
                    << row.alt_split_reads << "\t"
                    << row.alt_indel_reads << "\t"
                    << row.alt_left_clip_reads << "\t"
                    << row.alt_right_clip_reads << "\t"
                    << row.raw_cigar_insert_reads << "\t"
                    << row.max_raw_cigar_insert_len << "\t"
                    << row.full_context_input_reads << "\t"
                    << row.partial_context_input_reads << "\t"
                    << row.left_anchor_input_reads << "\t"
                    << row.right_anchor_input_reads << "\t"
                    << row.input_event_reads << "\t"
                    << row.event_consensus_len;
                if (include_insert_seq)
                {
                    out << "\t" << row.insert_seq;
                }
                out << "\t"
                    << row.left_flank_align_len << "\t"
                    << row.right_flank_align_len << "\t"
                    << row.ref_span_reads << "\t"
                    << row.support_qnames.size();
                if (include_support_qnames)
                {
                    out << "\t" << serialize_support_qnames(row.support_qnames);
                }
                out << "\t"
                    << row.best_te_identity << "\t"
                    << row.best_te_query_coverage << "\t"
                    << row.cross_family_margin << "\t"
                    << row.te_structure_path << "\t"
                    << row.te_structure_log_evidence << "\t"
                    << row.nonte_structure_log_evidence << "\t"
                    << row.artifact_structure_log_evidence << "\t"
                    << row.te_structure_path_confidence << "\t"
                    << row.polyA_posterior << "\t"
                    << row.transduction_posterior << "\t"
                    << row.te_core_coverage << "\t"
                    << row.unexplained_high_complexity_bp << "\t"
                    << row.te_posterior << "\t"
                    << row.non_te_posterior << "\t"
                    << row.artifact_posterior << "\t"
                    << row.lfdr << "\t"
                    << row.worst_case_lfdr << "\t"
                    << row.mechanistic_lower_log_bf_te_vs_artifact << "\t"
                    << row.mechanistic_lower_log_bf_te_vs_non_te << "\t"
                    << row.mechanistic_ref_conflict_signal << "\t"
                    << row.mechanistic_ambiguity_width << "\t"
                    << row.mechanistic_blocks << "\t"
                    << row.robust_mechanistic_lfdr << "\t"
                    << row.robust_mechanistic_worst_case_lfdr << "\t"
                    << row.robust_mechanistic_qc << "\t"
                    << row.conformal_null_p << "\t"
                    << row.conformal_by_threshold << "\t"
                    << row.conformal_dominated_nulls << "\t"
                    << row.conformal_null_count << "\t"
                    << row.conformal_qc << "\n";
            }
        }

        void write_coverage_candidates_tsv(const PipelineResult &result, const std::string &output_path)
        {
            std::ofstream out(output_path);
            if (!out.is_open())
            {
                throw std::runtime_error("Failed to open output file: " + output_path);
            }
            const bool include_support_qnames =
                env_flag_enabled("PLACER_DEBUG_OUTPUT_SUPPORT_QNAMES");
            const bool include_insert_seq = include_candidate_insert_seq_debug_output();
            out << "chrom\tcoverage_left\tcoverage_right\tpos\tbp_left\tbp_right"
                << "\tfamily\tsubfamily\tfinal_qc\tcandidate_retention_reason"
                << "\talt_struct_reads\talt_split_reads\talt_indel_reads"
                << "\talt_left_clip_reads\talt_right_clip_reads"
                << "\traw_cigar_insert_reads\tmax_raw_cigar_insert_len"
                << "\tfull_context_input_reads\tpartial_context_input_reads"
                << "\tleft_anchor_input_reads\tright_anchor_input_reads"
                << "\tinput_event_reads\tevent_consensus_len";
            if (include_insert_seq)
            {
                out << "\tinsert_seq";
            }
            out
                << "\tleft_flank_align_len\tright_flank_align_len"
                << "\tref_span_reads\tsupport_qname_count";
            if (include_support_qnames)
            {
                out << "\tsupport_qnames";
            }
            out << "\tbest_te_identity\tbest_te_query_coverage\tcross_family_margin"
                << "\tconformal_null_p\tconformal_qc\n";
            for (const auto &row : result.evidence_ledger)
            {
                if (row.coverage_left < 0 || row.coverage_right < 0)
                {
                    continue;
                }
                out << row.chrom << "\t"
                    << std::min(row.coverage_left, row.coverage_right) << "\t"
                    << std::max(row.coverage_left, row.coverage_right) << "\t"
                    << row.pos << "\t"
                    << row.bp_left << "\t"
                    << row.bp_right << "\t"
                    << row.family << "\t"
                    << row.subfamily << "\t"
                    << row.final_qc << "\t"
                    << row.candidate_retention_reason << "\t"
                    << row.alt_struct_reads << "\t"
                    << row.alt_split_reads << "\t"
                    << row.alt_indel_reads << "\t"
                    << row.alt_left_clip_reads << "\t"
                    << row.alt_right_clip_reads << "\t"
                    << row.raw_cigar_insert_reads << "\t"
                    << row.max_raw_cigar_insert_len << "\t"
                    << row.full_context_input_reads << "\t"
                    << row.partial_context_input_reads << "\t"
                    << row.left_anchor_input_reads << "\t"
                    << row.right_anchor_input_reads << "\t"
                    << row.input_event_reads << "\t"
                    << row.event_consensus_len;
                if (include_insert_seq)
                {
                    out << "\t" << row.insert_seq;
                }
                out << "\t"
                    << row.left_flank_align_len << "\t"
                    << row.right_flank_align_len << "\t"
                    << row.ref_span_reads << "\t"
                    << row.support_qnames.size();
                if (include_support_qnames)
                {
                    out << "\t" << serialize_support_qnames(row.support_qnames);
                }
                out << "\t"
                    << row.best_te_identity << "\t"
                    << row.best_te_query_coverage << "\t"
                    << row.cross_family_margin << "\t"
                    << row.conformal_null_p << "\t"
                    << row.conformal_qc << "\n";
            }
        }

        bool is_final_te_qc(const std::string &qc)
        {
            return qc.rfind("PASS_TE", 0) == 0;
        }

        bool has_mechanistic_certificate(const EvidenceLedgerRow &row)
        {
            return row.mechanistic_blocks != "NA";
        }

        double mechanistic_null_score(double bf_artifact, double bf_non_te, double ambiguity_width)
        {
            return std::min(bf_artifact, bf_non_te) - std::max(0.0, ambiguity_width);
        }

        bool is_empirical_null_control(const EvidenceLedgerRow &row)
        {
            if (!has_mechanistic_certificate(row) || is_final_te_qc(row.final_qc))
            {
                return false;
            }
            return true;
        }

        std::string null_control_kind(const EvidenceLedgerRow &row)
        {
            if (row.final_qc == "PASS_NONTE_INSERTION")
            {
                return "non_te_insert_control";
            }
            if (row.final_qc == "REFERENCE_OR_ARTIFACT")
            {
                return "artifact_reference_control";
            }
            if (row.final_qc == "TE_AMBIGUOUS")
            {
                return "ambiguous_te_control";
            }
            return "non_final_evidence_control";
        }

        void write_null_controls_tsv(const PipelineResult &result, const std::string &output_path)
        {
            EmpiricalNullTail null_tail;
            for (const auto &row : result.evidence_ledger)
            {
                if (!is_empirical_null_control(row))
                {
                    continue;
                }
                const double score = mechanistic_null_score(
                    row.mechanistic_lower_log_bf_te_vs_artifact,
                    row.mechanistic_lower_log_bf_te_vs_non_te,
                    row.mechanistic_ambiguity_width);
                if (std::isfinite(score))
                {
                    null_tail.add(score);
                }
            }

            std::ofstream out(output_path);
            if (!out.is_open())
            {
                throw std::runtime_error("Failed to open output file: " + output_path);
            }
            out << "row_type\tcontrol_kind\tchrom\tpos\tbp_left\tbp_right\tfamily\tsubfamily"
                << "\tfinal_qc\tcandidate_retention_reason\tcontrol_score"
                << "\tempirical_null_upper_tail_p\tnull_control_count"
                << "\trobust_mechanistic_worst_case_lfdr\tmechanistic_ref_conflict_signal"
                << "\tconformal_null_p\tconformal_by_threshold"
                << "\tconformal_dominated_nulls\tconformal_null_count\tconformal_qc\n";

            for (const auto &row : result.evidence_ledger)
            {
                if (!is_empirical_null_control(row))
                {
                    continue;
                }
                const double score = mechanistic_null_score(
                    row.mechanistic_lower_log_bf_te_vs_artifact,
                    row.mechanistic_lower_log_bf_te_vs_non_te,
                    row.mechanistic_ambiguity_width);
                out << "null_control\t"
                    << null_control_kind(row) << "\t"
                    << row.chrom << "\t"
                    << row.pos << "\t"
                    << row.bp_left << "\t"
                    << row.bp_right << "\t"
                    << row.family << "\t"
                    << row.subfamily << "\t"
                    << row.final_qc << "\t"
                    << row.candidate_retention_reason << "\t"
                    << score << "\t"
                    << null_tail.upper_tail_p(score) << "\t"
                    << null_tail.size() << "\t"
                    << row.robust_mechanistic_worst_case_lfdr << "\t"
                    << row.mechanistic_ref_conflict_signal << "\t"
                    << row.conformal_null_p << "\t"
                    << row.conformal_by_threshold << "\t"
                    << row.conformal_dominated_nulls << "\t"
                    << row.conformal_null_count << "\t"
                    << row.conformal_qc << "\n";
            }

            for (const auto &call : result.final_calls)
            {
                const double score = mechanistic_null_score(
                    call.mechanistic_lower_log_bf_te_vs_artifact,
                    call.mechanistic_lower_log_bf_te_vs_non_te,
                    call.mechanistic_ambiguity_width);
                out << "final_call_tail_query\tfinal_call\t"
                    << call.chrom << "\t"
                    << call.pos << "\t"
                    << call.bp_left << "\t"
                    << call.bp_right << "\t"
                    << call.family << "\t"
                    << call.subfamily << "\t"
                    << call.final_qc << "\t"
                    << "FINAL_PASS\t"
                    << score << "\t"
                    << null_tail.upper_tail_p(score) << "\t"
                    << null_tail.size() << "\t"
                    << call.robust_mechanistic_worst_case_lfdr << "\t"
                    << call.mechanistic_ref_conflict_signal << "\t"
                    << call.conformal_null_p << "\t"
                    << call.conformal_by_threshold << "\t"
                    << call.conformal_dominated_nulls << "\t"
                    << call.conformal_null_count << "\t"
                    << call.conformal_qc << "\n";
            }
        }

        int run_pipeline_once(const PipelineConfig &config)
        {
            try
            {
                const std::string scientific_path = safe_absolute_path("scientific.txt");
                const std::string evidence_ledger_path = safe_absolute_path("evidence_ledger.tsv");
                const std::string coverage_candidates_path = safe_absolute_path("coverage_candidates.tsv");
                const std::string null_controls_path = safe_absolute_path("null_controls.tsv");
                std::cerr << "[PLACER] run started\n"
                          << "  cwd=" << safe_current_path() << "\n"
                          << "  bam=" << config.bam_path << "\n"
                          << "  bam_region=" << (config.bam_region_scope.enabled ? config.bam_region_scope.chrom + ":" + std::to_string(config.bam_region_scope.start + 1) + "-" + (config.bam_region_scope.end > 0 ? std::to_string(config.bam_region_scope.end) : std::string("end")) : std::string("ALL")) << "\n"
                          << "  reference=" << config.reference_fasta_path << "\n"
                          << "  te_fasta=" << (config.te_fasta_path.empty() ? "NA" : config.te_fasta_path) << "\n"
                          << "  mode=" << (config.enable_parallel ? "parallel" : "streaming") << "\n"
                          << "  bam_threads=" << config.bam_threads << "\n"
                          << "  parallel_workers=" << config.parallel_workers << "\n"
                          << "  parallel_queue_max_tasks=" << config.parallel_queue_max_tasks << "\n"
                          << "  parallel_result_buffer_max=" << config.parallel_result_buffer_max << "\n"
                          << "  log_parallel_progress=" << (config.log_parallel_progress ? 1 : 0) << "\n"
                          << "  progress_interval=" << config.progress_interval << "\n"
                          << "  log_stage_bins=" << (config.log_stage_bins ? 1 : 0) << "\n"
                          << "  log_stage_components=" << (config.log_stage_components ? 1 : 0) << "\n"
                          << "  min_final_raw_cigar_insert_len_bp="
                          << config.min_final_raw_cigar_insert_len_bp << "\n"
                          << "  final_report_mode="
                          << final_report_mode_name(config.final_report_mode) << "\n"
                          << "  scientific_txt=" << scientific_path << "\n"
                          << "  evidence_ledger_tsv=" << evidence_ledger_path << "\n"
                          << "  coverage_candidates_tsv=" << coverage_candidates_path << "\n"
                          << "  null_controls_tsv=" << null_controls_path << std::endl;

                const auto construction_started = std::chrono::steady_clock::now();
                auto pipeline = build_default_pipeline(config);
                const double pipeline_construction_seconds = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - construction_started).count();
                PipelineResult result = pipeline->run();
                result.performance.pipeline_construction_seconds =
                    pipeline_construction_seconds;

                std::cerr << "[PLACER] pipeline finished\n"
                          << "  mode=" << (config.enable_parallel ? "parallel" : "streaming") << "\n"
                          << "  total_reads=" << result.total_reads << "\n"
                          << "  gate1_passed=" << result.gate1_passed << "\n"
                          << "  processed_bins=" << result.processed_bins << "\n"
                          << "  components=" << result.built_components << "\n"
                          << "  event_consensus_calls=" << result.event_consensus_calls << "\n"
                          << "  genotype_calls=" << result.genotype_calls << "\n"
                          << "  final_pass_calls=" << result.final_pass_calls << "\n"
                          << "  performance_pipeline_wall_seconds="
                          << result.performance.pipeline_wall_seconds << "\n"
                          << "  performance_pipeline_construction_seconds="
                          << result.performance.pipeline_construction_seconds << "\n"
                          << "  performance_bam_stream_wall_seconds="
                          << result.performance.bam_stream_wall_seconds << "\n"
                          << "  performance_bin_processing_wall_seconds="
                          << result.performance.bin_processing_wall_seconds << "\n"
                          << "  performance_finalization_seconds="
                          << result.performance.finalization_seconds << "\n"
                          << "  performance_gate1_preliminary_seconds="
                          << result.performance.gate1_preliminary_seconds << "\n"
                          << "  performance_exact_bin_scan_seconds="
                          << result.performance.exact_bin_scan_seconds << "\n"
                          << "  performance_component_build_seconds="
                          << result.performance.component_build_seconds << "\n"
                          << "  performance_local_fetch_seconds="
                          << result.performance.local_fetch_seconds << "\n"
                          << "  performance_projection_seconds="
                          << result.performance.projection_seconds << "\n"
                          << "  performance_local_component_refresh_seconds="
                          << result.performance.local_component_refresh_seconds << "\n"
                          << "  performance_fragment_extract_seconds="
                          << result.performance.fragment_extract_seconds << "\n"
                          << "  performance_signal_cache_build_seconds="
                          << result.performance.signal_cache_build_seconds << "\n"
                          << "  performance_hypothesis_summary_seconds="
                          << result.performance.hypothesis_summary_seconds << "\n"
                          << "  performance_validator_input_count_seconds="
                          << result.performance.validator_input_count_seconds << "\n"
                          << "  performance_event_consensus_seconds="
                          << result.performance.event_consensus_seconds << "\n"
                          << "  performance_event_segmentation_seconds="
                          << result.performance.event_segmentation_seconds << "\n"
                          << "  performance_te_alignment_seconds="
                          << result.performance.te_alignment_seconds << "\n"
                          << "  performance_te_blast_batches="
                          << result.performance.te_blast_batches << "\n"
                          << "  performance_te_blast_queries="
                          << result.performance.te_blast_queries << "\n"
                          << "  performance_te_blast_cache_hits="
                          << result.performance.te_blast_cache_hits << "\n"
                          << "  performance_te_blast_cache_misses="
                          << result.performance.te_blast_cache_misses << "\n"
                          << "  performance_te_blast_deduplicated_queries="
                          << result.performance.te_blast_deduplicated_queries << "\n"
                          << "  performance_segmentation_edit_distance_calls="
                          << result.performance.segmentation_edit_distance_calls << "\n"
                          << "  performance_segmentation_edit_distance_cache_hits="
                          << result.performance.segmentation_edit_distance_cache_hits << "\n"
                          << "  performance_segmentation_edit_distance_cache_misses="
                          << result.performance.segmentation_edit_distance_cache_misses << "\n"
                          << "  performance_segmentation_seed_bins_total="
                          << result.performance.segmentation_seed_bins_total << "\n"
                          << "  performance_segmentation_paired_searches="
                          << result.performance.segmentation_paired_searches << "\n"
                          << "  performance_segmentation_one_sided_searches="
                          << result.performance.segmentation_one_sided_searches << "\n"
                          << "  performance_segmentation_endpoint_slack_searches="
                          << result.performance.segmentation_endpoint_slack_searches << "\n"
                          << "  performance_segmentation_reverse_complement_retries="
                          << result.performance.segmentation_reverse_complement_retries << "\n"
                          << "  performance_segmentation_cache_hits="
                          << result.performance.segmentation_cache_hits << "\n"
                          << "  performance_segmentation_cache_misses="
                          << result.performance.segmentation_cache_misses << std::endl;

                double scientific_write_seconds = 0.0;
                double evidence_ledger_write_seconds = 0.0;
                double coverage_candidates_write_seconds = 0.0;
                double null_controls_write_seconds = 0.0;
                {
                    const auto started = std::chrono::steady_clock::now();
                    write_scientific_txt(result, "scientific.txt");
                    scientific_write_seconds = std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - started).count();
                }
                {
                    const auto started = std::chrono::steady_clock::now();
                    write_evidence_ledger_tsv(result, "evidence_ledger.tsv");
                    evidence_ledger_write_seconds = std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - started).count();
                }
                {
                    const auto started = std::chrono::steady_clock::now();
                    write_coverage_candidates_tsv(result, "coverage_candidates.tsv");
                    coverage_candidates_write_seconds = std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - started).count();
                }
                {
                    const auto started = std::chrono::steady_clock::now();
                    write_null_controls_tsv(result, "null_controls.tsv");
                    null_controls_write_seconds = std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - started).count();
                }
                std::cerr << "[PLACER] output write performance\n"
                          << "  performance_scientific_write_seconds="
                          << scientific_write_seconds << "\n"
                          << "  performance_evidence_ledger_write_seconds="
                          << evidence_ledger_write_seconds << "\n"
                          << "  performance_coverage_candidates_write_seconds="
                          << coverage_candidates_write_seconds << "\n"
                          << "  performance_null_controls_write_seconds="
                          << null_controls_write_seconds << std::endl;
                std::cerr << "[PLACER] wrote scientific.txt path=" << scientific_path << std::endl;
                std::cerr << "[PLACER] wrote evidence_ledger.tsv path="
                          << evidence_ledger_path << std::endl;
                std::cerr << "[PLACER] wrote coverage_candidates.tsv path="
                          << coverage_candidates_path << std::endl;
                std::cerr << "[PLACER] wrote null_controls.tsv path="
                          << null_controls_path << std::endl;
                return 0;
            }
            catch (const std::exception &ex)
            {
                std::cerr << "[PLACER] fatal: " << ex.what() << std::endl;
                return 2;
            }
        }

        struct BatchRunSpec
        {
            std::string sample_id;
            std::string bam_path;
            std::string run_dir;
        };

        std::vector<std::string> split_tab_line(const std::string &line)
        {
            std::vector<std::string> fields;
            std::string field;
            std::istringstream in(line);
            while (std::getline(in, field, '\t'))
            {
                fields.push_back(field);
            }
            if (!line.empty() && line.back() == '\t')
            {
                fields.emplace_back();
            }
            return fields;
        }

        int find_column_or_throw(const std::vector<std::string> &header, const std::string &name)
        {
            for (size_t i = 0; i < header.size(); ++i)
            {
                if (header[i] == name)
                {
                    return static_cast<int>(i);
                }
            }
            throw std::runtime_error("batch manifest missing column: " + name);
        }

        std::vector<BatchRunSpec> load_batch_manifest(const std::string &path)
        {
            std::ifstream in(path);
            if (!in.is_open())
            {
                throw std::runtime_error("failed to open batch manifest: " + path);
            }

            std::string line;
            if (!std::getline(in, line))
            {
                throw std::runtime_error("empty batch manifest: " + path);
            }
            const std::vector<std::string> header = split_tab_line(line);
            const int sample_col = find_column_or_throw(header, "sample_id");
            const int bam_col = find_column_or_throw(header, "bam");
            const int run_dir_col = find_column_or_throw(header, "run_dir");
            const int required_cols = 1 + std::max({sample_col, bam_col, run_dir_col});

            std::vector<BatchRunSpec> runs;
            int line_no = 1;
            while (std::getline(in, line))
            {
                ++line_no;
                if (line.empty())
                {
                    continue;
                }
                const std::vector<std::string> fields = split_tab_line(line);
                if (static_cast<int>(fields.size()) < required_cols)
                {
                    throw std::runtime_error("malformed batch manifest row " + std::to_string(line_no));
                }
                BatchRunSpec spec;
                spec.sample_id = fields[static_cast<size_t>(sample_col)];
                spec.bam_path = safe_absolute_path(fields[static_cast<size_t>(bam_col)]);
                spec.run_dir = safe_absolute_path(fields[static_cast<size_t>(run_dir_col)]);
                if (spec.sample_id.empty() || spec.bam_path.empty() || spec.run_dir.empty())
                {
                    throw std::runtime_error("empty required field in batch manifest row " + std::to_string(line_no));
                }
                runs.push_back(std::move(spec));
            }
            if (runs.empty())
            {
                throw std::runtime_error("batch manifest contains no runs: " + path);
            }
            return runs;
        }

        int run_batch_cli(int argc, char **argv)
        {
            if (argc != 5)
            {
                print_usage();
                return 1;
            }

            PipelineConfig base_config;
            base_config.reference_fasta_path = safe_absolute_path(argv[3]);
            base_config.te_fasta_path = safe_absolute_path(argv[4]);
            apply_environment_config(base_config);

            std::vector<BatchRunSpec> runs;
            try
            {
                runs = load_batch_manifest(argv[2]);
            }
            catch (const std::exception &ex)
            {
                std::cerr << "[PLACER][batch] fatal: " << ex.what() << std::endl;
                return 2;
            }

            std::cout << "sample_id\texit_code\telapsed_s\n";
            int failures = 0;
            const std::filesystem::path original_cwd = std::filesystem::current_path();
            for (const auto &run : runs)
            {
                const auto started = std::chrono::steady_clock::now();
                int exit_code = 2;
                try
                {
                    std::filesystem::create_directories(run.run_dir);
                    std::ofstream log(std::filesystem::path(run.run_dir) / "run.log", std::ios::app);
                    if (!log.is_open())
                    {
                        throw std::runtime_error("failed to open run log: " + run.run_dir);
                    }
                    log << "$ " << argv[0] << " batch " << argv[2] << " " << argv[3] << " " << argv[4]
                        << " sample_id=" << run.sample_id << "\n";

                    std::streambuf *old_cerr = std::cerr.rdbuf(log.rdbuf());
                    try
                    {
                        std::filesystem::current_path(run.run_dir);
                        PipelineConfig config = base_config;
                        config.bam_path = run.bam_path;
                        exit_code = run_pipeline_once(config);
                        std::filesystem::current_path(original_cwd);
                    }
                    catch (...)
                    {
                        std::filesystem::current_path(original_cwd);
                        std::cerr.rdbuf(old_cerr);
                        throw;
                    }
                    std::cerr.rdbuf(old_cerr);
                }
                catch (const std::exception &ex)
                {
                    std::filesystem::current_path(original_cwd);
                    std::cerr << "[PLACER][batch] sample_id=" << run.sample_id
                              << " fatal: " << ex.what() << std::endl;
                    exit_code = 2;
                }

                if (exit_code != 0)
                {
                    ++failures;
                }
                const auto elapsed = std::chrono::duration<double>(
                                         std::chrono::steady_clock::now() - started)
                                         .count();
                std::cout << run.sample_id << "\t" << exit_code << "\t" << elapsed << "\n";
            }
            return failures == 0 ? 0 : 2;
        }

    } // namespace
} // namespace placer

int main(int argc, char **argv)
{
    if (argc >= 2 && std::string(argv[1]) == "denovo")
    {
        return placer::run_denovo_cli(argc - 1, argv + 1);
    }
    if (argc >= 2 && std::string(argv[1]) == "batch")
    {
        return placer::run_batch_cli(argc, argv);
    }

    placer::BamRegionScope bam_region_scope;
    bool final_fdr_q_set = false;
    double final_fdr_q = 0.10;
    bool final_report_mode_set = false;
    placer::FinalReportMode final_report_mode = placer::FinalReportMode::TeCalibrated;
    bool min_final_raw_cigar_insert_len_bp_set = false;
    int32_t min_final_raw_cigar_insert_len_bp = 50;
    int argi = 1;
    while (argi < argc)
    {
        const std::string arg = argv[argi];
        if (arg == "--region")
        {
            if (bam_region_scope.enabled || argi + 1 >= argc)
            {
                placer::print_usage();
                return 1;
            }
            try
            {
                bam_region_scope = placer::parse_region_scope_or_throw(argv[argi + 1]);
            }
            catch (const std::exception &ex)
            {
                std::cerr << "[PLACER] invalid --region: " << ex.what() << std::endl;
                return 1;
            }
            argi += 2;
            continue;
        }
        if (arg == "--final-fdr-q")
        {
            if (argi + 1 >= argc)
            {
                placer::print_usage();
                return 1;
            }
            try
            {
                final_fdr_q = std::stod(argv[argi + 1]);
            }
            catch (const std::exception &)
            {
                std::cerr << "[PLACER] invalid --final-fdr-q: " << argv[argi + 1] << std::endl;
                return 1;
            }
            if (!std::isfinite(final_fdr_q) || final_fdr_q < 0.0 || final_fdr_q > 1.0)
            {
                std::cerr << "[PLACER] --final-fdr-q must be in [0,1]" << std::endl;
                return 1;
            }
            final_fdr_q_set = true;
            argi += 2;
            continue;
        }
        if (arg == "--final-report-mode")
        {
            if (argi + 1 >= argc)
            {
                placer::print_usage();
                return 1;
            }
            if (!placer::parse_final_report_mode(argv[argi + 1], final_report_mode))
            {
                std::cerr << "[PLACER] --final-report-mode must be legacy or te-calibrated"
                          << std::endl;
                return 1;
            }
            final_report_mode_set = true;
            argi += 2;
            continue;
        }
        if (arg == "--min-final-raw-cigar-insert-len-bp")
        {
            if (argi + 1 >= argc)
            {
                placer::print_usage();
                return 1;
            }
            try
            {
                size_t parsed_chars = 0;
                min_final_raw_cigar_insert_len_bp = std::stoi(argv[argi + 1], &parsed_chars);
                if (parsed_chars != std::string(argv[argi + 1]).size())
                {
                    throw std::invalid_argument("trailing characters");
                }
            }
            catch (const std::exception &)
            {
                std::cerr << "[PLACER] invalid --min-final-raw-cigar-insert-len-bp: "
                          << argv[argi + 1] << std::endl;
                return 1;
            }
            if (min_final_raw_cigar_insert_len_bp < 0)
            {
                std::cerr << "[PLACER] --min-final-raw-cigar-insert-len-bp must be >= 0" << std::endl;
                return 1;
            }
            min_final_raw_cigar_insert_len_bp_set = true;
            argi += 2;
            continue;
        }
        if (arg.rfind("--", 0) == 0)
        {
            placer::print_usage();
            return 1;
        }
        break;
    }

    if (argc - argi != 3)
    {
        placer::print_usage();
        return 1;
    }

    placer::PipelineConfig config;
    config.bam_path = argv[argi];
    config.reference_fasta_path = argv[argi + 1];
    config.te_fasta_path = argv[argi + 2];
    config.bam_region_scope = bam_region_scope;

    placer::apply_environment_config(config);
    if (final_fdr_q_set)
    {
        config.final_fdr_q = final_fdr_q;
    }
    if (final_report_mode_set)
    {
        config.final_report_mode = final_report_mode;
    }
    if (min_final_raw_cigar_insert_len_bp_set)
    {
        config.min_final_raw_cigar_insert_len_bp = min_final_raw_cigar_insert_len_bp;
    }
    return placer::run_pipeline_once(config);
}
