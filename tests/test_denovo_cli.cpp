#include "denovo.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

bool require_true(bool condition, const char* expr, const char* file, int line) {
    if (condition) {
        return true;
    }
    std::cerr << file << ":" << line << ": check failed: " << expr << std::endl;
    return false;
}

#define REQUIRE(expr) \
    do { \
        if (!require_true((expr), #expr, __FILE__, __LINE__)) { \
            return 1; \
        } \
    } while (0)

int write_parent_bam_with_pysam(
    const std::filesystem::path& work_dir,
    const std::filesystem::path& bam_path,
    const std::string& chrom,
    int32_t chrom_len,
    const std::string& qname,
    int32_t pos0,
    const std::string& seq,
    size_t softclip_len,
    size_t match_len) {
    const std::filesystem::path script_path = work_dir / "make_test_bam.py";
    {
        std::ofstream script(script_path);
        script << "import pysam\n"
               << "import sys\n"
               << "bam_path, chrom, chrom_len, qname, pos0, seq, softclip_len, match_len = sys.argv[1:9]\n"
               << "header = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': chrom, 'LN': int(chrom_len)}]}\n"
               << "with pysam.AlignmentFile(bam_path, 'wb', header=header) as out:\n"
               << "    seg = pysam.AlignedSegment()\n"
               << "    seg.query_name = qname\n"
               << "    seg.query_sequence = seq\n"
               << "    seg.flag = 0\n"
               << "    seg.reference_id = 0\n"
               << "    seg.reference_start = int(pos0)\n"
               << "    seg.mapping_quality = 60\n"
               << "    seg.cigartuples = [(4, int(softclip_len)), (0, int(match_len))]\n"
               << "    seg.query_qualities = pysam.qualitystring_to_array('I' * len(seq))\n"
               << "    out.write(seg)\n"
               << "pysam.index(bam_path)\n";
    }

    const std::string command =
        "python3 '" + script_path.string() + "' '" + bam_path.string() + "' '" + chrom + "' '" +
        std::to_string(chrom_len) + "' '" + qname + "' '" + std::to_string(pos0) + "' '" + seq +
        "' '" + std::to_string(softclip_len) + "' '" + std::to_string(match_len) + "'";
    return std::system(command.c_str());
}

}  // namespace

int main() {
    namespace fs = std::filesystem;

    const fs::path temp_dir = fs::temp_directory_path() / "placer_test_denovo_cli";
    fs::remove_all(temp_dir);
    fs::create_directories(temp_dir);

    const fs::path child_scientific = temp_dir / "child.scientific.txt";
    {
        std::ofstream out(child_scientific);
        out << "#PLACER streaming pipeline summary\n";
        out << "schema_version\t0.0.2\n";
        out << "#chrom\tpos\tte\tte_status\tconfidence\ttier\tsupport_reads\tte_bp_win_start\tte_bp_win_end\tte_ref_junc_min\tte_ref_junc_max\n";
        out << "chr1\t100\tLINE1HS\tTE_CANDIDATE\tHIGH\t1\t3\t80\t120\t95\t105\n";
        out << "chr2\t200\tALU\tTE_CANDIDATE\tHIGH\t2\t2\t195\t205\t-1\t-1\n";
        out << "chr3\t300\tNA\tNON_TE\tHIGH\t1\t5\t295\t305\t-1\t-1\n";
    }

    const fs::path parent_bams = temp_dir / "parents.fofn";
    {
        std::ofstream out(parent_bams);
        out << "/tmp/parent1.bam\n";
        out << "/tmp/parent2.bam\n";
    }

    std::vector<std::string> args = {
        "denovo",
        "--child-scientific", child_scientific.string(),
        "--parent-bam-list", parent_bams.string(),
        "--ref", "ref.fa",
        "--te", "te.fa",
        "--out-prefix", (temp_dir / "denovo_out").string()
    };

    std::vector<char*> argv;
    argv.reserve(args.size());
    for (auto& arg : args) {
        argv.push_back(arg.data());
    }

    placer::DenovoConfig config;
    std::string error_message;
    const bool ok = placer::parse_denovo_cli_args(
        static_cast<int>(argv.size()),
        argv.data(),
        config,
        error_message);
    REQUIRE(ok);
    REQUIRE(error_message.empty());
    REQUIRE(config.child_scientific_path == child_scientific.string());
    REQUIRE(config.parent_bam_list_path == parent_bams.string());
    REQUIRE(config.reference_fasta_path == "ref.fa");
    REQUIRE(config.te_fasta_path == "te.fa");

    int64_t total_rows = 0;
    const auto candidates = placer::load_denovo_child_candidates(config, &total_rows);
    REQUIRE(total_rows == 3);
    REQUIRE(candidates.size() == 2);
    REQUIRE(candidates[0].chrom == "chr1");
    REQUIRE(candidates[0].event_start == 95);
    REQUIRE(candidates[0].event_end == 105);
    REQUIRE(candidates[1].chrom == "chr2");
    REQUIRE(candidates[1].event_start == 195);
    REQUIRE(candidates[1].event_end == 205);

    const fs::path ref_fa = temp_dir / "ref.fa";
    {
        std::ofstream out(ref_fa);
        out << ">chr1\nACGTACGTACGT\n";
    }

    const fs::path te_fa = temp_dir / "te.fa";
    {
        std::ofstream out(te_fa);
        out << ">LINE1HS\nACGTACGTACGTACGTACGT\n";
        out << ">ALU\nTTTTCCCCAAAAGGGGTTTT\n";
    }

    config.reference_fasta_path = ref_fa.string();
    config.te_fasta_path = te_fa.string();
    config.dry_run = true;

    const placer::DenovoResult dry_run_result = placer::run_denovo(config);
    REQUIRE(dry_run_result.child_rows_total == 3);
    REQUIRE(dry_run_result.child_candidates_considered == 2);
    REQUIRE(dry_run_result.parent_bams == 2);
    REQUIRE(dry_run_result.calls_written == 2);
    REQUIRE(dry_run_result.implementation_status == "DRY_RUN");
    REQUIRE(dry_run_result.calls.size() == 2);
    REQUIRE(dry_run_result.calls[0].status == "DRY_RUN");
    REQUIRE(dry_run_result.calls[0].de_novo == "DRY_RUN");

    const std::string line1_seq =
        "ACGTTGCAAGTCCTGATCGATGCTAGCTTGACTGACCTGATCGTAGCTAGGCTAATCG";
    const std::string anchor_seq(50, 'C');
    const fs::path parent_bam = temp_dir / "parent1.bam";
    REQUIRE(write_parent_bam_with_pysam(
        temp_dir,
        parent_bam,
        "chr1",
        1000,
        "parent_read_1",
        100,
        line1_seq.substr(0, 50) + anchor_seq,
        50,
        50) == 0);
    {
        std::ofstream out(parent_bams);
        out << parent_bam.string() << "\n";
    }

    const fs::path ref_fa_real = temp_dir / "ref_real.fa";
    {
        std::ofstream out(ref_fa_real);
        out << ">chr1\n";
        out << std::string(1200, 'C') << "\n";
        out << ">chr2\n";
        out << std::string(1200, 'G') << "\n";
    }

    const fs::path te_fa_real = temp_dir / "te_real.fa";
    {
        std::ofstream out(te_fa_real);
        out << ">LINE1HS\n" << line1_seq << "\n";
        out << ">ALU\nTTGGAACCTTGGAACCTTGGAACCTTGGAACCTTGGAACC\n";
    }

    config.reference_fasta_path = ref_fa_real.string();
    config.te_fasta_path = te_fa_real.string();
    config.parent_bam_list_path = parent_bams.string();
    config.parent_bam_paths.clear();
    config.dry_run = false;

    const placer::DenovoResult scanned = placer::run_denovo(config);
    REQUIRE(scanned.implementation_status == "TARGETED_PARENT_SCAN_WITH_SPLIT_FRAGMENT_CLASSIFICATION");
    REQUIRE(scanned.parent_bams == 1);
    REQUIRE(scanned.calls.size() == 2);
    REQUIRE(scanned.parent_veto_calls == 1);
    REQUIRE(scanned.denovo_pass_calls == 1);
    REQUIRE(scanned.review_calls == 0);
    REQUIRE(scanned.calls[0].status == "PARENT_VETO");
    REQUIRE(scanned.calls[0].de_novo == "0");
    REQUIRE(scanned.calls[0].parent_summary.exact_te_reads == 1);
    REQUIRE(scanned.calls[0].parent_summary.total_support_reads == 1);
    REQUIRE(scanned.calls[1].status == "DENOVO_PASS");
    REQUIRE(scanned.calls[1].de_novo == "1");

    fs::remove_all(temp_dir);
    return 0;
}
