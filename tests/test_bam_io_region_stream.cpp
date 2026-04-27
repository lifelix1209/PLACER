#include "bam_io.h"

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

int write_test_bam(const std::filesystem::path& work_dir, const std::filesystem::path& bam_path) {
    const std::filesystem::path script_path = work_dir / "make_region_test_bam.py";
    {
        std::ofstream script(script_path);
        script << "import pysam\n"
               << "import sys\n"
               << "bam_path = sys.argv[1]\n"
               << "header = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': ["
               << "{'SN': 'chr1', 'LN': 1000}, {'SN': 'chr2', 'LN': 1000}]}\n"
               << "records = [('chr1', 10, 'chr1_in'), ('chr1', 150, 'chr1_out'), ('chr2', 20, 'chr2_in')]\n"
               << "with pysam.AlignmentFile(bam_path, 'wb', header=header) as out:\n"
               << "    for chrom, pos, qname in records:\n"
               << "        seg = pysam.AlignedSegment()\n"
               << "        seg.query_name = qname\n"
               << "        seg.query_sequence = 'A' * 50\n"
               << "        seg.flag = 0\n"
               << "        seg.reference_id = out.get_tid(chrom)\n"
               << "        seg.reference_start = pos\n"
               << "        seg.mapping_quality = 60\n"
               << "        seg.cigartuples = [(0, 50)]\n"
               << "        seg.query_qualities = pysam.qualitystring_to_array('I' * 50)\n"
               << "        out.write(seg)\n"
               << "pysam.index(bam_path)\n";
    }

    const std::string command = "python3 '" + script_path.string() + "' '" + bam_path.string() + "'";
    return std::system(command.c_str());
}

}  // namespace

int main() {
    namespace fs = std::filesystem;

    const fs::path temp_dir = fs::temp_directory_path() / "placer_test_bam_io_region_stream";
    fs::remove_all(temp_dir);
    fs::create_directories(temp_dir);

    const fs::path bam_path = temp_dir / "input.bam";
    REQUIRE(write_test_bam(temp_dir, bam_path) == 0);

    placer::BamRegionScope scope;
    scope.enabled = true;
    scope.chrom = "chr1";
    scope.start = 0;
    scope.end = 120;

    auto reader = placer::make_bam_reader(bam_path.string(), 1, scope);
    REQUIRE(reader != nullptr);
    REQUIRE(reader->is_valid());
    REQUIRE(reader->can_fetch());

    std::vector<std::string> streamed;
    const int64_t total = reader->stream(
        [&](placer::BamRecordPtr&& record) {
            placer::ReadView view(record.get());
            streamed.emplace_back(view.qname());
        });

    REQUIRE(total == 1);
    REQUIRE(streamed.size() == 1);
    REQUIRE(streamed[0] == "chr1_in");

    std::vector<std::string> fetched;
    const bool fetch_ok = reader->fetch(
        "chr1",
        0,
        1000,
        [&](placer::BamRecordPtr&& record) {
            placer::ReadView view(record.get());
            fetched.emplace_back(view.qname());
            return true;
        });

    REQUIRE(fetch_ok);
    REQUIRE(fetched.size() == 1);
    REQUIRE(fetched[0] == "chr1_in");

    fs::remove_all(temp_dir);
    return 0;
}
