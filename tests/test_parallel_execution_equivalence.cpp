#include <cassert>
#include <cstdio>
#include <fstream>
#include <memory>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

#define private public
#include "pipeline.h"
#undef private

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) {
        std::fprintf(stderr, "test_parallel_execution_equivalence failed: %s\n", message.c_str());
        std::fflush(stderr);
        std::exit(1);
    }
}

std::string make_temp_path(const std::string& stem, const std::string& ext) {
    const long long pid = static_cast<long long>(::getpid());
    return "/tmp/" + stem + "_" + std::to_string(pid) + ext;
}

void write_fasta(
    const std::string& path,
    const std::vector<std::pair<std::string, std::string>>& entries) {
    std::ofstream out(path);
    assert(out.is_open());
    for (const auto& entry : entries) {
        out << ">" << entry.first << "\n";
        out << entry.second << "\n";
    }
}

class ReplayBamReader : public placer::BamStreamReader {
public:
    explicit ReplayBamReader(std::vector<placer::BamRecordPtr> records)
        : records_(std::move(records)) {}

    bool is_valid() const override { return true; }
    std::unique_ptr<placer::BamStreamReader> clone(
        int32_t decompression_threads) const override {
        (void)decompression_threads;
        std::vector<placer::BamRecordPtr> copies;
        copies.reserve(records_.size());
        for (const auto& record : records_) {
            copies.push_back(placer::BamRecordPtr(bam_dup1(record.get())));
        }
        return std::make_unique<ReplayBamReader>(std::move(copies));
    }
    const std::string& bam_path() const override { return bam_path_; }
    int32_t chromosome_count() const override { return 1; }
    std::string chromosome_name(int32_t tid) const override { return tid == 0 ? "chr1" : ""; }
    bool can_fetch() const override { return true; }

    bool fetch(
        const std::string& chrom,
        int32_t start,
        int32_t end,
        const placer::FetchRecordHandler& handler) const override {
        if (chrom != "chr1") {
            return false;
        }
        for (const auto& record : records_) {
            placer::ReadView view(record.get());
            const int32_t ref_end = view.pos() + 140;
            if (ref_end <= start || view.pos() >= end) {
                continue;
            }
            if (!handler(placer::BamRecordPtr(bam_dup1(record.get())))) {
                return false;
            }
        }
        return true;
    }

    int64_t stream(
        const placer::RecordHandler& handler,
        const placer::ProgressHandler& progress_handler,
        int64_t progress_interval) override {
        int64_t count = 0;
        for (const auto& record : records_) {
            handler(placer::BamRecordPtr(bam_dup1(record.get())));
            ++count;
            if (progress_handler && progress_interval > 0 && (count % progress_interval) == 0) {
                progress_handler(count, 0);
            }
        }
        if (progress_handler) {
            progress_handler(count, 0);
        }
        return count;
    }

private:
    std::string bam_path_ = "replay.bam";
    std::vector<placer::BamRecordPtr> records_;
};

}  // namespace

int main() {
    using namespace placer;

    const std::string ref_path = make_temp_path("placer_parallel_exact_ref", ".fa");
    const std::string te_path = make_temp_path("placer_parallel_exact_te", ".fa");
    write_fasta(ref_path, {{"chr1", std::string(50000, 'A')}});
    write_fasta(te_path, {{"UnknownTE#LTR/Unknown", std::string(400, 'A')}});

    const std::vector<uint32_t> cigar = {
        static_cast<uint32_t>(bam_cigar_gen(70, BAM_CMATCH)),
        static_cast<uint32_t>(bam_cigar_gen(80, BAM_CINS)),
        static_cast<uint32_t>(bam_cigar_gen(70, BAM_CMATCH)),
    };
    const std::vector<uint32_t> background_cigar = {
        static_cast<uint32_t>(bam_cigar_gen(220, BAM_CMATCH)),
    };

    const auto make_records = [&](
                                  const std::vector<int32_t>& positions,
                                  const std::string& prefix,
                                  const std::vector<uint32_t>& record_cigar) {
        std::vector<BamRecordPtr> out;
        out.reserve(positions.size());
        const size_t seq_len = static_cast<size_t>(
            bam_cigar2qlen(
                static_cast<int>(record_cigar.size()),
                record_cigar.data()));
        for (size_t i = 0; i < positions.size(); ++i) {
            BamRecordPtr rec(bam_init1());
            const std::string qname = prefix + "_" + std::to_string(i);
            const std::string seq(seq_len, 'A');
            const std::string qual(seq_len, 'I');
            require(
                bam_set1(
                    rec.get(),
                    static_cast<size_t>(qname.size()),
                    qname.c_str(),
                    0,
                    0,
                    positions[i],
                    60,
                    static_cast<size_t>(record_cigar.size()),
                    record_cigar.data(),
                    -1,
                    -1,
                    0,
                    seq_len,
                    seq.c_str(),
                    qual.c_str(),
                    0) >= 0,
                "bam_set1 failed while building replay records");
            out.push_back(std::move(rec));
        }
        return out;
    };

    const std::vector<BamRecordPtr> records = make_records(
        {1000, 1040, 1080, 1120, 1160, 1200, 1240, 1280, 1320, 1360, 1400, 1440},
        "read",
        cigar);
    const std::vector<BamRecordPtr> spread_records = make_records(
        {1000, 13000, 25000, 37000},
        "spread",
        background_cigar);
    const std::vector<uint32_t> right_clip_cigar = {
        static_cast<uint32_t>(bam_cigar_gen(140, BAM_CMATCH)),
        static_cast<uint32_t>(bam_cigar_gen(120, BAM_CSOFT_CLIP)),
    };
    const std::vector<uint32_t> left_clip_cigar = {
        static_cast<uint32_t>(bam_cigar_gen(120, BAM_CSOFT_CLIP)),
        static_cast<uint32_t>(bam_cigar_gen(140, BAM_CMATCH)),
    };
    std::vector<BamRecordPtr> boundary_records = make_records(
        {9860, 9870},
        "boundary_right",
        right_clip_cigar);
    std::vector<BamRecordPtr> boundary_left_records = make_records(
        {10000, 10010},
        "boundary_left",
        left_clip_cigar);
    boundary_records.insert(
        boundary_records.end(),
        std::make_move_iterator(boundary_left_records.begin()),
        std::make_move_iterator(boundary_left_records.end()));

    const auto clone_records = [](const std::vector<BamRecordPtr>& src) {
        std::vector<BamRecordPtr> out;
        out.reserve(src.size());
        for (const auto& record : src) {
            out.push_back(BamRecordPtr(bam_dup1(record.get())));
        }
        return out;
    };

    PipelineConfig streaming_config;
    PipelineConfig parallel_one_config;
    PipelineConfig parallel_two_config;
    streaming_config.reference_fasta_path = ref_path;
    streaming_config.te_fasta_path = te_path;
    streaming_config.te_kmer_size = 9;
    streaming_config.te_kmer_sizes_csv = "9,11";
    parallel_one_config = streaming_config;
    parallel_one_config.enable_parallel = true;
    parallel_one_config.parallel_workers = 1;
    parallel_one_config.parallel_queue_max_tasks = 2;
    parallel_one_config.parallel_result_buffer_max = 2;
    parallel_two_config = streaming_config;
    parallel_two_config.enable_parallel = true;
    parallel_two_config.parallel_workers = 2;
    parallel_two_config.parallel_queue_max_tasks = 4;
    parallel_two_config.parallel_result_buffer_max = 4;

    Pipeline streaming(
        streaming_config,
        std::make_unique<ReplayBamReader>(clone_records(records)));
    Pipeline parallel_one(
        parallel_one_config,
        std::make_unique<ReplayBamReader>(clone_records(records)));
    Pipeline parallel_two(
        parallel_two_config,
        std::make_unique<ReplayBamReader>(clone_records(records)));

    PipelineResult streaming_result = streaming.run();
    PipelineResult parallel_one_result = parallel_one.run();
    PipelineResult parallel_two_result = parallel_two.run();

    const auto assert_same = [](const PipelineResult& lhs, const PipelineResult& rhs) {
        require(lhs.total_reads == rhs.total_reads, "total_reads mismatch");
        require(lhs.gate1_passed == rhs.gate1_passed, "gate1_passed mismatch");
        require(lhs.processed_bins == rhs.processed_bins, "processed_bins mismatch");
        require(lhs.built_components == rhs.built_components, "built_components mismatch");
        require(lhs.event_consensus_calls == rhs.event_consensus_calls, "event_consensus_calls mismatch");
        require(lhs.genotype_calls == rhs.genotype_calls, "genotype_calls mismatch");
        require(lhs.final_pass_calls == rhs.final_pass_calls, "final_pass_calls mismatch");
        require(lhs.final_calls.size() == rhs.final_calls.size(), "final_calls size mismatch");
        for (size_t i = 0; i < lhs.final_calls.size(); ++i) {
            require(lhs.final_calls[i].chrom == rhs.final_calls[i].chrom, "final_calls chrom mismatch");
            require(lhs.final_calls[i].pos == rhs.final_calls[i].pos, "final_calls pos mismatch");
            require(
                lhs.final_calls[i].final_qc == rhs.final_calls[i].final_qc,
                "final_calls final_qc mismatch");
        }
    };

    assert_same(streaming_result, parallel_one_result);
    assert_same(streaming_result, parallel_two_result);
    assert_same(parallel_one_result, parallel_two_result);

    Pipeline spread_parallel_one(
        parallel_one_config,
        std::make_unique<ReplayBamReader>(clone_records(spread_records)));
    Pipeline spread_parallel_two(
        parallel_two_config,
        std::make_unique<ReplayBamReader>(clone_records(spread_records)));

    PipelineResult spread_parallel_one_result = spread_parallel_one.run();
    PipelineResult spread_parallel_two_result = spread_parallel_two.run();

    require(spread_parallel_one_result.processed_bins >= 2, "parallel path did not split spread input into multiple bins");
    require(
        spread_parallel_one_result.total_reads == spread_parallel_two_result.total_reads,
        "spread total_reads mismatch between worker counts");
    require(
        spread_parallel_one_result.gate1_passed == spread_parallel_two_result.gate1_passed,
        "spread gate1_passed mismatch between worker counts");
    require(
        spread_parallel_one_result.processed_bins == spread_parallel_two_result.processed_bins,
        "spread processed_bins mismatch between worker counts");
    require(
        spread_parallel_one_result.built_components == spread_parallel_two_result.built_components,
        "spread built_components mismatch between worker counts");
    require(
        spread_parallel_one_result.event_consensus_calls ==
            spread_parallel_two_result.event_consensus_calls,
        "spread event_consensus_calls mismatch between worker counts");
    require(
        spread_parallel_one_result.final_pass_calls ==
            spread_parallel_two_result.final_pass_calls,
        "spread final_pass_calls mismatch between worker counts");
    require(
        spread_parallel_one_result.final_calls.size() ==
            spread_parallel_two_result.final_calls.size(),
        "spread final_calls size mismatch between worker counts");

    Pipeline boundary_streaming(
        streaming_config,
        std::make_unique<ReplayBamReader>(clone_records(boundary_records)));
    Pipeline boundary_parallel_one(
        parallel_one_config,
        std::make_unique<ReplayBamReader>(clone_records(boundary_records)));
    Pipeline boundary_parallel_two(
        parallel_two_config,
        std::make_unique<ReplayBamReader>(clone_records(boundary_records)));

    PipelineResult boundary_streaming_result = boundary_streaming.run();
    PipelineResult boundary_parallel_one_result = boundary_parallel_one.run();
    PipelineResult boundary_parallel_two_result = boundary_parallel_two.run();

    assert_same(boundary_streaming_result, boundary_parallel_one_result);
    assert_same(boundary_streaming_result, boundary_parallel_two_result);
    assert_same(boundary_parallel_one_result, boundary_parallel_two_result);

    std::remove(ref_path.c_str());
    std::remove((ref_path + ".fai").c_str());
    std::remove(te_path.c_str());
    std::remove((te_path + ".fai").c_str());
    return 0;
}
