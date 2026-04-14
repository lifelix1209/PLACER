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
    write_fasta(ref_path, {{"chr1", std::string(5000, 'A')}});
    write_fasta(te_path, {{"UnknownTE#LTR/Unknown", std::string(400, 'A')}});

    std::vector<BamRecordPtr> records;
    const std::vector<uint32_t> cigar = {
        static_cast<uint32_t>(bam_cigar_gen(70, BAM_CMATCH)),
        static_cast<uint32_t>(bam_cigar_gen(80, BAM_CINS)),
        static_cast<uint32_t>(bam_cigar_gen(70, BAM_CMATCH)),
    };
    for (int i = 0; i < 12; ++i) {
        BamRecordPtr rec(bam_init1());
        const std::string qname = "read_" + std::to_string(i);
        const std::string seq(220, 'A');
        const std::string qual(220, 'I');
        assert(
            bam_set1(
                rec.get(),
                static_cast<size_t>(qname.size() + 1),
                qname.c_str(),
                0,
                0,
                1000 + (i * 40),
                60,
                static_cast<size_t>(cigar.size()),
                cigar.data(),
                -1,
                -1,
                0,
                static_cast<size_t>(seq.size()),
                seq.c_str(),
                qual.c_str(),
                0) >= 0);
        records.push_back(std::move(rec));
    }

    const auto clone_records = [&]() {
        std::vector<BamRecordPtr> out;
        out.reserve(records.size());
        for (const auto& record : records) {
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

    Pipeline streaming(streaming_config, std::make_unique<ReplayBamReader>(clone_records()));
    Pipeline parallel_one(parallel_one_config, std::make_unique<ReplayBamReader>(clone_records()));
    Pipeline parallel_two(parallel_two_config, std::make_unique<ReplayBamReader>(clone_records()));

    PipelineResult streaming_result = streaming.run();
    PipelineResult parallel_one_result = parallel_one.run();
    PipelineResult parallel_two_result = parallel_two.run();

    const auto assert_same = [](const PipelineResult& lhs, const PipelineResult& rhs) {
        assert(lhs.total_reads == rhs.total_reads);
        assert(lhs.gate1_passed == rhs.gate1_passed);
        assert(lhs.processed_bins == rhs.processed_bins);
        assert(lhs.built_components == rhs.built_components);
        assert(lhs.event_consensus_calls == rhs.event_consensus_calls);
        assert(lhs.genotype_calls == rhs.genotype_calls);
        assert(lhs.final_pass_calls == rhs.final_pass_calls);
        assert(lhs.final_calls.size() == rhs.final_calls.size());
        for (size_t i = 0; i < lhs.final_calls.size(); ++i) {
            assert(lhs.final_calls[i].chrom == rhs.final_calls[i].chrom);
            assert(lhs.final_calls[i].pos == rhs.final_calls[i].pos);
            assert(lhs.final_calls[i].final_qc == rhs.final_calls[i].final_qc);
        }
    };

    assert_same(streaming_result, parallel_one_result);
    assert_same(streaming_result, parallel_two_result);
    assert_same(parallel_one_result, parallel_two_result);

    std::remove(ref_path.c_str());
    std::remove((ref_path + ".fai").c_str());
    std::remove(te_path.c_str());
    std::remove((te_path + ".fai").c_str());
    return 0;
}
