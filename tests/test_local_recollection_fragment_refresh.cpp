#include <cassert>
#include <cstdio>
#include <fstream>
#include <limits>
#include <memory>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

#include <htslib/sam.h>

#define private public
#include "pipeline.h"
#undef private

namespace {

using BamPtr = placer::BamRecordPtr;

std::string make_temp_path(const std::string& stem, const std::string& ext) {
    const long long pid = static_cast<long long>(::getpid());
    return "/tmp/" + stem + "_" + std::to_string(pid) + ext;
}

std::string build_sequence(size_t len, uint32_t seed) {
    static const char kBases[] = {'A', 'C', 'G', 'T'};
    uint32_t state = seed;
    std::string out;
    out.reserve(len);
    for (size_t i = 0; i < len; ++i) {
        state = (state * 1664525u) + 1013904223u;
        out.push_back(kBases[(state >> 24) & 3u]);
    }
    return out;
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

BamPtr make_record(
    const std::string& qname,
    int32_t pos0,
    int32_t mapq,
    const std::vector<uint32_t>& cigar,
    const std::string& seq) {
    BamPtr record(bam_init1());
    assert(record != nullptr);

    const std::string qual(seq.size(), 'I');
    const int ret = bam_set1(
        record.get(),
        static_cast<size_t>(qname.size() + 1),
        qname.c_str(),
        0,
        0,
        pos0,
        mapq,
        static_cast<size_t>(cigar.size()),
        cigar.data(),
        -1,
        -1,
        0,
        static_cast<size_t>(seq.size()),
        seq.c_str(),
        qual.c_str(),
        0);
    assert(ret >= 0);
    return record;
}

BamPtr clone_record(const bam1_t* record) {
    return BamPtr(bam_dup1(record));
}

class FakeBamReader : public placer::BamStreamReader {
public:
    explicit FakeBamReader(std::vector<BamPtr> fetch_records)
        : fetch_records_(std::move(fetch_records)) {}

    bool is_valid() const override { return true; }

    std::unique_ptr<placer::BamStreamReader> clone(
        int32_t decompression_threads) const override {
        (void)decompression_threads;
        std::vector<BamPtr> copies;
        copies.reserve(fetch_records_.size());
        for (const auto& record : fetch_records_) {
            copies.push_back(clone_record(record.get()));
        }
        return std::make_unique<FakeBamReader>(std::move(copies));
    }

    const std::string& bam_path() const override { return bam_path_; }

    int32_t chromosome_count() const override { return 1; }

    std::string chromosome_name(int32_t tid) const override {
        return tid == 0 ? std::string("chr1") : std::string();
    }

    bool can_fetch() const override { return true; }

    bool fetch(
        const std::string& chrom,
        int32_t start,
        int32_t end,
        const placer::FetchRecordHandler& record_handler) const override {
        (void)start;
        (void)end;
        if (chrom != "chr1" || !record_handler) {
            return false;
        }
        for (const auto& record : fetch_records_) {
            if (!record_handler(clone_record(record.get()))) {
                return false;
            }
        }
        return true;
    }

    int64_t stream(
        const placer::RecordHandler& record_handler,
        const placer::ProgressHandler& progress_handler = nullptr,
        int64_t progress_interval = 100000) override {
        (void)record_handler;
        (void)progress_handler;
        (void)progress_interval;
        return 0;
    }

private:
    std::string bam_path_ = "fake.bam";
    std::vector<BamPtr> fetch_records_;
};

}  // namespace

int main() {
    using namespace placer;

    const std::string ref_path = make_temp_path("placer_local_refresh_ref", ".fa");
    const std::string te_path = make_temp_path("placer_local_refresh_te", ".fa");

    const std::string reference = build_sequence(900, 17u);
    const std::string te_insert = build_sequence(80, 101u);
    write_fasta(ref_path, {{"chr1", reference}});
    write_fasta(te_path, {{"SubGypsyRefresh#LTR/Gypsy", te_insert}});

    const int32_t ref_pos0 = 180;
    const std::string left = reference.substr(static_cast<size_t>(ref_pos0), 70);
    const std::string right = reference.substr(static_cast<size_t>(ref_pos0 + 70), 70);
    const std::string read_seq = left + te_insert + right;
    const std::vector<uint32_t> cigar = {
        bam_cigar_gen(70, BAM_CMATCH),
        bam_cigar_gen(80, BAM_CINS),
        bam_cigar_gen(70, BAM_CMATCH),
    };

    std::vector<BamPtr> owned_bin_records;
    owned_bin_records.push_back(make_record("alt_read_bin_only", ref_pos0, 60, cigar, read_seq));

    std::vector<const bam1_t*> bin_records;
    std::vector<BamPtr> fetch_records;
    for (const auto& record : owned_bin_records) {
        bin_records.push_back(record.get());
        fetch_records.push_back(clone_record(record.get()));
    }
    fetch_records.push_back(make_record("alt_read_local_only", ref_pos0, 60, cigar, read_seq));

    PipelineConfig config;
    config.reference_fasta_path = ref_path;
    config.te_fasta_path = te_path;
    config.te_kmer_size = 9;
    config.te_kmer_sizes_csv = "9,11";
    config.event_consensus_poa_min_reads = 2;

    Pipeline pipeline(
        config,
        std::make_unique<FakeBamReader>(std::move(fetch_records)));

    PipelineResult result;
    pipeline.process_bin_records(
        std::move(bin_records),
        0,
        0,
        std::numeric_limits<int32_t>::min(),
        std::numeric_limits<int32_t>::max(),
        result);
    finalize_final_calls(result);

    assert(result.event_consensus_calls == 1);
    assert(result.genotype_calls == 1);
    assert(result.final_pass_calls == 1);
    assert(result.final_calls.size() == 1);
    const FinalCall& call = result.final_calls.front();
    assert(call.chrom == "chr1");
    assert(call.family == "Gypsy");
    assert(call.subfamily == "SubGypsyRefresh");
    assert(call.insert_len == 80);
    assert(call.support_reads == 2);
    assert(call.alt_struct_reads == 2);
    assert(call.ref_span_reads == 0);
    assert(call.gq >= 20);
    assert(call.best_te_identity >= 0.80);
    assert(call.best_te_query_coverage >= 0.80);
    assert(call.event_consensus_len >= 200);
    assert(call.final_qc.rfind("PASS_FINAL_TE_CALL", 0) == 0);

    std::remove(ref_path.c_str());
    std::remove((ref_path + ".fai").c_str());
    std::remove(te_path.c_str());
    return 0;
}
