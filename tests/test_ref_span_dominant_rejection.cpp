#include <cassert>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <htslib/sam.h>

#define private public
#include "pipeline.h"
#undef private

namespace {

using BamPtr = placer::BamRecordPtr;

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

placer::ReadReferenceSpan make_span(int32_t tid, int32_t start, int32_t end) {
    placer::ReadReferenceSpan span;
    span.valid = true;
    span.tid = tid;
    span.start = start;
    span.end = end;
    return span;
}

}  // namespace

int main() {
    using namespace placer;

    const std::string reference = build_sequence(700, 17u);
    const std::string te_insert = build_sequence(80, 101u);
    const int32_t ref_pos0 = 150;
    const std::string left = reference.substr(static_cast<size_t>(ref_pos0), 70);
    const std::string right = reference.substr(static_cast<size_t>(ref_pos0 + 70), 70);
    const std::string alt_read_seq = left + te_insert + right;
    const std::vector<uint32_t> alt_cigar = {
        bam_cigar_gen(70, BAM_CMATCH),
        bam_cigar_gen(80, BAM_CINS),
        bam_cigar_gen(70, BAM_CMATCH),
    };

    const int32_t ref_read_pos0 = 120;
    const std::string ref_read_seq = reference.substr(static_cast<size_t>(ref_read_pos0), 160);
    const std::vector<uint32_t> ref_cigar = {
        bam_cigar_gen(160, BAM_CMATCH),
    };

    std::vector<BamPtr> owned_records;
    owned_records.push_back(make_record("alt_read_1", ref_pos0, 60, alt_cigar, alt_read_seq));
    owned_records.push_back(make_record("alt_read_2", ref_pos0, 60, alt_cigar, alt_read_seq));
    for (int i = 0; i < 10; ++i) {
        owned_records.push_back(make_record(
            "ref_read_" + std::to_string(i),
            ref_read_pos0,
            60,
            ref_cigar,
            ref_read_seq));
    }

    std::vector<const bam1_t*> bin_records;
    std::vector<BamPtr> fetch_records;
    std::vector<ReadReferenceSpan> read_spans;
    for (size_t i = 0; i < owned_records.size(); ++i) {
        bin_records.push_back(owned_records[i].get());
        fetch_records.push_back(clone_record(owned_records[i].get()));
        if (i < 2) {
            read_spans.push_back(make_span(0, ref_pos0, ref_pos0 + 140));
        } else {
            read_spans.push_back(make_span(0, ref_read_pos0, ref_read_pos0 + 160));
        }
    }

    PipelineConfig config;
    config.event_consensus_poa_min_reads = 2;
    Pipeline pipeline(
        config,
        std::make_unique<FakeBamReader>(std::move(fetch_records)));

    const auto components = pipeline.component_module_.build(
        bin_records,
        "chr1",
        0,
        0,
        config.bin_size);
    assert(components.size() == 1);

    const auto fragments = pipeline.ins_fragment_module_.extract(components.front(), bin_records);
    assert(fragments.size() == 2);

    const EventReadEvidence event_evidence = pipeline.collect_event_read_evidence(
        components.front(),
        bin_records,
        read_spans,
        fragments);
    assert(event_evidence.alt_struct_reads == 2);
    assert(event_evidence.ref_span_reads == 10);

    const EventConsensus event_consensus = pipeline.build_event_consensus(
        components.front(),
        bin_records,
        fragments,
        event_evidence);
    assert(event_consensus.qc_pass);

    const GenotypeCall genotype = pipeline.genotype_call(
        components.front(),
        event_evidence);
    assert(genotype.genotype == "0/0");
    assert(genotype.gq < 20);

    std::vector<const bam1_t*> process_records;
    process_records.reserve(owned_records.size());
    for (const auto& record : owned_records) {
        process_records.push_back(record.get());
    }

    PipelineResult result;
    pipeline.process_bin_records(std::move(process_records), 0, 0, result);
    finalize_final_calls(result);

    assert(result.event_consensus_calls == 1);
    assert(result.genotype_calls == 1);
    assert(result.final_pass_calls == 0);
    assert(result.final_calls.empty());

    return 0;
}
