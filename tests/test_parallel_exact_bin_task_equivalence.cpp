#include <cassert>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <htslib/sam.h>

#define private public
#include "pipeline.h"
#undef private

#include "parallel_executor_internal.h"

namespace {

using placer::BamRecordPtr;
using placer::BufferedRecord;

BamRecordPtr make_record(
    const std::string& qname,
    int32_t tid,
    int32_t pos0,
    const std::vector<uint32_t>& cigar,
    const std::string& seq) {
    BamRecordPtr record(bam_init1());
    assert(record != nullptr);
    const std::string qual(seq.size(), 'I');
    const int ret = bam_set1(
        record.get(),
        static_cast<size_t>(qname.size() + 1),
        qname.c_str(),
        0,
        tid,
        pos0,
        60,
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

BufferedRecord make_buffered(
    const std::string& qname,
    int32_t tid,
    int32_t pos0,
    int32_t ref_span) {
    const std::string seq(static_cast<size_t>(ref_span), 'A');
    const std::vector<uint32_t> cigar = {
        static_cast<uint32_t>(bam_cigar_gen(ref_span, BAM_CMATCH))};
    BufferedRecord rec;
    rec.record = make_record(qname, tid, pos0, cigar, seq);
    rec.ref_end = pos0 + ref_span;
    return rec;
}

}  // namespace

int main() {
    using namespace placer;

    std::vector<BufferedRecord> bin0;
    std::vector<BufferedRecord> bin1;
    std::vector<BufferedRecord> bin2;

    bin0.push_back(make_buffered("bin0_keep", 0, 3000, 2500));
    bin0.push_back(make_buffered("bin0_drop", 0, 0, 1000));
    bin1.push_back(make_buffered("bin1_keep_a", 0, 15000, 2000));
    bin1.push_back(make_buffered("bin1_keep_b", 0, 18000, 2000));
    bin2.push_back(make_buffered("bin2_keep_a", 0, 20500, 1000));
    bin2.push_back(make_buffered("bin2_keep_b", 0, 24000, 1000));

    std::vector<ExactBinSnapshot> snapshots;
    snapshots.push_back(ExactBinSnapshot{
        0, 0, std::make_shared<const std::vector<BufferedRecord>>(std::move(bin0))});
    snapshots.push_back(ExactBinSnapshot{
        0, 1, std::make_shared<const std::vector<BufferedRecord>>(std::move(bin1))});
    snapshots.push_back(ExactBinSnapshot{
        0, 2, std::make_shared<const std::vector<BufferedRecord>>(std::move(bin2))});

    ExactBinTask task = build_exact_bin_task(
        0,
        2,
        snapshots,
        10000,
        2,
        5000);

    assert(task.tid == 0);
    assert(task.bin_index == 2);
    assert(task.bin_start == 20000);
    assert(task.bin_end == 30000);
    assert(task.records.size() == 5);
    assert(task.records.front().record->core.pos == 3000);
    assert(task.records[1].record->core.pos == 15000);
    assert(task.records[2].record->core.pos == 18000);
    assert(task.records[3].record->core.pos == 20500);
    assert(task.records[4].record->core.pos == 24000);
    assert(task.current_bin_records == 2);
    assert(task.context_bin_records == 3);
    return 0;
}
