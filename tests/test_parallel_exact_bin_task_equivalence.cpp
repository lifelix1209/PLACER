#include <cstdio>
#include <cstdlib>
#include <deque>
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

void require(bool condition, const std::string& message) {
    if (!condition) {
        std::fprintf(
            stderr,
            "test_parallel_exact_bin_task_equivalence failed: %s\n",
            message.c_str());
        std::fflush(stderr);
        std::exit(1);
    }
}

BamRecordPtr make_record(
    const std::string& qname,
    int32_t tid,
    int32_t pos0,
    const std::vector<uint32_t>& cigar,
    const std::string& seq) {
    BamRecordPtr record(bam_init1());
    require(record != nullptr, "bam_init1 returned null");
    const std::string qual(seq.size(), 'I');
    const int ret = bam_set1(
        record.get(),
        static_cast<size_t>(qname.size()),
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
    require(ret >= 0, "bam_set1 failed");
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

    require(task.tid == 0, "unexpected task tid");
    require(task.bin_index == 2, "unexpected task bin index");
    require(task.bin_start == 20000, "unexpected task bin start");
    require(task.bin_end == 30000, "unexpected task bin end");
    require(task.records.size() == 4, "unexpected task record count");
    require(task.records.front().record->core.pos == 15000, "unexpected first record");
    require(task.records[1].record->core.pos == 18000, "unexpected second record");
    require(task.records[2].record->core.pos == 20500, "unexpected third record");
    require(task.records[3].record->core.pos == 24000, "unexpected fourth record");
    require(task.current_bin_records == 2, "unexpected current-bin record count");
    require(task.context_bin_records == 2, "unexpected context-bin record count");

    std::deque<ExactBinSnapshot> retained_snapshots;
    size_t next_owner_offset = 0;

    std::vector<BufferedRecord> owner_bin;
    std::vector<BufferedRecord> forward_bin_1;
    std::vector<BufferedRecord> forward_bin_2;
    owner_bin.push_back(make_buffered("owner", 0, 9880, 200));
    forward_bin_1.push_back(make_buffered("forward_1", 0, 10020, 200));
    forward_bin_2.push_back(make_buffered("forward_2", 0, 20020, 200));

    retained_snapshots.push_back(ExactBinSnapshot{
        0,
        0,
        std::make_shared<const std::vector<BufferedRecord>>(std::move(owner_bin)),
    });
    std::vector<ExactBinTask> materialized = materialize_ready_exact_bin_tasks(
        retained_snapshots,
        next_owner_offset,
        0,
        false,
        10000,
        2,
        5000);
    require(materialized.empty(), "owner bin materialized without forward context");

    retained_snapshots.push_back(ExactBinSnapshot{
        0,
        1,
        std::make_shared<const std::vector<BufferedRecord>>(std::move(forward_bin_1)),
    });
    materialized = materialize_ready_exact_bin_tasks(
        retained_snapshots,
        next_owner_offset,
        1,
        false,
        10000,
        2,
        5000);
    require(materialized.empty(), "owner bin materialized after only one forward bin");

    retained_snapshots.push_back(ExactBinSnapshot{
        0,
        2,
        std::make_shared<const std::vector<BufferedRecord>>(std::move(forward_bin_2)),
    });
    materialized = materialize_ready_exact_bin_tasks(
        retained_snapshots,
        next_owner_offset,
        2,
        false,
        10000,
        2,
        5000);
    require(materialized.size() == 1, "expected one ready task after symmetric context completed");
    require(materialized.front().bin_index == 0, "unexpected owner bin after delayed materialization");
    require(materialized.front().records.size() == 3, "symmetric task did not include forward records");
    require(materialized.front().records[0].record->core.pos == 9880, "unexpected owner record pos");
    require(materialized.front().records[1].record->core.pos == 10020, "missing first forward record");
    require(materialized.front().records[2].record->core.pos == 20020, "missing second forward record");
    require(
        materialized.front().current_bin_records == 1,
        "unexpected owner current-bin record count after delayed materialization");
    require(
        materialized.front().context_bin_records == 2,
        "unexpected owner context-bin record count after delayed materialization");
    return 0;
}
