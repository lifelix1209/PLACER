#include <algorithm>
#include <cassert>
#include <memory>
#include <string>
#include <vector>

#include <htslib/sam.h>

#include "pipeline.h"
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
    assert(ret >= 0);
    return record;
}

BufferedRecord make_insert_read(
    const std::string& qname,
    int32_t pos0,
    int32_t match_before,
    int32_t ins_len,
    int32_t match_after) {
    const std::vector<uint32_t> cigar = {
        static_cast<uint32_t>(bam_cigar_gen(match_before, BAM_CMATCH)),
        static_cast<uint32_t>(bam_cigar_gen(ins_len, BAM_CINS)),
        static_cast<uint32_t>(bam_cigar_gen(match_after, BAM_CMATCH)),
    };
    BufferedRecord out;
    out.record = make_record(
        qname,
        0,
        pos0,
        cigar,
        std::string(static_cast<size_t>(match_before + ins_len + match_after), 'A'));
    out.ref_end = pos0 + match_before + match_after;
    return out;
}

std::vector<const bam1_t*> task_record_ptrs(const placer::ExactBinTask& task) {
    std::vector<const bam1_t*> out;
    for (const auto& ref : task.records) {
        if (ref.record != nullptr) {
            out.push_back(ref.record);
        }
    }
    return out;
}

std::vector<placer::ComponentCall> owned_components(
    const placer::ExactBinTask& task) {
    const auto records = task_record_ptrs(task);
    placer::LinearBinComponentModule module;
    auto components = module.build(records, "chr1", 0, task.bin_start, task.bin_end);
    components.erase(
        std::remove_if(
            components.begin(),
            components.end(),
            [&](const placer::ComponentCall& component) {
                return component.anchor_pos < task.bin_start ||
                    component.anchor_pos >= task.bin_end;
            }),
        components.end());
    return components;
}

}  // namespace

int main() {
    std::vector<BufferedRecord> bin0;
    std::vector<BufferedRecord> bin1;
    bin0.push_back(make_insert_read("context_read", 9900, 120, 100, 100));
    bin1.push_back(make_insert_read("owner_read", 10000, 20, 100, 100));

    std::vector<placer::ExactBinSnapshot> snapshots;
    snapshots.push_back(placer::ExactBinSnapshot{
        0,
        0,
        std::make_shared<const std::vector<BufferedRecord>>(std::move(bin0)),
    });
    snapshots.push_back(placer::ExactBinSnapshot{
        0,
        1,
        std::make_shared<const std::vector<BufferedRecord>>(std::move(bin1)),
    });

    const placer::ExactBinTask owner0 = placer::build_exact_bin_task(
        0,
        0,
        snapshots,
        10000,
        2,
        5000);
    const auto owned0 = owned_components(owner0);
    assert(owned0.empty());

    const placer::ExactBinTask owner1 = placer::build_exact_bin_task(
        0,
        1,
        snapshots,
        10000,
        2,
        5000);
    const auto owned1 = owned_components(owner1);
    assert(owned1.size() == 1);
    assert(owned1.front().anchor_pos == 10020);
    assert(owned1.front().read_indices.size() == 2);
    assert(owned1.front().insertion_read_indices.size() == 2);
    return 0;
}
