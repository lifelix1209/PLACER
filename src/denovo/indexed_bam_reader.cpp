#include "denovo.h"

#include <algorithm>
#include <iostream>
#include <string>

#include <htslib/hts.h>
#include <htslib/sam.h>

namespace placer {

IndexedBamReader::IndexedBamReader(std::string bam_path, int32_t decompression_threads)
    : bam_path_(std::move(bam_path)) {
    file_ = hts_open(bam_path_.c_str(), "r");
    if (file_ == nullptr) {
        std::cerr << "[DENOVO] failed to open BAM: " << bam_path_ << '\n';
        return;
    }

    if (decompression_threads > 1) {
        hts_set_threads(file_, decompression_threads);
    }

    header_ = sam_hdr_read(file_);
    if (header_ == nullptr) {
        std::cerr << "[DENOVO] failed to read BAM header: " << bam_path_ << '\n';
        return;
    }

    index_ = sam_index_load(file_, bam_path_.c_str());
    valid_ = true;
}

IndexedBamReader::~IndexedBamReader() {
    if (index_ != nullptr) {
        hts_idx_destroy(index_);
    }
    if (header_ != nullptr) {
        bam_hdr_destroy(header_);
    }
    if (file_ != nullptr) {
        hts_close(file_);
    }
}

int32_t IndexedBamReader::chromosome_count() const {
    return header_ ? header_->n_targets : 0;
}

std::string IndexedBamReader::chromosome_name(int32_t tid) const {
    if (header_ == nullptr || tid < 0 || tid >= header_->n_targets) {
        return "";
    }
    return header_->target_name[tid];
}

bool IndexedBamReader::fetch(
    const std::string& chrom,
    int32_t start,
    int32_t end,
    const std::function<bool(const bam1_t*)>& handler) const {
    if (!valid_ || index_ == nullptr || handler == nullptr) {
        return false;
    }

    const int32_t tid = sam_hdr_name2tid(header_, chrom.c_str());
    if (tid < 0) {
        return false;
    }

    const int32_t query_start = std::max<int32_t>(0, start);
    const int32_t query_end = std::max(query_start + 1, end);

    hts_itr_t* iter = sam_itr_queryi(index_, tid, query_start, query_end);
    if (iter == nullptr) {
        return false;
    }

    bam1_t* record = bam_init1();
    if (record == nullptr) {
        sam_itr_destroy(iter);
        return false;
    }

    int rc = 0;
    bool keep_iterating = true;
    while ((rc = sam_itr_next(file_, iter, record)) >= 0) {
        if (!handler(record)) {
            keep_iterating = false;
            break;
        }
    }

    bam_destroy1(record);
    sam_itr_destroy(iter);
    return keep_iterating && rc >= -1;
}

}  // namespace placer
