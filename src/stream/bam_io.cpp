#include "bam_io.h"

#include <cstring>
#include <iostream>
#include <memory>
#include <string>
#include <utility>

#include <htslib/sam.h>

namespace placer {
namespace {

constexpr char kBamSeqChars[] = "=ACMGRSVTWYHKDBN";

class HtslibBamStreamReader final : public BamStreamReader {
public:
    explicit HtslibBamStreamReader(std::string bam_path, int32_t decompression_threads)
        : bam_path_(std::move(bam_path)) {
        file_ = hts_open(bam_path_.c_str(), "r");
        if (!file_) {
            std::cerr << "[BAM_IO] failed to open BAM: " << bam_path_ << '\n';
            return;
        }

        if (decompression_threads > 1) {
            hts_set_threads(file_, decompression_threads);
        }

        header_ = sam_hdr_read(file_);
        if (!header_) {
            std::cerr << "[BAM_IO] failed to read BAM header: " << bam_path_ << '\n';
            return;
        }

        valid_ = true;
    }

    ~HtslibBamStreamReader() override {
        if (header_ != nullptr) {
            bam_hdr_destroy(header_);
        }
        if (file_ != nullptr) {
            hts_close(file_);
        }
    }

    bool is_valid() const override { return valid_; }

    const std::string& bam_path() const override { return bam_path_; }

    int32_t chromosome_count() const override {
        return header_ ? header_->n_targets : 0;
    }

    std::string chromosome_name(int32_t tid) const override {
        if (!header_ || tid < 0 || tid >= header_->n_targets) {
            return "";
        }
        return header_->target_name[tid];
    }

    int64_t stream(
        const RecordHandler& record_handler,
        const ProgressHandler& progress_handler,
        int64_t progress_interval) override {
        if (!valid_ || !record_handler) {
            return -1;
        }

        int64_t processed = 0;
        int64_t last_progress = 0;
        int32_t current_tid = -1;

        while (true) {
            BamRecordPtr record(bam_init1());
            if (!record) {
                std::cerr << "[BAM_IO] failed to allocate BAM record\n";
                return -1;
            }

            const int rc = sam_read1(file_, header_, record.get());
            if (rc < 0) {
                break;
            }

            if (record->core.flag & BAM_FSECONDARY) {
                continue;
            }
            if (record->core.flag & BAM_FUNMAP) {
                continue;
            }

            current_tid = record->core.tid;
            record_handler(std::move(record));
            ++processed;

            if (progress_handler && progress_interval > 0 &&
                processed - last_progress >= progress_interval) {
                if (!progress_handler(processed, current_tid)) {
                    break;
                }
                last_progress = processed;
            }
        }

        if (progress_handler) {
            progress_handler(processed, current_tid);
        }

        return processed;
    }

private:
    std::string bam_path_;
    htsFile* file_ = nullptr;
    bam_hdr_t* header_ = nullptr;
    bool valid_ = false;
};

}  // namespace

int32_t ReadView::tid() const {
    return record_ ? record_->core.tid : -1;
}

int32_t ReadView::pos() const {
    return record_ ? record_->core.pos : -1;
}

int32_t ReadView::mapq() const {
    return record_ ? record_->core.qual : 0;
}

uint16_t ReadView::flag() const {
    return record_ ? record_->core.flag : 0;
}

int32_t ReadView::seq_len() const {
    return record_ ? record_->core.l_qseq : 0;
}

const uint8_t* ReadView::seq_encoded() const {
    return record_ ? bam_get_seq(record_) : nullptr;
}

int32_t ReadView::n_cigar() const {
    return record_ ? record_->core.n_cigar : 0;
}

const uint32_t* ReadView::cigar() const {
    return record_ ? bam_get_cigar(record_) : nullptr;
}

std::string_view ReadView::qname() const {
    if (!record_) {
        return std::string_view();
    }
    const char* qn = bam_get_qname(record_);
    return qn ? std::string_view(qn) : std::string_view();
}

bool ReadView::has_tag(const char tag[2]) const {
    if (!record_ || !tag) {
        return false;
    }
    return bam_aux_get(record_, tag) != nullptr;
}

bool ReadView::has_sa_tag() const {
    return has_tag("SA");
}

bool ReadView::has_md_tag() const {
    return has_tag("MD");
}

bool ReadView::get_int_tag(const char tag[2], int64_t& value) const {
    if (!record_ || !tag) {
        return false;
    }
    uint8_t* aux = bam_aux_get(record_, tag);
    if (!aux) {
        return false;
    }
    value = bam_aux2i(aux);
    return true;
}

bool ReadView::get_string_tag(const char tag[2], std::string& value) const {
    value.clear();
    if (!record_ || !tag) {
        return false;
    }
    uint8_t* aux = bam_aux_get(record_, tag);
    if (!aux) {
        return false;
    }
    const char* str = bam_aux2Z(aux);
    if (!str) {
        return false;
    }
    value.assign(str);
    return true;
}

std::string ReadView::decode_subsequence(int32_t start, int32_t length) const {
    std::string decoded;
    if (!record_) {
        return decoded;
    }

    const int32_t n = record_->core.l_qseq;
    if (n <= 0 || length <= 0) {
        return decoded;
    }

    if (start < 0) {
        start = 0;
    }
    if (start >= n) {
        return decoded;
    }

    int32_t end = start + length;
    if (end > n) {
        end = n;
    }
    const int32_t out_len = end - start;
    if (out_len <= 0) {
        return decoded;
    }

    decoded.resize(static_cast<size_t>(out_len));
    const uint8_t* packed = bam_get_seq(record_);
    for (int32_t i = 0; i < out_len; ++i) {
        decoded[static_cast<size_t>(i)] = kBamSeqChars[bam_seqi(packed, start + i)];
    }
    return decoded;
}

std::string ReadView::decode_sequence() const {
    std::string decoded;
    if (!record_) {
        return decoded;
    }

    const int32_t n = record_->core.l_qseq;
    if (n <= 0) {
        return decoded;
    }

    decoded.resize(static_cast<size_t>(n));
    const uint8_t* packed = bam_get_seq(record_);
    for (int32_t i = 0; i < n; ++i) {
        decoded[static_cast<size_t>(i)] = kBamSeqChars[bam_seqi(packed, i)];
    }
    return decoded;
}

std::unique_ptr<BamStreamReader> make_bam_reader(
    const std::string& bam_path,
    int32_t decompression_threads) {
    return std::make_unique<HtslibBamStreamReader>(bam_path, decompression_threads);
}

}  // namespace placer
