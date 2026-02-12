#ifndef PLACER_BAM_IO_H
#define PLACER_BAM_IO_H

#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <string_view>

#include <htslib/sam.h>

namespace placer {

struct BamRecordDeleter {
    void operator()(bam1_t* b) const {
        if (b != nullptr) {
            bam_destroy1(b);
        }
    }
};

using BamRecordPtr = std::unique_ptr<bam1_t, BamRecordDeleter>;

class ReadView {
public:
    explicit ReadView(const bam1_t* record) : record_(record) {}

    const bam1_t* raw() const { return record_; }

    int32_t tid() const;
    int32_t pos() const;
    int32_t mapq() const;
    uint16_t flag() const;

    int32_t seq_len() const;
    const uint8_t* seq_encoded() const;

    int32_t n_cigar() const;
    const uint32_t* cigar() const;

    std::string_view qname() const;

    bool has_tag(const char tag[2]) const;
    bool has_sa_tag() const;
    bool has_md_tag() const;
    bool get_int_tag(const char tag[2], int64_t& value) const;
    bool get_string_tag(const char tag[2], std::string& value) const;

    // Decode only [start, start+length) in read coordinates (0-based).
    std::string decode_subsequence(int32_t start, int32_t length) const;
    std::string decode_sequence() const;

private:
    const bam1_t* record_ = nullptr;
};

using RecordHandler = std::function<void(BamRecordPtr&&)>;
using ProgressHandler = std::function<bool(int64_t processed, int32_t current_tid)>;

class BamStreamReader {
public:
    virtual ~BamStreamReader() = default;

    virtual bool is_valid() const = 0;
    virtual const std::string& bam_path() const = 0;
    virtual int32_t chromosome_count() const = 0;
    virtual std::string chromosome_name(int32_t tid) const = 0;

    virtual int64_t stream(
        const RecordHandler& record_handler,
        const ProgressHandler& progress_handler = nullptr,
        int64_t progress_interval = 100000) = 0;
};

std::unique_ptr<BamStreamReader> make_bam_reader(
    const std::string& bam_path,
    int32_t decompression_threads = 2);

}  // namespace placer

#endif  // PLACER_BAM_IO_H
