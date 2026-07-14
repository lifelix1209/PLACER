constexpr int32_t kSoftClipSignalMin = 20;
constexpr int32_t kLongInsertionSignalMin = 50;
constexpr int32_t kLargeIndelEvidenceMin = 40;
constexpr double kSaHintWeight = 0.35;
constexpr int32_t kHistBinSize = 20;
constexpr int32_t kHistPadding = 120;
constexpr double kPeakMinWeight = 1.2;
constexpr int32_t kPeakMergeDistance = 100;
constexpr int32_t kWindowMergeGap = 60;
// Allow component windows to cross bin boundaries so insertion evidence carried
// by long reads (whose starts may fall in adjacent bins) is not dropped.
constexpr int32_t kWindowBinSlackBp = 5000;
// Merge current bin with a small number of previous bins so split support from
// long reads can be assembled in one component.
constexpr int32_t kCrossBinContextBins = 2;
constexpr int32_t kWindowMaxExpandBins = 15;
constexpr double kWindowDropRatio = 0.50;
constexpr int32_t kMinWindowSpan = 120;
constexpr double kReadAssignMinScore = 0.8;
constexpr double kReadAmbiguousRatio = 0.75;
constexpr double kReadAssignDecayBp = 150.0;
constexpr int32_t kLocalEventFetchSlackBp = 1000;
constexpr int32_t kEventConsensusFlankBp = 80;
constexpr int32_t kEventConsensusMinAnchorBp = 50;
constexpr int32_t kEventSegmentationMinFlankAlignBp = 50;
constexpr double kEventSegmentationMinFlankIdentity = 0.90;
constexpr int32_t kEventSegmentationBreakpointSlackBp = 200;
constexpr int32_t kEventSegmentationEndpointSlackBp = 32;
constexpr int32_t kEventSegmentationMinInsertBp = 1;
constexpr int32_t kEventSegmentationMaxFlankQueryBp = 120;
constexpr double kEventSegmentationUniquenessMargin = 0.02;
constexpr int32_t kEventSegmentationSeedK = 11;
constexpr int32_t kEventSegmentationSeedBinBp = 8;
constexpr int32_t kEventSegmentationSeedTopBins = 8;

enum class EvidenceKind : uint8_t {
    kSoftClip = 0,
    kIndel = 1,
    kSAHint = 2
};

struct EvidencePoint {
    size_t read_index = 0;
    int32_t pos = -1;
    double weight = 0.0;
    EvidenceKind kind = EvidenceKind::kSoftClip;
    int32_t signal_len = 0;
    uint8_t class_mask = 0;
};

struct ReadSignalSummary {
    std::string read_id;
    bool is_reverse = false;
    bool has_sa_or_supp = false;
    int32_t max_soft_clip = 0;
    int32_t max_ins = 0;
    uint8_t class_mask = 0;
};

struct EvidenceBundle {
    std::vector<EvidencePoint> points;
    std::vector<ReadSignalSummary> read_summaries;
};

struct ReadWindowAssignmentScore {
    bool valid = false;
    int32_t specificity_rank = std::numeric_limits<int32_t>::max();
    double specificity_score = 0.0;
    double total_score = 0.0;
};

struct LocalFetchedReads {
    std::vector<BamRecordPtr> owned_records;
    std::vector<const bam1_t*> records;
    std::vector<ReadReferenceSpan> read_spans;
};

std::vector<int32_t> robust_local_split_insertion_positions(
    const ReadView& read,
    const std::string& chrom);

ComponentCall build_local_fragment_component(
    const PipelineConfig& config,
    const ComponentCall& component,
    const std::vector<const bam1_t*>& local_records) {
    ComponentCall local_component = component;
    local_component.read_indices.clear();
    local_component.soft_clip_read_indices.clear();
    local_component.split_sa_read_indices.clear();
    local_component.insertion_read_indices.clear();

    local_component.read_indices.reserve(local_records.size());
    local_component.soft_clip_read_indices.reserve(local_records.size());
    local_component.split_sa_read_indices.reserve(local_records.size());
    local_component.insertion_read_indices.reserve(local_records.size());

    for (size_t idx = 0; idx < local_records.size(); ++idx) {
        const bam1_t* record = local_records[idx];
        if (!record) {
            continue;
        }

        ReadView read(record);
        if (read.tid() != component.tid) {
            continue;
        }

        local_component.read_indices.push_back(idx);

        const uint16_t flag = read.flag();
        if (!robust_local_split_insertion_positions(read, component.chrom).empty() ||
            ((flag & BAM_FSUPPLEMENTARY) != 0)) {
            local_component.split_sa_read_indices.push_back(idx);
        }

        const uint32_t* cigar = read.cigar();
        const int32_t n_cigar = read.n_cigar();
        if (!cigar || n_cigar <= 0) {
            continue;
        }

        int32_t max_soft_clip = 0;
        int32_t max_insertion = 0;
        for (int32_t ci = 0; ci < n_cigar; ++ci) {
            const int op = bam_cigar_op(cigar[ci]);
            const int32_t len = static_cast<int32_t>(bam_cigar_oplen(cigar[ci]));
            if (op == BAM_CSOFT_CLIP) {
                max_soft_clip = std::max(max_soft_clip, len);
            } else if (op == BAM_CINS) {
                max_insertion = std::max(max_insertion, len);
            }
        }

        if (max_soft_clip >= config.min_soft_clip_for_seq_extract) {
            local_component.soft_clip_read_indices.push_back(idx);
        }
        if (read.mapq() == kInsertionCandidateRequiredMapq &&
            max_insertion >= config.min_long_ins_for_seq_extract) {
            local_component.insertion_read_indices.push_back(idx);
        }
    }

    return local_component;
}

struct DensityPeak {
    int32_t pos = -1;
    double weight = 0.0;
};

struct CandidateWindow {
    int32_t start = -1;
    int32_t end = -1;
    int32_t center = -1;
    double peak_weight = 0.0;
};

struct EventFlankPlacement {
    int32_t query_start = -1;
    int32_t query_end = -1;
    int32_t ref_start = -1;
    int32_t ref_end = -1;
    int32_t align_len = 0;
    int32_t breakpoint_delta = std::numeric_limits<int32_t>::max();
    int32_t endpoint_offset = std::numeric_limits<int32_t>::max();
    double identity = 0.0;
};

struct EventFlankSearchResult {
    std::vector<EventFlankPlacement> candidates;
};

enum class SaHintSide : uint8_t {
    kNone = 0,
    kLeft = 1,
    kRight = 2
};

struct LocalCigarOp {
    int32_t len = 0;
    char op = '\0';
};

struct LocalSAEntry {
    std::string rname;
    int32_t pos = 0;  // 1-based in SA tag
    char strand = '+';
    std::string cigar;
};

struct LocalQueryInterval {
    int32_t qstart = -1;
    int32_t qend = -1;
    int32_t leading_s = 0;
    int32_t trailing_s = 0;
};

struct LocalNormalizedAlignment {
    std::string chrom;
    int32_t ref_start = -1;
    int32_t ref_end = -1;
    int32_t qstart = -1;
    int32_t qend = -1;
    bool is_reverse = false;
};

SaHintSide choose_sa_hint_side(int32_t leading_soft, int32_t trailing_soft) {
    const bool left = leading_soft >= kSoftClipSignalMin;
    const bool right = trailing_soft >= kSoftClipSignalMin;
    if (left && !right) {
        return SaHintSide::kLeft;
    }
    if (right && !left) {
        return SaHintSide::kRight;
    }
    if (left && right) {
        if (leading_soft > trailing_soft) {
            return SaHintSide::kLeft;
        }
        if (trailing_soft > leading_soft) {
            return SaHintSide::kRight;
        }
    }
    return SaHintSide::kNone;
}

bool local_consumes_query(char op) {
    return op == 'M' || op == '=' || op == 'X' || op == 'I' || op == 'S';
}

bool local_consumes_ref(char op) {
    return op == 'M' || op == '=' || op == 'X' || op == 'D' || op == 'N';
}

bool parse_local_int32(const char* text, int32_t& out) {
    if (!text || !*text) {
        return false;
    }
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (end == text || (end && *end != '\0')) {
        return false;
    }
    out = static_cast<int32_t>(value);
    return true;
}

std::vector<LocalSAEntry> parse_local_sa_tag_z(const std::string& sa_z) {
    std::vector<LocalSAEntry> out;
    size_t start = 0;
    while (start < sa_z.size()) {
        const size_t rec_end = sa_z.find(';', start);
        const size_t rec_len =
            (rec_end == std::string::npos) ? (sa_z.size() - start) : (rec_end - start);
        if (rec_len == 0) {
            break;
        }

        const std::string rec = sa_z.substr(start, rec_len);
        std::vector<std::string> fields;
        fields.reserve(6);
        size_t field_start = 0;
        while (field_start <= rec.size()) {
            const size_t comma = rec.find(',', field_start);
            if (comma == std::string::npos) {
                fields.push_back(rec.substr(field_start));
                break;
            }
            fields.push_back(rec.substr(field_start, comma - field_start));
            field_start = comma + 1;
        }

        if (fields.size() >= 4) {
            LocalSAEntry entry;
            entry.rname = fields[0];
            parse_local_int32(fields[1].c_str(), entry.pos);
            entry.strand = fields[2].empty() ? '+' : fields[2][0];
            entry.cigar = fields[3];
            if (!entry.rname.empty() && entry.pos > 0 && !entry.cigar.empty()) {
                out.push_back(std::move(entry));
            }
        }

        if (rec_end == std::string::npos) {
            break;
        }
        start = rec_end + 1;
    }
    return out;
}

bool parse_local_cigar_ops(const std::string& cigar, std::vector<LocalCigarOp>& ops) {
    ops.clear();
    int32_t len = 0;
    bool have_len = false;
    for (char c : cigar) {
        if (std::isdigit(static_cast<unsigned char>(c))) {
            have_len = true;
            len = (len * 10) + (c - '0');
            continue;
        }
        if (!have_len || len <= 0) {
            return false;
        }
        ops.push_back(LocalCigarOp{len, c});
        len = 0;
        have_len = false;
    }
    return !have_len && !ops.empty();
}

bool local_cigar_to_query_interval(
    const std::vector<LocalCigarOp>& ops,
    int32_t read_len,
    LocalQueryInterval& out,
    int32_t& ref_aligned_len) {
    out = LocalQueryInterval{};
    ref_aligned_len = 0;
    if (read_len <= 0 || ops.empty()) {
        return false;
    }

    int32_t first = 0;
    while (first < static_cast<int32_t>(ops.size()) &&
           ops[static_cast<size_t>(first)].op == 'H') {
        ++first;
    }
    int32_t last = static_cast<int32_t>(ops.size()) - 1;
    while (last >= 0 && ops[static_cast<size_t>(last)].op == 'H') {
        --last;
    }
    if (first > last) {
        return false;
    }
    if (ops[static_cast<size_t>(first)].op == 'S') {
        out.leading_s = ops[static_cast<size_t>(first)].len;
    }
    if (ops[static_cast<size_t>(last)].op == 'S') {
        out.trailing_s = ops[static_cast<size_t>(last)].len;
    }

    int32_t q_aligned = 0;
    for (const auto& op : ops) {
        if (local_consumes_query(op.op) && op.op != 'S') {
            q_aligned += op.len;
        }
        if (local_consumes_ref(op.op)) {
            ref_aligned_len += op.len;
        }
    }
    if (q_aligned <= 0 || ref_aligned_len <= 0) {
        return false;
    }

    int32_t qstart = 0;
    if (out.leading_s > 0) {
        qstart = out.leading_s;
    } else if (out.trailing_s > 0) {
        qstart = read_len - out.trailing_s - q_aligned;
    } else if (q_aligned == read_len) {
        qstart = 0;
    } else {
        return false;
    }
    if (qstart < 0) {
        return false;
    }
    out.qstart = qstart;
    out.qend = qstart + q_aligned;
    return out.qend > out.qstart && out.qend <= read_len;
}

bool local_bam_to_query_interval(
    const ReadView& read,
    LocalQueryInterval& out) {
    out = LocalQueryInterval{};
    const int32_t read_len = read.seq_len();
    const uint32_t* cigar = read.cigar();
    const int32_t n_cigar = read.n_cigar();
    if (read_len <= 0 || !cigar || n_cigar <= 0) {
        return false;
    }

    int first = -1;
    int last = -1;
    for (int32_t i = 0; i < n_cigar; ++i) {
        if (bam_cigar_op(cigar[i]) != BAM_CHARD_CLIP) {
            first = i;
            break;
        }
    }
    for (int32_t i = n_cigar - 1; i >= 0; --i) {
        if (bam_cigar_op(cigar[i]) != BAM_CHARD_CLIP) {
            last = i;
            break;
        }
    }
    if (first < 0 || last < 0) {
        return false;
    }
    if (bam_cigar_op(cigar[first]) == BAM_CSOFT_CLIP) {
        out.leading_s = static_cast<int32_t>(bam_cigar_oplen(cigar[first]));
    }
    if (bam_cigar_op(cigar[last]) == BAM_CSOFT_CLIP) {
        out.trailing_s = static_cast<int32_t>(bam_cigar_oplen(cigar[last]));
    }

    int32_t q_aligned = 0;
    for (int32_t i = 0; i < n_cigar; ++i) {
        const int op = bam_cigar_op(cigar[i]);
        const int32_t len = static_cast<int32_t>(bam_cigar_oplen(cigar[i]));
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF || op == BAM_CINS) {
            q_aligned += len;
        }
    }
    if (q_aligned <= 0) {
        return false;
    }

    out.qstart = out.leading_s;
    out.qend = read_len - out.trailing_s;
    if (out.qend <= out.qstart || out.qend > read_len) {
        out.qend = out.qstart + q_aligned;
    }
    return out.qend > out.qstart && out.qend <= read_len;
}

int32_t local_compute_ref_end(const ReadView& read) {
    const uint32_t* cigar = read.cigar();
    const int32_t n_cigar = read.n_cigar();
    if (!cigar || n_cigar <= 0) {
        return read.pos();
    }
    int32_t ref_pos = read.pos();
    for (int32_t i = 0; i < n_cigar; ++i) {
        const int op = bam_cigar_op(cigar[i]);
        const int32_t len = static_cast<int32_t>(bam_cigar_oplen(cigar[i]));
        if ((bam_cigar_type(op) & 2) != 0) {
            ref_pos += len;
        }
    }
    return ref_pos;
}

bool local_normalized_primary_alignment(
    const ReadView& read,
    const std::string& chrom,
    LocalNormalizedAlignment& out) {
    out = LocalNormalizedAlignment{};
    LocalQueryInterval q_interval;
    if (!local_bam_to_query_interval(read, q_interval)) {
        return false;
    }
    out.chrom = chrom;
    out.ref_start = read.pos();
    out.ref_end = local_compute_ref_end(read);
    out.qstart = q_interval.qstart;
    out.qend = q_interval.qend;
    out.is_reverse = (read.flag() & BAM_FREVERSE) != 0;
    return out.ref_end > out.ref_start && out.qend > out.qstart;
}

bool local_normalized_sa_alignment(
    const LocalSAEntry& entry,
    int32_t read_len,
    LocalNormalizedAlignment& out) {
    out = LocalNormalizedAlignment{};
    std::vector<LocalCigarOp> ops;
    LocalQueryInterval q_interval;
    int32_t ref_aligned_len = 0;
    if (!parse_local_cigar_ops(entry.cigar, ops) ||
        !local_cigar_to_query_interval(ops, read_len, q_interval, ref_aligned_len)) {
        return false;
    }
    out.chrom = entry.rname;
    out.ref_start = entry.pos - 1;
    out.ref_end = out.ref_start + ref_aligned_len;
    out.qstart = q_interval.qstart;
    out.qend = q_interval.qend;
    out.is_reverse = (entry.strand == '-');
    return out.ref_end > out.ref_start && out.qend > out.qstart;
}

std::vector<int32_t> robust_local_split_insertion_positions(
    const ReadView& read,
    const std::string& chrom) {
    std::vector<int32_t> positions;
    std::string sa_z;
    if (!read.get_string_tag("SA", sa_z) || sa_z.empty()) {
        return positions;
    }

    LocalNormalizedAlignment primary;
    if (!local_normalized_primary_alignment(read, chrom, primary)) {
        return positions;
    }

    for (const LocalSAEntry& entry : parse_local_sa_tag_z(sa_z)) {
        if (entry.rname != chrom) {
            continue;
        }
        LocalNormalizedAlignment mate;
        if (!local_normalized_sa_alignment(entry, read.seq_len(), mate)) {
            continue;
        }
        if (mate.is_reverse != primary.is_reverse) {
            continue;
        }
        const LocalNormalizedAlignment* left = &primary;
        const LocalNormalizedAlignment* right = &mate;
        if (mate.qstart < primary.qstart) {
            left = &mate;
            right = &primary;
        }
        const int32_t query_gap = right->qstart - left->qend;
        const int32_t ref_gap = std::max(0, right->ref_start - left->ref_end);
        const int32_t insertion_len = query_gap - ref_gap;
        if (insertion_len >= kLongInsertionSignalMin) {
            positions.push_back(left->ref_end);
        }
    }

    std::sort(positions.begin(), positions.end());
    positions.erase(std::unique(positions.begin(), positions.end()), positions.end());
    return positions;
}

void append_record_ptrs(
    const std::vector<BufferedRecord>& src,
    std::vector<const bam1_t*>& dst,
    int32_t min_ref_end) {
    dst.reserve(dst.size() + src.size());
    for (const auto& rec : src) {
        if (!rec.record) {
            continue;
        }
        if (rec.ref_end < min_ref_end) {
            continue;
        }
        dst.push_back(rec.record.get());
    }
}
