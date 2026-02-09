#include "local_realign.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <thread>
#include <future>
#include <execution>
#include <cstdint>
#include <immintrin.h>
#include <emmintrin.h>
#include <cstring>

#if defined(__AVX2__)
#include <avx2intrin.h>
#endif

namespace placer {

// ============================================================================
// Forward declarations for helper functions
// ============================================================================

static int8_t n_base_penalty(char c) {
    // N bases get severe penalty - alignment should avoid them
    switch (c) {
        case 'N': case 'n': return -10;  // Very bad match
        default: return 0;
    }
}

AlignmentResult affine_gap_align(
    std::string_view query,
    std::string_view target,
    const RealignConfig& config);

// ============================================================================
// BitpackedSeq Implementation
// ============================================================================

// ============================================================================
// BitpackedSeq Implementation
// ============================================================================

void BitpackedSeq::encode(std::string_view ascii_seq) {
    size_t n = ascii_seq.size();
    size_ = n;
    data_.resize((n + 15) / 16);

    size_t out_idx = 0;
    uint64_t word = 0;
    int bits = 0;

    for (size_t i = 0; i < n; ++i) {
        uint8_t code = ascii_to_2bit(ascii_seq[i]);
        word |= (static_cast<uint64_t>(code & MASK) << (bits * 2));
        bits++;

        if (bits == 32) {
            data_[out_idx++] = word;
            word = 0;
            bits = 0;
        }
    }

    if (bits > 0) {
        data_[out_idx++] = word;
    }
}

// ============================================================================
// SIMD-Accelerated Operations
// ============================================================================

#if defined(__AVX2__)
int simd_hamming_avx2(const char* a, const char* b, size_t len) {
    int mismatches = 0;
    const __m256i* av = reinterpret_cast<const __m256i*>(a);
    const __m256i* bv = reinterpret_cast<const __m256i*>(b);
    size_t simd_len = len / 32;

    for (size_t i = 0; i < simd_len; ++i) {
        __m256i xm = _mm256_xor_si256(av[i], bv[i]);
        // Count bits set (popcount of xor = mismatches)
        xm = _mm256_cmpeq_epi8(xm, _mm256_setzero_si256());
        // Invert: now 0=match, 255=mismatch
        xm = _mm256_andnot_si256(xm, _mm256_set1_epi32(-1));
        mismatches += __builtin_popcount(_mm256_movemask_epi8(xm));
    }

    // Handle remaining bytes
    for (size_t i = simd_len * 32; i < len; ++i) {
        if (a[i] != b[i]) ++mismatches;
    }

    return mismatches;
}
#else
int simd_hamming_avx2(const char*, const char*, size_t len) {
    return -1; // Not available
}
#endif

int LocalRealigner::hamming_distance(std::string_view a, std::string_view b) {
    if (a.size() != b.size()) return -1;
    if (a.empty()) return 0;

#if defined(__AVX2__)
    int result = simd_hamming_avx2(a.data(), b.data(), a.size());
    if (result >= 0) return result;
#endif

    // Scalar fallback
    int distance = 0;
    const char* pa = a.data();
    const char* pb = b.data();
    const char* end = pa + a.size();

    // Process 8 bytes at a time for better cache utilization
    while (pa + 8 <= end) {
        distance += (pa[0] != pb[0]) + (pa[1] != pb[1]) + (pa[2] != pb[2]) +
                    (pa[3] != pb[3]) + (pa[4] != pb[4]) + (pa[5] != pb[5]) +
                    (pa[6] != pb[6]) + (pa[7] != pb[7]);
        pa += 8;
        pb += 8;
    }

    // Remaining bytes
    while (pa < end) {
        if (*pa != *pb) ++distance;
        ++pa;
        ++pb;
    }

    return distance;
}

// ============================================================================
// Affine Gap Alignment (Needleman-Wunsch with affine penalties)
// ============================================================================

AlignmentResult affine_gap_align(
    std::string_view query,
    std::string_view target,
    const RealignConfig& config) {

    AlignmentResult result;
    const int8_t match = config.match;
    const int8_t mismatch = config.mismatch;
    const int8_t gap_open = config.gap_open;
    const int8_t gap_extend = config.gap_extend;

    int qlen = static_cast<int>(query.size());
    int tlen = static_cast<int>(target.size());

    if (qlen == 0 || tlen == 0) return result;

    // Use smaller dimension for columns to save memory
    bool swap = qlen > tlen;
    std::string_view q = swap ? target : query;
    std::string_view t = swap ? query : target;
    int m = qlen, n = tlen;
    if (swap) std::swap(m, n);

    // Score matrices: M=match/mismatch, X=gap in t (deletion), Y=gap in q (insertion)
    const int NEG_INF = INT_MIN / 4;
    std::vector<int> M_prev(m + 1, NEG_INF), M_curr(m + 1, NEG_INF);
    std::vector<int> X_prev(m + 1, NEG_INF), X_curr(m + 1, NEG_INF);
    std::vector<int> Y_prev(m + 1, NEG_INF), Y_curr(m + 1, NEG_INF);

    // Initialize first row
    M_prev[0] = 0;
    X_prev[0] = NEG_INF;
    Y_prev[0] = NEG_INF;
    for (int i = 1; i <= m; ++i) {
        M_prev[i] = NEG_INF;
        X_prev[i] = gap_open + i * gap_extend;
        Y_prev[i] = NEG_INF;
    }

    int best_score = INT_MIN;
    int best_i = 0, best_j = 0;

    for (int j = 1; j <= n; ++j) {
        char tc = t[j - 1];

        M_curr[0] = NEG_INF;
        X_curr[0] = NEG_INF;
        Y_curr[0] = gap_open + j * gap_extend;

        for (int i = 1; i <= m; ++i) {
            char qc = q[i - 1];

            // Handle N bases - treat as soft match with penalty
            bool is_n_q = (qc == 'N' || qc == 'n');
            bool is_n_t = (tc == 'N' || tc == 'N');
            int8_t n_penalty = (is_n_q || is_n_t) ? n_base_penalty('N') : 0;

            int8_t match_score = (qc == tc || is_n_q || is_n_t) ? match : mismatch;
            int score_m = M_prev[i - 1] + match_score + n_penalty;
            int score_x = X_prev[i - 1] + match_score + n_penalty;  // Extend alignment into gap
            int score_y = Y_prev[i - 1] + match_score + n_penalty;

            // M = max of all three coming from diagonal
            M_curr[i] = std::max({score_m, score_x, score_y});

            // X = gap in t (deletion from query): max from M (open) or X (extend)
            X_curr[i] = std::max(
                M_prev[i] + gap_open + gap_extend,
                X_prev[i] + gap_extend
            );

            // Y = gap in q (insertion to target): max from M (open) or Y (extend)
            Y_curr[i] = std::max(
                M_curr[i - 1] + gap_open + gap_extend,
                Y_curr[i - 1] + gap_extend
            );

            int cell_score = std::max({M_curr[i], X_curr[i], Y_curr[i]});
            if (cell_score > best_score) {
                best_score = cell_score;
                if (swap) {
                    best_i = j - 1;
                    best_j = i - 1;
                } else {
                    best_i = i - 1;
                    best_j = j - 1;
                }
            }
        }

        std::swap(M_prev, M_curr);
        std::swap(X_prev, X_curr);
        std::swap(Y_prev, Y_curr);
    }

    result.score = static_cast<double>(best_score);
    result.normalized_score = best_score / static_cast<double>(std::min(qlen, tlen));

    // Estimate matches/mismatches from score
    int estimated_matches = best_score > 0 ?
        static_cast<int>(best_score / (match - std::abs(mismatch) * 0.5)) : 0;
    estimated_matches = std::max(0, std::min(estimated_matches, qlen));

    result.matches = estimated_matches;
    result.mismatches = qlen - estimated_matches;
    result.identity = static_cast<float>(estimated_matches) / qlen;
    result.has_indel = (qlen != tlen);
    result.query_start = 0;
    result.query_end = qlen;
    result.target_start = 0;
    result.target_end = tlen;

    // Count gaps
    result.gap_opens = (result.has_indel ? 1 : 0);

    // Gap extensions approximated
    result.gap_extensions = std::abs(qlen - tlen);

    return result;
}

// ============================================================================
// GenomeAccessor Implementation
// ============================================================================

GenomeAccessor::GenomeAccessor(std::string_view fasta_path)
    : fasta_path_(fasta_path) {
    load_fai(fasta_path);
}

GenomeAccessor::GenomeAccessor(GenomeAccessor&& other) noexcept
    : fasta_path_(std::move(other.fasta_path_))
    , index_(std::move(other.index_)) {}

GenomeAccessor& GenomeAccessor::operator=(GenomeAccessor&& other) noexcept {
    if (this != &other) {
        fasta_path_ = std::move(other.fasta_path_);
        index_ = std::move(other.index_);
    }
    return *this;
}

bool GenomeAccessor::load_fai(std::string_view fasta_path) {
    std::string fai_path(fasta_path);
    fai_path += ".fai";

    std::ifstream file(fai_path);
    if (!file.is_open()) {
        return false;
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        IndexEntry entry;
        std::istringstream iss(line);
        iss >> entry.name >> entry.length >> entry.offset >> entry.line_len >> entry.line_blen;
        index_.push_back(entry);
    }

    return !index_.empty();
}

std::optional<SeqView> GenomeAccessor::fetch_region(int chrom_tid, int32_t start, int32_t end) const {
    if (chrom_tid < 0 || chrom_tid >= static_cast<int>(index_.size())) {
        return std::nullopt;
    }

    const auto& idx = index_[chrom_tid];
    int64_t file_offset = idx.offset + start + (start / idx.line_blen);
    int64_t fetch_len = (end - start) + ((end - start) / idx.line_blen) + 2;

    std::string chunk(static_cast<size_t>(fetch_len), '\0');
    std::ifstream file(fasta_path_);
    file.seekg(file_offset);
    file.read(chunk.data(), fetch_len);

    // Remove newlines and convert to uppercase
    static thread_local std::string cleaned;
    cleaned.clear();
    cleaned.reserve(end - start);

    for (char c : chunk) {
        if (c == '\n' || c == '\r') continue;
        cleaned.push_back(std::toupper(c));
        if (static_cast<int>(cleaned.size()) >= end - start) break;
    }

    return SeqView(cleaned.data(), cleaned.size());
}

std::optional<SeqView> GenomeAccessor::fetch(int chrom_tid, int32_t start, int32_t end) const {
    // Validate range
    if (chrom_tid < 0) return std::nullopt;
    if (start < 0) start = 0;

    int64_t chrom_len = 0;
    if (chrom_tid < static_cast<int>(index_.size())) {
        chrom_len = static_cast<int64_t>(index_[chrom_tid].length);
    }
    if (end > chrom_len) end = static_cast<int32_t>(chrom_len);
    if (start >= end) return std::nullopt;

    return fetch_region(chrom_tid, start, end);
}

std::optional<BitpackedSeq> GenomeAccessor::fetch_packed(int chrom_tid) const {
    auto view = fetch(chrom_tid, 0, static_cast<int32_t>(index_[chrom_tid].length));
    if (!view.has_value()) return std::nullopt;

    BitpackedSeq seq(view->sv());
    return seq;
}

// ============================================================================
// LocalRealigner Implementation
// ============================================================================

LocalRealigner::LocalRealigner(RealignConfig config)
    : config_(std::move(config)) {}

// ============================================================================
// Alignment Methods
// ============================================================================

AlignmentResult LocalRealigner::align_sequences(
    std::string_view query,
    std::string_view target,
    const RealignConfig& config) {

    AlignmentResult result;

    if (query.empty() || target.empty()) {
        return result;
    }

    // Use affine gap alignment for proper indel handling
    result = affine_gap_align(query, target, config);

    // Check if alignment meets minimum thresholds
    if (!result.is_valid()) {
        // Return invalid result
        result.score = -1e9;
        return result;
    }

    // Validate against config thresholds
    if (result.identity < config.min_identity ||
        result.normalized_score < config.min_normalized_score) {
        result.score = -1e9;
        return result;
    }

    // Calculate equivalent MAPQ based on delta score and quality
    if (result.matches + result.mismatches > 0) {
        float id = result.identity;
        result.mapq_equiv = static_cast<float>(-10.0 * std::log10(1.0 - id + 1e-6));
        result.mapq_equiv = std::min(result.mapq_equiv, 60.0f);
    }

    return result;
}

std::vector<uint32_t> LocalRealigner::find_seeds(
    std::string_view seq,
    int kmer_size,
    int stride) {

    std::vector<uint32_t> seeds;
    if (seq.size() < static_cast<size_t>(kmer_size)) return seeds;

    uint32_t kmer = 0;
    for (int i = 0; i < kmer_size; ++i) {
        uint8_t code = BitpackedSeq::ascii_to_2bit(seq[i]);
        kmer = (kmer << 2) | (code & 0x03);
    }
    seeds.push_back(kmer);

    uint32_t mask = (1u << (kmer_size * 2)) - 1;
    for (size_t i = kmer_size; i < seq.size(); ++i) {
        if ((i - kmer_size) % static_cast<size_t>(stride) != 0) continue;

        // Rolling hash
        uint8_t outgoing = BitpackedSeq::ascii_to_2bit(seq[i - kmer_size]);
        uint8_t incoming = BitpackedSeq::ascii_to_2bit(seq[i]);

        kmer = ((kmer << 2) | (incoming & 0x03)) & mask;
        seeds.push_back(kmer);
    }

    return seeds;
}

// ============================================================================
// Locus Management
// ============================================================================

void LocalRealigner::populate_locus_set(
    Component& component,
    const std::vector<ReadSketch>& reads) {

    component.locus_set.clear();

    for (size_t read_idx : component.read_indices) {
        if (read_idx >= reads.size()) continue;

        const auto& read = reads[read_idx];

        // Primary alignment
        if (read.pos >= component.start - config_.search_window &&
            read.pos <= component.end + config_.search_window) {
            component.locus_set.push_back({
                read.tid,
                read.pos,
                static_cast<double>(read.mapq),
                1,
                0x01
            });
        }

        // SA splits
        if (read.has_sa) {
            for (const auto& [sa_tid, sa_pos] : read.sa_targets) {
                if (sa_tid != component.chrom_tid) continue;

                int32_t window_start = component.start - config_.search_window;
                int32_t window_end = component.end + config_.search_window;

                if (sa_pos >= window_start && sa_pos <= window_end) {
                    component.locus_set.push_back({
                        sa_tid,
                        sa_pos,
                        static_cast<double>(read.mapq) * 0.8,
                        1,
                        0x02
                    });
                }
            }
        }
    }

    // Sort by position
    std::sort(component.locus_set.begin(), component.locus_set.end(),
        [](const LocusCandidate& a, const LocusCandidate& b) {
            return a.pos < b.pos;
        });

    // Merge nearby loci
    std::vector<LocusCandidate> merged;
    const int merge_distance = 50;

    for (const auto& locus : component.locus_set) {
        if (merged.empty()) {
            merged.push_back(locus);
        } else {
            LocusCandidate& last = merged.back();
            if (std::abs(locus.pos - last.pos) <= merge_distance) {
                last.score = std::max(last.score, locus.score);
                last.support_reads += locus.support_reads;
                last.evidence_mask |= locus.evidence_mask;
            } else {
                merged.push_back(locus);
            }
        }
    }

    component.locus_set = std::move(merged);

    // Limit count
    if (component.locus_set.size() > config_.max_locus_per_component) {
        std::partial_sort(
            component.locus_set.begin(),
            component.locus_set.begin() + config_.max_locus_per_component,
            component.locus_set.end(),
            [](const LocusCandidate& a, const LocusCandidate& b) {
                return a.score > b.score;
            }
        );
        component.locus_set.resize(config_.max_locus_per_component);
    }
}

std::vector<LocusCandidate> LocalRealigner::rank_loci(
    std::vector<LocusCandidate>&& loci,
    const RealignConfig& config) {

    (void)config;  // unused in this implementation

    std::sort(loci.begin(), loci.end(),
        [](const LocusCandidate& a, const LocusCandidate& b) {
            // Score by evidence strength (use __builtin_popcount for C++17)
            double a_strength = a.score * __builtin_popcount(a.evidence_mask);
            double b_strength = b.score * __builtin_popcount(b.evidence_mask);
            return a_strength > b_strength;
        }
    );

    return std::move(loci);
}

// ============================================================================
// Component Processing
// ============================================================================

PlaceabilityReport LocalRealigner::realign_component(
    Component& component,
    const std::vector<ReadSketch>& reads,
    const GenomeAccessor& genome) {

    PlaceabilityReport report;

    // Populate locus set
    populate_locus_set(component, reads);

    if (component.locus_set.empty()) {
        // Use centroid as default
        component.locus_set.push_back({
            component.chrom_tid,
            component.centroid,
            1.0,
            component.read_count,
            0x04
        });
    }

    std::vector<LocusEvidence> all_evidence;

    // Process each read
    for (size_t read_idx : component.read_indices) {
        if (read_idx >= reads.size()) continue;
        if (all_evidence.size() >= config_.max_reads_per_component) break;

        const auto& read = reads[read_idx];

        // Extract flanks (simplified - would use genome accessor in production)
        // For now, skip actual alignment and collect evidence from loci
        for (const auto& locus : component.locus_set) {
            LocusEvidence ev;
            ev.read_idx = read_idx;
            ev.locus_pos = locus.pos;

            // Placeholder scoring based on distance
            double dist = std::abs(static_cast<double>(read.pos) - locus.pos);
            ev.normalized_score = std::max(0.0, 1.0 - dist / config_.search_window);
            ev.total_score = ev.normalized_score * 100;

            ev.is_reverse = read.flag & 0x10;

            // Evidence bits
            if (ev.normalized_score > config_.min_normalized_score) {
                ev.evidence_bits |= 0x01;
            }

            all_evidence.push_back(ev);
        }
    }

    // Calculate placeability
    report = calculate_placeability(all_evidence);

    // Strand counts
    for (const auto& ev : all_evidence) {
        if (ev.is_reverse) ++report.reverse_count;
        else ++report.forward_count;
    }
    report.strand_balanced = (report.forward_count > 0 && report.reverse_count > 0);

    // Determine tier
    report.tier = determine_tier_(report);

    return report;
}

// ============================================================================
// Utility Methods
// ============================================================================

std::pair<SeqView, SeqView> LocalRealigner::extract_flanks_view(
    const ReadSketch& read,
    int flank_length,
    const GenomeAccessor& genome) {

    // Upstream flank
    auto up_opt = genome.fetch(read.tid, std::max<int32_t>(0, read.pos - flank_length), read.pos);
    SeqView up_flank = up_opt.value_or(SeqView());

    // Downstream flank
    auto down_opt = genome.fetch(read.tid, read.end_pos,
        std::min<int32_t>(read.end_pos + flank_length,
            static_cast<int32_t>(genome.get_index(read.tid).length)));
    SeqView down_flank = down_opt.value_or(SeqView());

    return {up_flank, down_flank};
}

PlaceabilityReport LocalRealigner::calculate_placeability(
    const std::vector<LocusEvidence>& evidence) {

    PlaceabilityReport report;

    if (evidence.empty()) {
        return report;
    }

    // Collect scores per locus (threshold = 0.3 as default)
    constexpr double kMinScoreThreshold = 0.3;
    std::vector<double> locus_scores;
    for (const auto& ev : evidence) {
        if (ev.normalized_score >= kMinScoreThreshold) {
            locus_scores.push_back(ev.normalized_score);
        }
    }

    if (locus_scores.empty()) {
        return report;
    }

    std::sort(locus_scores.begin(), locus_scores.end(), std::greater<double>());

    report.best_score = locus_scores[0];
    report.best_normalized = locus_scores[0];

    if (locus_scores.size() > 1) {
        report.second_score = locus_scores[1];
        report.second_normalized = locus_scores[1];
        report.delta_score = report.best_score - report.second_score;
    } else {
        report.delta_score = report.best_score;
    }

    // Confidence based on delta and consistency
    report.confidence = static_cast<float>(
        std::min(1.0, report.delta_score / 30.0) *
        (report.best_normalized > 0.5 ? 1.0 : report.best_normalized * 2.0)
    );

    return report;
}

// ============================================================================
// Thread Pool for Batch Processing
// ============================================================================

namespace {

class ThreadPool {
public:
    explicit ThreadPool(size_t num_threads) {
        size_t n = num_threads > 0 ? num_threads : std::thread::hardware_concurrency();
        for (size_t i = 0; i < n; ++i) {
            workers_.emplace_back([this] {
                while (true) {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(queue_mutex_);
                        condition_.wait(lock, [this] {
                            return stop_ || !tasks_.empty();
                        });
                        if (stop_ && tasks_.empty()) return;
                        task = std::move(tasks_.front());
                        tasks_.pop();
                    }
                    task();
                }
            });
        }
    }

    ~ThreadPool() {
        {
            std::unique_lock<std::mutex> lock(queue_mutex_);
            stop_ = true;
        }
        condition_.notify_all();
        for (auto& t : workers_) {
            if (t.joinable()) t.join();
        }
    }

    template<typename F>
    auto enqueue(F&& f) -> std::future<decltype(f())> {
        using R = decltype(f());
        auto task = std::make_shared<std::packaged_task<R()>>(std::forward<F>(f));
        {
            std::unique_lock<std::mutex> lock(queue_mutex_);
            tasks_.emplace([task]() { (*task)(); });
        }
        condition_.notify_one();
        return task->get_future();
    }

    size_t size() const { return workers_.size(); }

private:
    std::vector<std::thread> workers_;
    std::queue<std::function<void()>> tasks_;
    std::mutex queue_mutex_;
    std::condition_variable condition_;
    bool stop_ = false;
};

}  // anonymous namespace

// ============================================================================
// Internal Methods
// ============================================================================

AlignmentResult LocalRealigner::dispatch_align_(
    std::string_view query,
    std::string_view target) {

    if (config_.use_simd && query.size() >= 32) {
        return simd_align_(query, target);
    }
    return scalar_align_(query, target);
}

AlignmentResult LocalRealigner::simd_align_(
    std::string_view query,
    std::string_view target) {

#if defined(__AVX2__)
    // SIMD-accelerated seed finding
    // For now, use affine gap which is already optimized
    return scalar_align_(query, target);
#else
    return scalar_align_(query, target);
#endif
}

AlignmentResult LocalRealigner::scalar_align_(
    std::string_view query,
    std::string_view target) {

    return align_sequences(query, target, config_);
}

PlaceabilityReport LocalRealigner::process_component_(
    Component& component,
    const std::vector<ReadSketch>& reads,
    const GenomeAccessor& genome) {

    return realign_component(component, reads, genome);
}

std::vector<PlaceabilityReport> LocalRealigner::realign_batch(
    std::vector<Component>& components,
    const std::vector<ReadSketch>& reads,
    const GenomeAccessor& genome) {

    std::vector<PlaceabilityReport> reports;
    reports.reserve(components.size());

    if (components.empty()) return reports;

    // Use thread pool for parallel processing
    ThreadPool pool(config_.num_threads);
    std::vector<std::future<PlaceabilityReport>> futures;
    futures.reserve(components.size());

    for (auto& comp : components) {
        futures.push_back(pool.enqueue([this, &comp, &reads, &genome] {
            return realign_component(comp, reads, genome);
        }));
    }

    // Collect results
    for (auto& f : futures) {
        if (f.valid()) {
            reports.push_back(f.get());
        }
    }

    return reports;
}

}  // namespace placer
