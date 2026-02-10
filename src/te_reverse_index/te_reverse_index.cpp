#include "te_reverse_index.h"
#include "gate1.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <iostream>
#include <fstream>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <limits>
#include <deque>

namespace placer {

// ============================================================================
// GenomeKmerIndex Implementation
// ============================================================================

GenomeKmerIndex::GenomeKmerIndex(const TEReverseIndexConfig& config)
    : config_(config) {
    kmer_mask_ = (1ULL << (2 * config_.kmer_size)) - 1;
}

GenomeKmerIndex::~GenomeKmerIndex() = default;

// ============================================================================
// extract_minimizer：O(n) 单调队列实现
//
// [修正] INVALID_KMER 不入队，遇到 INVALID 时清空队列
// （包含 N 的窗口不应产出 minimizer）
// ============================================================================

uint64_t GenomeKmerIndex::extract_minimizer(
    const std::string& seq, size_t& min_pos) const {

    min_pos = 0;

    if (seq.size() < config_.kmer_size) {
        return 0;
    }

    const size_t k = config_.kmer_size;
    const size_t w = config_.minimizer_window;
    const size_t n_kmers = seq.size() - k + 1;

    constexpr uint64_t INVALID_KMER = std::numeric_limits<uint64_t>::max();

    // 第一步：预计算所有 k-mer 值
    std::vector<uint64_t> kmer_values(n_kmers, INVALID_KMER);

    uint64_t current_kmer = 0;
    size_t valid_bases = 0;

    for (size_t i = 0; i < seq.size(); ++i) {
        uint8_t code = char_to_2bit(seq[i]);
        if (code > 3) {
            valid_bases = 0;
            current_kmer = 0;
        } else {
            current_kmer = ((current_kmer << 2) | code) & kmer_mask_;
            valid_bases++;
        }

        if (valid_bases >= k) {
            kmer_values[i - k + 1] = current_kmer;
        }
    }

    // 第二步：单调队列
    // [修正] INVALID_KMER 不入队；遇到 INVALID 清空队列
    uint64_t global_min = INVALID_KMER;
    size_t global_min_pos = 0;

    std::deque<std::pair<uint64_t, size_t>> dq;

    for (size_t i = 0; i < n_kmers; ++i) {
        uint64_t val = kmer_values[i];

        if (val == INVALID_KMER) {
            // [修正] 含 N 的 kmer：清空队列
            // 任何包含此位置的窗口都不应产出 minimizer
            dq.clear();
            continue;
        }

        // 维护单调递增队列（队首最小）
        while (!dq.empty() && dq.back().first >= val) {
            dq.pop_back();
        }
        dq.push_back({val, i});

        // 移除超出窗口的元素
        while (!dq.empty() && i >= w && dq.front().second <= i - w) {
            dq.pop_front();
        }

        // 窗口已满且队列非空时，队首是当前窗口的 minimizer
        if (i >= w - 1 && !dq.empty()) {
            if (dq.front().first < global_min) {
                global_min = dq.front().first;
                global_min_pos = dq.front().second;
            }
        }
    }

    min_pos = global_min_pos;
    return (global_min == INVALID_KMER) ? 0 : global_min;
}

// ============================================================================
// add_sequence
//
// [修正] Reservoir sampling：total_count 在每次看到 kmer 时都递增，
// 不仅仅在 reservoir 满了之后。否则替换概率严重偏置。
// ============================================================================

void GenomeKmerIndex::add_sequence(
    const std::string& seq, int32_t chrom_tid, int64_t offset) {

    if (seq.size() < config_.kmer_size) {
        return;
    }

    if (chrom_tid >= static_cast<int32_t>(chrom_starts_.size())) {
        chrom_starts_.resize(chrom_tid + 1, offset);
    }

    const size_t max_positions =
        static_cast<size_t>(config_.max_genome_hits_per_kmer);

    // RNG：用 chrom_tid + offset 做 seed，保证可复现
    std::mt19937_64 rng(
        static_cast<uint64_t>(chrom_tid) * 2654435761ULL +
        static_cast<uint64_t>(offset));

    uint64_t current_kmer = 0;
    size_t valid_bases = 0;

    for (size_t i = 0; i < seq.size(); ++i) {
        uint8_t code = char_to_2bit(seq[i]);
        if (code > 3) {
            valid_bases = 0;
            current_kmer = 0;
            continue;
        }

        current_kmer = ((current_kmer << 2) | code) & kmer_mask_;
        valid_bases++;

        if (valid_bases < config_.kmer_size) {
            continue;
        }

        size_t kmer_start = i - config_.kmer_size + 1;

        // [修正] 使用 uint64_t 编码，避免符号扩展问题
        uint64_t encoded_pos =
            (static_cast<uint64_t>(chrom_tid) << 40) |
            static_cast<uint64_t>(offset + static_cast<int64_t>(kmer_start));

        auto& positions = kmer_to_positions_[current_kmer];

        // [修正] Reservoir sampling：
        // total_count 在每次看到 kmer 时都递增（不管 reservoir 是否已满）
        auto& total_count = kmer_total_counts_[current_kmer];
        total_count++;

        if (total_count <= max_positions) {
            // reservoir 未满：直接 push
            positions.push_back(static_cast<int64_t>(encoded_pos));
        } else {
            // reservoir 已满：以 max_positions / total_count 概率替换
            std::uniform_int_distribution<size_t> dist(0, total_count - 1);
            size_t j = dist(rng);
            if (j < max_positions) {
                positions[j] = static_cast<int64_t>(encoded_pos);
            }
        }
    }

    total_genome_bases_ += seq.size();
}

// ============================================================================
// build_from_fasta
// ============================================================================

bool GenomeKmerIndex::build_from_fasta(const std::string& fasta_path) {
    std::ifstream file(fasta_path);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open " << fasta_path << "\n";
        return false;
    }

    std::string line;
    std::string current_seq;
    int32_t chrom_tid = -1;

    while (std::getline(file, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (!current_seq.empty() && chrom_tid >= 0) {
                add_sequence(current_seq, chrom_tid, 0);
            }

            chrom_tid++;

            std::string header = line.substr(1);
            size_t space_pos = header.find_first_of(" \t");
            if (space_pos != std::string::npos) {
                header = header.substr(0, space_pos);
            }
            chrom_names_.push_back(header);

            current_seq.clear();
        } else {
            current_seq += line;
        }
    }

    if (!current_seq.empty() && chrom_tid >= 0) {
        add_sequence(current_seq, chrom_tid, 0);
    }

    std::cout << "GenomeKmerIndex built: "
              << total_genome_bases_ << " bases, "
              << kmer_to_positions_.size() << " unique k-mers, "
              << chrom_names_.size() << " chromosomes\n";

    return true;
}

void GenomeKmerIndex::build_from_sequences(
    const std::vector<std::string>& sequences) {

    for (size_t i = 0; i < sequences.size(); ++i) {
        add_sequence(sequences[i], static_cast<int32_t>(i), 0);
    }

    total_genome_bases_ = 0;
    for (const auto& seq : sequences) {
        total_genome_bases_ += seq.size();
    }

    std::cout << "GenomeKmerIndex built from " << sequences.size()
              << " sequences: " << total_genome_bases_ << " bases, "
              << kmer_to_positions_.size() << " unique k-mers\n";
}

// ============================================================================
// query
//
// 投票 key = probe 对齐起点 = genome_kmer_pos - kmer_offset_in_probe
// sort + 线性计数
// [修正] 全程使用 uint64_t 编码，避免 signed/unsigned 混用问题
// ============================================================================

std::vector<std::pair<int32_t, int32_t>> GenomeKmerIndex::query(
    const std::string& sequence) const {

    std::vector<std::pair<int32_t, int32_t>> results;

    if (sequence.size() < config_.kmer_size) {
        return results;
    }

    const size_t k = config_.kmer_size;

    // [修正] 使用 uint64_t 统一编码
    std::vector<uint64_t> origin_hits;
    origin_hits.reserve(64);

    uint64_t current_kmer = 0;
    size_t valid_bases = 0;

    for (size_t i = 0; i < sequence.size(); ++i) {
        uint8_t code = char_to_2bit(sequence[i]);
        if (code > 3) {
            valid_bases = 0;
            current_kmer = 0;
            continue;
        }

        current_kmer = ((current_kmer << 2) | code) & kmer_mask_;
        valid_bases++;

        if (valid_bases < k) continue;

        size_t kmer_offset = i - k + 1;

        auto it = kmer_to_positions_.find(current_kmer);
        if (it == kmer_to_positions_.end()) continue;

        for (int64_t encoded : it->second) {
            // 解码：(tid << 40) | position
            uint64_t uencoded = static_cast<uint64_t>(encoded);
            uint64_t tid = uencoded >> 40;
            uint64_t pos = uencoded & ((1ULL << 40) - 1);

            // 计算 probe 对齐起点
            uint64_t origin_pos = (pos >= kmer_offset) ? (pos - kmer_offset) : 0;
            if (origin_pos == 0 && pos < kmer_offset) continue;

            // 重新编码：(tid << 40) | origin_pos
            uint64_t encoded_origin = (tid << 40) | origin_pos;
            origin_hits.push_back(encoded_origin);
        }
    }

    if (origin_hits.empty()) return results;

    std::sort(origin_hits.begin(), origin_hits.end());

    uint64_t prev = origin_hits[0];
    int count = 1;

    auto emit = [&](uint64_t enc, int cnt) {
        if (cnt >= config_.min_genome_hits) {
            uint64_t tid = enc >> 40;
            uint64_t pos = enc & ((1ULL << 40) - 1);
            results.emplace_back(static_cast<int32_t>(tid), static_cast<int32_t>(pos));
        }
    };

    for (size_t i = 1; i < origin_hits.size(); ++i) {
        if (origin_hits[i] == prev) {
            count++;
        } else {
            emit(prev, count);
            prev = origin_hits[i];
            count = 1;
        }
    }
    emit(prev, count);

    std::sort(results.begin(), results.end());

    return results;
}

// ============================================================================
// query_hit_count
// ============================================================================

int GenomeKmerIndex::query_hit_count(const std::string& sequence) const {

    if (sequence.size() < config_.kmer_size) {
        return 0;
    }

    std::unordered_set<int64_t> unique_positions;

    uint64_t current_kmer = 0;
    size_t valid_bases = 0;

    for (size_t i = 0; i < sequence.size(); ++i) {
        uint8_t code = char_to_2bit(sequence[i]);
        if (code > 3) {
            valid_bases = 0;
            current_kmer = 0;
            continue;
        }

        current_kmer = ((current_kmer << 2) | code) & kmer_mask_;
        valid_bases++;

        if (valid_bases < config_.kmer_size) continue;

        auto it = kmer_to_positions_.find(current_kmer);
        if (it != kmer_to_positions_.end()) {
            for (int64_t encoded : it->second) {
                unique_positions.insert(encoded);
            }
        }
    }

    return static_cast<int>(unique_positions.size());
}

// ============================================================================
// TEReverseIndex Implementation
// ============================================================================

TEReverseIndex::TEReverseIndex(const TEReverseIndexConfig& config)
    : config_(config),
      gate1_(nullptr, config.gate1_config) {
    genome_index_ = std::make_unique<GenomeKmerIndex>(config);
}

TEReverseIndex::~TEReverseIndex() = default;

bool TEReverseIndex::initialize(const std::string& genome_fasta) {
    return genome_index_->build_from_fasta(genome_fasta);
}

// ============================================================================
// collect_hits_from_probes：统一的 probe → Hit 收集器
// ============================================================================

std::vector<Hit> TEReverseIndex::collect_hits_from_probes(
    const std::vector<ProbeFragment>& probes,
    size_t read_idx) const {

    std::vector<Hit> hits;

    for (const auto& probe : probes) {
        std::string seq_str(probe.sequence);
        auto genome_hits = genome_index_->query(seq_str);

        for (const auto& [tid, pos] : genome_hits) {
            // query() 返回的 pos 是 probe 在 genome 上的对齐起点
            // 再减去 probe 在 read 中的偏移，得到 read 的起点估计
            int32_t read_origin = pos - probe.read_offset;
            if (read_origin >= 0) {
                hits.push_back({tid, read_origin, read_idx, probe.read_offset});
            }
        }
    }

    return hits;
}

// ============================================================================
// cluster_hits：结构化 Hit + 每个 cluster 独立 start/end
//
// [修正] 在 RescuedLocus 中额外存储 genome_hit_positions（genome 坐标），
// 用于 merge 时正确计算 position 和 density。
// ============================================================================

std::vector<RescuedLocus> TEReverseIndex::cluster_hits(
    const std::vector<Hit>& hits) {

    std::vector<RescuedLocus> loci;

    if (hits.empty()) return loci;

    // 按 (tid, pos) 排序
    std::vector<size_t> order(hits.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
        [&hits](size_t a, size_t b) {
            if (hits[a].tid != hits[b].tid)
                return hits[a].tid < hits[b].tid;
            return hits[a].pos < hits[b].pos;
        });

    struct Cluster {
        int32_t tid;
        int32_t start;
        int32_t end;
        std::vector<size_t> hit_indices;
    };

    std::vector<Cluster> clusters;

    {
        size_t first = order[0];
        Cluster c;
        c.tid = hits[first].tid;
        c.start = hits[first].pos;
        c.end = hits[first].pos;
        c.hit_indices.push_back(first);
        clusters.push_back(std::move(c));
    }

    for (size_t oi = 1; oi < order.size(); ++oi) {
        size_t idx = order[oi];
        const auto& h = hits[idx];
        auto& cur = clusters.back();

        if (h.tid != cur.tid ||
            h.pos - cur.end > config_.locus_cluster_radius) {
            Cluster c;
            c.tid = h.tid;
            c.start = h.pos;
            c.end = h.pos;
            c.hit_indices.push_back(idx);
            clusters.push_back(std::move(c));
        } else {
            cur.end = std::max(cur.end, h.pos);
            cur.hit_indices.push_back(idx);
        }
    }

    // 转换为 RescuedLocus
    for (const auto& cluster : clusters) {
        RescuedLocus locus;
        locus.chrom_tid = cluster.tid;

        // [修正] 收集 genome 坐标用于 position/density 计算
        std::vector<int32_t> genome_positions;
        genome_positions.reserve(cluster.hit_indices.size());

        std::unordered_set<size_t> unique_reads;

        for (size_t idx : cluster.hit_indices) {
            genome_positions.push_back(hits[idx].pos);
            unique_reads.insert(hits[idx].read_idx);
            locus.probe_read_indices.push_back(hits[idx].read_idx);
            locus.probe_positions.push_back(hits[idx].probe_offset);
        }

        // [修正] 存储 genome 坐标列表，供 merge 使用
        locus.genome_hit_positions = std::move(genome_positions);

        // 位置取中位数
        auto& gpos = locus.genome_hit_positions;
        std::vector<int32_t> sorted_gpos = gpos;
        std::sort(sorted_gpos.begin(), sorted_gpos.end());
        locus.position = sorted_gpos[sorted_gpos.size() / 2];

        locus.support_reads = static_cast<int>(unique_reads.size());
        locus.total_hits = static_cast<int>(cluster.hit_indices.size());

        // density 使用 cluster 自己的 span
        int32_t span = cluster.end - cluster.start + 1;
        locus.hit_density = (span > 0) ?
            static_cast<double>(locus.total_hits) / span : 0.0;

        // 存储 cluster span 供 merge 使用
        locus.cluster_start = cluster.start;
        locus.cluster_end = cluster.end;

        evaluate_locus(locus);
        loci.push_back(std::move(locus));
    }

    return loci;
}

// ============================================================================
// evaluate_locus
//
// [修正] consistency_score 基于 genome_hit_positions 的方差
// （而非 probe_positions，后者是 read 内偏移，语义不同）
// ============================================================================

void TEReverseIndex::evaluate_locus(RescuedLocus& locus) const {

    double read_score = std::min(
        static_cast<double>(locus.support_reads) / 5.0, 1.0);

    double density_score = std::min(locus.hit_density * 10.0, 1.0);

    // [修正] consistency_score 基于 genome_hit_positions 的集中度
    double consistency_score = 0.5;
    if (locus.genome_hit_positions.size() >= 2) {
        double sum = 0.0, sum_sq = 0.0;
        for (int32_t p : locus.genome_hit_positions) {
            double dp = static_cast<double>(p);
            sum += dp;
            sum_sq += dp * dp;
        }
        double n = static_cast<double>(locus.genome_hit_positions.size());
        double mean = sum / n;
        double variance = (sum_sq / n) - (mean * mean);
        double stddev = std::sqrt(std::max(variance, 0.0));

        // stddev 越小越集中越好
        consistency_score = std::exp(-stddev / 200.0);
        consistency_score = std::clamp(consistency_score, 0.0, 1.0);
    } else if (locus.genome_hit_positions.size() == 1) {
        consistency_score = 0.7;
    }

    locus.placeability_score = read_score * 0.4 +
                               density_score * 0.3 +
                               consistency_score * 0.3;

    locus.passed_filter = (
        locus.support_reads >= config_.min_locus_reads &&
        locus.hit_density >= config_.min_hit_density
    );
}

// ============================================================================
// filter_loci
// ============================================================================

std::vector<RescuedLocus> TEReverseIndex::filter_loci(
    std::vector<RescuedLocus>& loci) {

    std::vector<RescuedLocus> filtered;

    for (auto& locus : loci) {
        if (!locus.passed_filter) {
            stats_.loci_discarded++;
            continue;
        }

        bool merged = false;
        for (auto& existing : filtered) {
            if (existing.chrom_tid == locus.chrom_tid &&
                std::abs(existing.position - locus.position) <
                    config_.locus_cluster_radius) {

                merge_loci(existing, locus);
                merged = true;
                stats_.loci_merged++;
                break;
            }
        }

        if (!merged) {
            filtered.push_back(std::move(locus));
        }
    }

    return filtered;
}

// ============================================================================
// merge_loci
//
// [修正] 完全重写：
// 1. position 和 density 基于 genome_hit_positions（genome 坐标），
//    不再错误地使用 probe_positions（read 内偏移）
// 2. probe_read_indices 和 probe_positions 保持对应关系
// 3. 合并后重新 evaluate
// ============================================================================

void TEReverseIndex::merge_loci(
    RescuedLocus& target, const RescuedLocus& source) const {

    // 1. 合并 genome_hit_positions
    target.genome_hit_positions.insert(
        target.genome_hit_positions.end(),
        source.genome_hit_positions.begin(),
        source.genome_hit_positions.end());

    // 2. 合并 probe_read_indices 和 probe_positions（保持对应关系）
    //    使用 (read_idx, probe_offset) 对去重
    struct ReadProbe {
        size_t read_idx;
        int32_t probe_offset;
        bool operator==(const ReadProbe& o) const {
            return read_idx == o.read_idx && probe_offset == o.probe_offset;
        }
    };
    struct ReadProbeHash {
        size_t operator()(const ReadProbe& rp) const {
            return std::hash<size_t>{}(rp.read_idx) ^
                   (std::hash<int32_t>{}(rp.probe_offset) << 16);
        }
    };

    std::unordered_set<ReadProbe, ReadProbeHash> seen;
    std::vector<size_t> merged_read_indices;
    std::vector<int32_t> merged_probe_positions;

    auto add_pair = [&](size_t ridx, int32_t poff) {
        ReadProbe rp{ridx, poff};
        if (seen.insert(rp).second) {
            merged_read_indices.push_back(ridx);
            merged_probe_positions.push_back(poff);
        }
    };

    for (size_t i = 0; i < target.probe_read_indices.size(); ++i) {
        int32_t poff = (i < target.probe_positions.size()) ?
            target.probe_positions[i] : 0;
        add_pair(target.probe_read_indices[i], poff);
    }
    for (size_t i = 0; i < source.probe_read_indices.size(); ++i) {
        int32_t poff = (i < source.probe_positions.size()) ?
            source.probe_positions[i] : 0;
        add_pair(source.probe_read_indices[i], poff);
    }

    target.probe_read_indices = std::move(merged_read_indices);
    target.probe_positions = std::move(merged_probe_positions);

    // 3. 统计 unique reads
    std::unordered_set<size_t> unique_reads(
        target.probe_read_indices.begin(),
        target.probe_read_indices.end());
    target.support_reads = static_cast<int>(unique_reads.size());

    // 4. total_hits
    target.total_hits = static_cast<int>(target.genome_hit_positions.size());

    // 5. [修正] position：基于 genome_hit_positions 的中位数
    if (!target.genome_hit_positions.empty()) {
        std::vector<int32_t> sorted_gpos = target.genome_hit_positions;
        std::sort(sorted_gpos.begin(), sorted_gpos.end());
        target.position = sorted_gpos[sorted_gpos.size() / 2];
    }

    // 6. [修正] cluster span 和 density：基于 genome_hit_positions
    if (!target.genome_hit_positions.empty()) {
        auto [min_it, max_it] = std::minmax_element(
            target.genome_hit_positions.begin(),
            target.genome_hit_positions.end());
        target.cluster_start = *min_it;
        target.cluster_end = *max_it;
        int32_t span = target.cluster_end - target.cluster_start + 1;
        target.hit_density = (span > 0) ?
            static_cast<double>(target.total_hits) / span : 0.0;
    }

    // 7. 重新 evaluate
    evaluate_locus(target);
}

// ============================================================================
// rescue_batch
// ============================================================================

std::vector<ReadRescueResult> TEReverseIndex::rescue_batch(
    const std::vector<ReadSketch>& reads,
    const std::vector<std::vector<ProbeFragment>>& all_probes,
    size_t start, size_t end) {

    std::vector<ReadRescueResult> results;

    for (size_t i = start; i < end && i < reads.size(); ++i) {
        ReadRescueResult result;
        result.read_idx = i;

        const auto* probes =
            (i < all_probes.size()) ? &all_probes[i] : nullptr;
        if (!probes || probes->empty()) {
            continue;
        }

        auto hits = collect_hits_from_probes(*probes, i);

        result.num_valid_hits = 0;
        for (const auto& h : hits) {
            result.genome_hits.emplace_back(h.tid, h.pos);
            if (h.pos >= 0) result.num_valid_hits++;
        }

        result.rescued =
            (result.num_valid_hits >= config_.min_probe_hits);

        if (result.rescued) {
            stats_.reads_rescued++;
        }

        results.push_back(std::move(result));
    }

    return results;
}

// ============================================================================
// rescue_reads
//
// [修正] 删除粘贴破坏的行
// ============================================================================

std::vector<RescuedLocus> TEReverseIndex::rescue_reads(
    const std::vector<ReadSketch>& reads,
    const std::vector<Component>& existing_components,
    const std::vector<bool>& has_sa_flags) {

    stats_.total_reads_processed = reads.size();

    std::vector<Hit> all_hits;

    for (size_t i = 0; i < reads.size(); ++i) {
        if (i < has_sa_flags.size() && has_sa_flags[i]) {
            stats_.reads_with_sa++;
            continue;
        }

        const auto& read = reads[i];

        // [修正] 使用成员 gate1_，不在循环内构造
        std::vector<ProbeFragment> probes;
        gate1_.extract_probes(read, probes);

        auto hits = collect_hits_from_probes(probes, i);
        all_hits.insert(all_hits.end(), hits.begin(), hits.end());
    }

    auto clustered = cluster_hits(all_hits);
    auto filtered = filter_loci(clustered);

    stats_.total_rescued_loci = filtered.size();

    return filtered;
}

// ============================================================================
// rescue_with_evidence
//
// [修正] stats 计数与 rescue_reads 一致：使用赋值而非 +=
// ============================================================================

std::vector<RescuedLocus> TEReverseIndex::rescue_with_evidence(
    const std::vector<ReadSketch>& reads,
    const std::vector<std::vector<ProbeFragment>>& all_probes) {

    stats_.total_reads_processed = reads.size();  // [修正] 赋值，与 rescue_reads 一致

    std::vector<Hit> all_hits;

    for (size_t i = 0; i < reads.size(); ++i) {
        const auto& read = reads[i];

        if (read.has_sa) {
            stats_.reads_with_sa++;
            continue;
        }

        const auto* probes =
            (i < all_probes.size()) ? &all_probes[i] : nullptr;
        if (!probes || probes->empty()) continue;

        auto hits = collect_hits_from_probes(*probes, i);
        all_hits.insert(all_hits.end(), hits.begin(), hits.end());
    }

    auto clustered = cluster_hits(all_hits);
    auto filtered = filter_loci(clustered);

    stats_.total_rescued_loci = filtered.size();

    return filtered;
}

// ============================================================================
// integrate_with_components
// ============================================================================

std::vector<int32_t> TEReverseIndex::integrate_with_components(
    std::vector<Component>& components,
    const std::vector<ReadSketch>& reads) {

    std::vector<int32_t> updated_ids;

    std::vector<bool> has_sa(reads.size(), false);
    for (size_t i = 0; i < reads.size(); ++i) {
        has_sa[i] = reads[i].has_sa;
    }

    auto rescued = rescue_reads(reads, components, has_sa);

    for (const auto& locus : rescued) {
        bool integrated = false;

        for (size_t i = 0; i < components.size(); ++i) {
            auto& comp = components[i];

            for (const auto& existing : comp.locus_set) {
                if (existing.chrom_tid == locus.chrom_tid &&
                    std::abs(existing.pos - locus.position) <
                        config_.locus_cluster_radius) {

                    bool exists = false;
                    for (const auto& ex : comp.locus_set) {
                        if (ex.chrom_tid == locus.chrom_tid &&
                            std::abs(ex.pos - locus.position) <
                                config_.locus_cluster_radius) {
                            exists = true;
                            break;
                        }
                    }

                    if (!exists) {
                        LocusCandidate candidate;
                        candidate.chrom_tid = locus.chrom_tid;
                        candidate.pos = locus.position;
                        candidate.score = locus.placeability_score;
                        candidate.support_reads = locus.support_reads;
                        comp.locus_set.push_back(candidate);
                        updated_ids.push_back(comp.id);
                    }

                    integrated = true;
                    break;
                }
            }

            if (integrated) break;
        }
    }

    return updated_ids;
}

void TEReverseIndex::add_external_loci(
    const std::vector<std::pair<int32_t, int32_t>>& loci) {
    for (const auto& [tid, pos] : loci) {
        (void)tid;
        (void)pos;
    }
}

// ============================================================================
// Utility: cluster_positions
// ============================================================================

std::vector<std::vector<std::pair<int32_t, int32_t>>> cluster_positions(
    const std::vector<std::pair<int32_t, int32_t>>& positions,
    int32_t radius) {

    std::vector<std::vector<std::pair<int32_t, int32_t>>> clusters;

    if (positions.empty()) return clusters;

    auto sorted = positions;
    std::sort(sorted.begin(), sorted.end(),
        [](const auto& a, const auto& b) {
            if (a.first != b.first) return a.first < b.first;
            return a.second < b.second;
        });

    std::vector<std::pair<int32_t, int32_t>> current_cluster;
    int32_t current_tid = sorted[0].first;
    int32_t cluster_end = sorted[0].second;

    for (const auto& [tid, pos] : sorted) {
        if (tid != current_tid || pos - cluster_end > radius) {
            if (!current_cluster.empty()) {
                clusters.push_back(std::move(current_cluster));
                current_cluster.clear();
            }
            current_tid = tid;
            cluster_end = pos;
        }
        current_cluster.push_back({tid, pos});
        cluster_end = std::max(cluster_end, pos);
    }

    if (!current_cluster.empty()) {
        clusters.push_back(std::move(current_cluster));
    }

    return clusters;
}

}  // namespace placer


