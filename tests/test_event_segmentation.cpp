#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include <unistd.h>

#define private public
#include "pipeline.h"
#undef private

namespace {

std::string make_temp_path(const std::string& stem, const std::string& ext) {
    const long long pid = static_cast<long long>(::getpid());
    return "/tmp/" + stem + "_" + std::to_string(pid) + ext;
}

std::string build_reference(size_t len) {
    static const char kBases[] = {'A', 'C', 'G', 'T'};
    uint32_t state = 17u;
    std::string out;
    out.reserve(len);
    for (size_t i = 0; i < len; ++i) {
        state = (state * 1103515245u) + 12345u;
        out.push_back(kBases[(state >> 16) & 3u]);
    }
    return out;
}

char complement_base(char c) {
    switch (c) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        case 'a': return 'T';
        case 'c': return 'G';
        case 'g': return 'C';
        case 't': return 'A';
        default: return 'N';
    }
}

std::string reverse_complement(std::string seq) {
    std::string out;
    out.reserve(seq.size());
    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        out.push_back(complement_base(*it));
    }
    return out;
}

bool reference_edit_identity_if_at_least(
    const std::string& lhs,
    const std::string& rhs,
    int32_t max_edits,
    double& identity_out) {
    identity_out = 0.0;
    const int32_t n = static_cast<int32_t>(lhs.size());
    const int32_t m = static_cast<int32_t>(rhs.size());
    if (n <= 0 || m <= 0 || max_edits < 0) {
        return false;
    }
    if (std::abs(n - m) > max_edits) {
        return false;
    }

    std::vector<int32_t> prev(static_cast<size_t>(m + 1), 0);
    std::vector<int32_t> curr(static_cast<size_t>(m + 1), 0);
    for (int32_t j = 0; j <= m; ++j) {
        prev[static_cast<size_t>(j)] = j;
    }
    for (int32_t i = 1; i <= n; ++i) {
        curr[0] = i;
        for (int32_t j = 1; j <= m; ++j) {
            const int32_t sub_cost =
                prev[static_cast<size_t>(j - 1)] +
                ((lhs[static_cast<size_t>(i - 1)] == rhs[static_cast<size_t>(j - 1)]) ? 0 : 1);
            const int32_t del_cost = prev[static_cast<size_t>(j)] + 1;
            const int32_t ins_cost = curr[static_cast<size_t>(j - 1)] + 1;
            curr[static_cast<size_t>(j)] = std::min({sub_cost, del_cost, ins_cost});
        }
        std::swap(prev, curr);
    }

    const int32_t dist = prev[static_cast<size_t>(m)];
    if (dist > max_edits) {
        return false;
    }
    const double denom = static_cast<double>(std::max(n, m));
    identity_out = std::clamp(1.0 - (static_cast<double>(dist) / denom), 0.0, 1.0);
    return true;
}

std::string deterministic_sequence(size_t len, uint32_t seed) {
    static const char kBases[] = {'A', 'C', 'G', 'T', 'N'};
    std::string out;
    out.reserve(len);
    uint32_t state = seed;
    for (size_t i = 0; i < len; ++i) {
        state = (state * 1664525u) + 1013904223u;
        out.push_back(kBases[(state >> 24) % 5u]);
    }
    return out;
}

void write_fasta(const std::string& path, const std::string& seq) {
    std::ofstream out(path);
    assert(out.is_open());
    out << ">chr1\n";
    out << seq << '\n';
}

placer::EventSegmentation run_segmentation(
    const std::string& fasta_path,
    const std::string& consensus_seq,
    int32_t bp_left,
    int32_t bp_right) {
    placer::PipelineConfig cfg;
    cfg.reference_fasta_path = fasta_path;
    placer::Pipeline pipeline(cfg, nullptr);

    placer::ComponentCall component;
    component.chrom = "chr1";
    component.tid = 0;

    placer::EventReadEvidence evidence;
    evidence.bp_left = bp_left;
    evidence.bp_right = bp_right;

    placer::EventConsensus consensus;
    consensus.consensus_seq = consensus_seq;
    consensus.consensus_len = static_cast<int32_t>(consensus_seq.size());
    consensus.input_event_reads = 3;
    consensus.qc_pass = true;
    consensus.qc_reason = "PASS_EVENT_CONSENSUS";

    return pipeline.segment_event_consensus(component, evidence, consensus);
}

}  // namespace

int main() {
    using namespace placer;

    const std::string fasta_path = make_temp_path("placer_event_segmentation", ".fa");

    {
        const std::vector<std::pair<std::string, std::string>> edge_cases = {
            {"", "ACGT"},
            {"ACGT", ""},
            {"ACGT", "ACGT"},
            {"ACGT", "ACGA"},
            {"ACGT", "ACGTA"},
            {"ACGTA", "ACGT"},
            {"NNNNACGT", "NNNTACGT"},
            {std::string(128, 'A'), std::string(128, 'A')},
            {std::string(128, 'A'), std::string(128, 'C')},
        };
        for (const auto& example : edge_cases) {
            for (int32_t max_edits = -1; max_edits <= 20; ++max_edits) {
                double expected_identity = 0.0;
                double fast_identity = 0.0;
                const bool expected = reference_edit_identity_if_at_least(
                    example.first,
                    example.second,
                    max_edits,
                    expected_identity);
                const bool fast = Pipeline::edit_identity_if_at_least_fast(
                    example.first,
                    example.second,
                    max_edits,
                    fast_identity);
                assert(fast == expected);
                assert(std::abs(fast_identity - expected_identity) < 1e-12);
            }
        }

        for (int32_t lhs_len = 1; lhs_len <= 128; ++lhs_len) {
            for (int32_t delta = -3; delta <= 3; ++delta) {
                const int32_t rhs_len = lhs_len + delta;
                if (rhs_len <= 0 || rhs_len > 128) {
                    continue;
                }
                const std::string lhs = deterministic_sequence(
                    static_cast<size_t>(lhs_len),
                    static_cast<uint32_t>(17 + lhs_len * 13 + delta));
                const std::string rhs = deterministic_sequence(
                    static_cast<size_t>(rhs_len),
                    static_cast<uint32_t>(91 + lhs_len * 7 + delta * 19));
                for (int32_t max_edits = 0; max_edits <= 20; ++max_edits) {
                    double expected_identity = 0.0;
                    double fast_identity = 0.0;
                    const bool expected = reference_edit_identity_if_at_least(
                        lhs,
                        rhs,
                        max_edits,
                        expected_identity);
                    const bool fast = Pipeline::edit_identity_if_at_least_fast(
                        lhs,
                        rhs,
                        max_edits,
                        fast_identity);
                    assert(fast == expected);
                    assert(std::abs(fast_identity - expected_identity) < 1e-12);
                }
            }
        }
    }

    {
        const std::string reference = build_reference(600);
        write_fasta(fasta_path, reference);

        const int32_t bp = 220;
        const std::string left = reference.substr(150, 70);
        const std::string insert = "TTGGAACCTTGGAACCTTGGAACC";
        const std::string right = reference.substr(220, 70);
        const EventSegmentation segmentation = run_segmentation(
            fasta_path,
            left + insert + right,
            bp,
            bp);
        assert(segmentation.pass);
        assert(segmentation.qc_reason == "PASS_EVENT_SEGMENTATION");
        assert(segmentation.left_flank_seq == left);
        assert(segmentation.insert_seq == insert);
        assert(segmentation.right_flank_seq == right);
        assert(segmentation.left_ref_start == 150);
        assert(segmentation.left_ref_end == 220);
        assert(segmentation.right_ref_start == 220);
        assert(segmentation.right_ref_end == 290);
        assert(segmentation.left_flank_align_len == 70);
        assert(segmentation.right_flank_align_len == 70);
        assert(segmentation.left_flank_identity >= 0.99);
        assert(segmentation.right_flank_identity >= 0.99);

        PipelineConfig cfg;
        cfg.reference_fasta_path = fasta_path;
        Pipeline pipeline(cfg, nullptr);
        ComponentCall component;
        component.chrom = "chr1";
        component.tid = 0;
        EventReadEvidence evidence;
        evidence.bp_left = bp;
        evidence.bp_right = bp;
        EventConsensus consensus;
        consensus.consensus_seq = left + insert + right;
        consensus.consensus_len = static_cast<int32_t>(consensus.consensus_seq.size());
        consensus.input_event_reads = 3;
        consensus.qc_pass = true;
        consensus.qc_reason = "PASS_EVENT_CONSENSUS";
        const EventSegmentation instrumented =
            pipeline.segment_event_consensus(component, evidence, consensus);
        assert(instrumented.pass);
        const EventSegmentationSearchStats stats =
            pipeline.last_segmentation_stats();
        assert(stats.paired_searches >= 2);
        assert(stats.seed_bins_total > 0);
        assert(stats.edit_distance_calls > 0);
        assert(stats.edit_distance_cache_misses > 0);
        assert(stats.edit_distance_cache_hits >= 0);
    }

    {
        const std::string reference = build_reference(640);
        write_fasta(fasta_path, reference);

        PipelineConfig cfg;
        cfg.reference_fasta_path = fasta_path;
        Pipeline pipeline(cfg, nullptr);

        ComponentCall component;
        component.chrom = "chr1";
        component.tid = 0;

        const int32_t bp = 260;
        const std::string left = reference.substr(190, 70);
        const std::string insert = "AACCTTGGAACCTTGGAACC";
        const std::string right = reference.substr(260, 70);

        EventReadEvidence evidence;
        evidence.bp_left = bp;
        evidence.bp_right = bp;

        EventConsensus consensus;
        consensus.consensus_seq = left + insert + right;
        consensus.consensus_len = static_cast<int32_t>(consensus.consensus_seq.size());
        consensus.input_event_reads = 3;
        consensus.qc_pass = true;
        consensus.qc_reason = "PASS_EVENT_CONSENSUS";

        std::unordered_map<std::string, EventSegmentation> cache;
        const EventSegmentation first = pipeline.segment_event_consensus_cached(
            component,
            evidence,
            consensus,
            cache);
        assert(first.pass);
        assert(cache.size() == 1);

        cache.begin()->second.insert_seq = "CACHED_INSERT_SENTINEL";
        cache.begin()->second.qc_reason = "CACHED_SEGMENTATION_SENTINEL";
        const EventSegmentation second = pipeline.segment_event_consensus_cached(
            component,
            evidence,
            consensus,
            cache);
        assert(cache.size() == 1);
        assert(second.insert_seq == "CACHED_INSERT_SENTINEL");
        assert(second.qc_reason == "CACHED_SEGMENTATION_SENTINEL");
    }

    {
        std::string reference = build_reference(700);
        const std::string repeated_left = reference.substr(120, 70);
        reference.replace(250, 70, repeated_left);
        write_fasta(fasta_path, reference);

        const int32_t bp = 320;
        const std::string insert = "AACCGGTTAACCGGTTAACCGGTT";
        const std::string right = reference.substr(320, 70);
        const EventSegmentation segmentation = run_segmentation(
            fasta_path,
            repeated_left + insert + right,
            bp,
            bp);
        assert(segmentation.pass);
        assert(segmentation.qc_reason == "PASS_EVENT_SEGMENTATION");
        assert(segmentation.left_ref_start == 250);
        assert(segmentation.left_ref_end == 320);
        assert(segmentation.right_ref_start == 320);
        assert(segmentation.right_ref_end == 390);
    }

    {
        const std::string reference = build_reference(720);
        write_fasta(fasta_path, reference);

        const int32_t bp = 260;
        const std::string left = reference.substr(190, 70);
        const std::string insert = "CCGTTACCGTTACCGTTACCGTTAACCGGTTA";
        const std::string right = reference.substr(260, 70);
        const EventSegmentation segmentation = run_segmentation(
            fasta_path,
            std::string("GG") + left + insert + right,
            bp,
            bp);
        assert(segmentation.pass);
        assert(segmentation.qc_reason == "PASS_EVENT_SEGMENTATION");
        assert(segmentation.left_flank_seq.size() >= left.size());
        assert(segmentation.left_flank_seq.compare(
            segmentation.left_flank_seq.size() - left.size(),
            left.size(),
            left) == 0);
        assert(segmentation.insert_seq == insert);
        assert(segmentation.right_flank_seq == right);
        assert(segmentation.left_ref_end == 260);
        assert(segmentation.right_ref_start == 260);
        assert(segmentation.right_ref_end == 330);
    }

    {
        const std::string reference = build_reference(760);
        write_fasta(fasta_path, reference);

        const int32_t bp = 300;
        const std::string left = reference.substr(230, 70);
        const std::string insert = "TTAACCGGTTAACCGGTTCCAATTAACCGGTT";
        const std::string right = reference.substr(300, 70);
        const EventSegmentation segmentation = run_segmentation(
            fasta_path,
            left + insert + right + std::string("ACGTACGTACGTACGTACGT"),
            bp,
            bp);
        assert(segmentation.pass);
        assert(segmentation.qc_reason == "PASS_EVENT_SEGMENTATION");
        assert(segmentation.left_flank_seq == left);
        assert(segmentation.insert_seq == insert);
        assert(segmentation.right_flank_seq.size() >= right.size());
        assert(segmentation.right_flank_seq.compare(0, right.size(), right) == 0);
        assert(segmentation.left_ref_start == 230);
        assert(segmentation.left_ref_end == 300);
        assert(segmentation.right_ref_start == 300);
        assert(segmentation.right_ref_end >= 370);
    }

    {
        const std::string reference = build_reference(760);
        write_fasta(fasta_path, reference);

        const int32_t bp = 300;
        const std::string left = reference.substr(230, 70);
        const std::string insert = "TTAACCGGTTAACCGGTTCCAATTAACCGGTT";
        const EventSegmentation segmentation = run_segmentation(
            fasta_path,
            left + insert,
            bp,
            bp);
        assert(segmentation.pass);
        assert(segmentation.qc_reason == "PASS_EVENT_SEGMENTATION_ONE_SIDED_LEFT");
        assert(segmentation.left_flank_seq == left);
        assert(segmentation.insert_seq == insert);
        assert(segmentation.left_ref_start == 230);
        assert(segmentation.left_ref_end == 300);
        assert(segmentation.right_flank_align_len == 0);
        assert(segmentation.right_ref_start == 300);
        assert(segmentation.right_ref_end == 300);
    }

    {
        const std::string reference = build_reference(760);
        write_fasta(fasta_path, reference);

        const int32_t bp = 300;
        const std::string insert = "CCAATTAACCGGTTAACCGGTTCCAATTAA";
        const std::string right = reference.substr(300, 70);
        const EventSegmentation segmentation = run_segmentation(
            fasta_path,
            insert + right,
            bp,
            bp);
        assert(segmentation.pass);
        assert(segmentation.qc_reason == "PASS_EVENT_SEGMENTATION_ONE_SIDED_RIGHT");
        assert(segmentation.insert_seq == insert);
        assert(segmentation.right_flank_seq == right);
        assert(segmentation.left_flank_align_len == 0);
        assert(segmentation.left_ref_start == 300);
        assert(segmentation.left_ref_end == 300);
        assert(segmentation.right_ref_start == 300);
        assert(segmentation.right_ref_end == 370);
    }

    {
        const std::string reference = build_reference(1000);
        write_fasta(fasta_path, reference);

        const int32_t bp = 420;
        const int32_t displaced_right_seed = bp + 235;
        const std::string left = reference.substr(340, 80);
        const std::string insert = "AACCGGTTAACCGGTTCCTTAACCGGTTAACCGGTTCCAATTAACCGG";
        const std::string right = reference.substr(420, 80);
        const EventSegmentation segmentation = run_segmentation(
            fasta_path,
            left + insert + right,
            bp,
            displaced_right_seed);
        assert(segmentation.pass);
        assert(segmentation.qc_reason == "PASS_EVENT_SEGMENTATION_TETHERED_FLANK");
        assert(segmentation.insert_seq == insert);
        assert(segmentation.left_ref_end == bp);
        assert(segmentation.right_ref_start == bp);
        assert(segmentation.left_flank_align_len >= 80);
        assert(segmentation.right_flank_align_len >= 80);
    }

    {
        const std::string reference = build_reference(1000);
        write_fasta(fasta_path, reference);

        const int32_t bp = 520;
        const int32_t displaced_left_seed = bp - 235;
        const std::string left = reference.substr(440, 80);
        const std::string insert = "TTAACCGGTTAACCGGTTCCAATTAACCGGTTAACCGGTTCCTTAACC";
        const std::string right = reference.substr(520, 80);
        const EventSegmentation segmentation = run_segmentation(
            fasta_path,
            left + insert + right,
            displaced_left_seed,
            bp);
        assert(segmentation.pass);
        assert(segmentation.qc_reason == "PASS_EVENT_SEGMENTATION_TETHERED_FLANK");
        assert(segmentation.insert_seq == insert);
        assert(segmentation.left_ref_end == bp);
        assert(segmentation.right_ref_start == bp);
        assert(segmentation.left_flank_align_len >= 80);
        assert(segmentation.right_flank_align_len >= 80);
    }

    {
        std::string reference = build_reference(1200);
        const int32_t bp = 420;
        const std::string repeated_right = reference.substr(420, 50);
        reference.replace(470, 50, repeated_right);
        write_fasta(fasta_path, reference);

        const int32_t displaced_right_seed = bp + 260;
        const std::string left = reference.substr(340, 80);
        const std::string insert = "CCAATTAACCGGTTAACCGGTTCCAATTAACCGGTTAACCGGTTCCAATTAACCGGTTAA";
        const EventSegmentation segmentation = run_segmentation(
            fasta_path,
            left + insert + repeated_right,
            bp,
            displaced_right_seed);
        assert(segmentation.pass);
        assert(segmentation.qc_reason == "PASS_EVENT_SEGMENTATION_ONE_SIDED_LEFT");
        assert(segmentation.left_ref_end == bp);
        assert(segmentation.right_flank_align_len == 0);
    }

    {
        const std::string reference = build_reference(820);
        write_fasta(fasta_path, reference);

        const int32_t bp = 360;
        const std::string left = reference.substr(280, 80);
        const std::string insert = "AACCGGTTAACCGGTTCCTTAACCGGTTAACC";
        const std::string right = reference.substr(360, 80);
        const EventSegmentation segmentation = run_segmentation(
            fasta_path,
            reverse_complement(left + insert + right),
            bp,
            bp);
        assert(segmentation.pass);
        assert(segmentation.qc_reason == "PASS_EVENT_SEGMENTATION");
        assert(segmentation.left_flank_seq == left);
        assert(segmentation.insert_seq == insert);
        assert(segmentation.right_flank_seq == right);
        assert(segmentation.left_ref_start == 280);
        assert(segmentation.left_ref_end == 360);
        assert(segmentation.right_ref_start == 360);
        assert(segmentation.right_ref_end == 440);
    }

    {
        std::string reference = build_reference(900);
        const std::string repeated_left = reference.substr(180, 80);
        reference.replace(420, 80, repeated_left);
        write_fasta(fasta_path, reference);

        const int32_t bp = 500;
        const std::string insert = "AACCGGTTAACCGGTTAACCGGTTAACCGGTT";
        const std::string right = reference.substr(500, 80);
        const EventSegmentation segmentation = run_segmentation(
            fasta_path,
            repeated_left + insert + right,
            bp,
            bp);
        assert(segmentation.pass);
        assert(segmentation.left_ref_start == 420);
        assert(segmentation.left_ref_end == 500);
        assert(segmentation.right_ref_start == 500);
    }

    {
        const std::string reference = build_reference(900);
        write_fasta(fasta_path, reference);

        PipelineConfig cfg;
        cfg.reference_fasta_path = fasta_path;
        Pipeline pipeline(cfg, nullptr);

        ComponentCall component;
        component.chrom = "chr1";
        component.tid = 0;

        EventReadEvidence evidence;
        evidence.bp_left = 520;
        evidence.bp_right = 520;
        evidence.alt_struct_reads = 12;
        evidence.ref_span_reads = 0;

        EventConsensus consensus;
        consensus.consensus_seq = std::string(240, 'T') + std::string(120, 'G');
        consensus.consensus_len = static_cast<int32_t>(consensus.consensus_seq.size());
        consensus.input_event_reads = 6;
        consensus.left_anchor_input_reads = 3;
        consensus.right_anchor_input_reads = 3;
        consensus.partial_context_input_reads = 6;
        consensus.qc_pass = true;
        consensus.qc_reason = "PASS_EVENT_CONSENSUS";

        const EventSegmentation segmentation = pipeline.segment_event_consensus(
            component,
            evidence,
            consensus);
        assert(segmentation.pass);
        assert(segmentation.qc_reason == "PASS_EVENT_SEGMENTATION_UNPLACED_INSERT");
        assert(segmentation.insert_seq == consensus.consensus_seq);
        assert(segmentation.left_ref_end == 520);
        assert(segmentation.right_ref_start == 520);
        assert(segmentation.left_flank_align_len == 0);
        assert(segmentation.right_flank_align_len == 0);
    }

    {
        PipelineConfig cfg;
        cfg.min_soft_clip_for_seq_extract = 20;
        cfg.te_softclip_entropy_min = 1.0;
        cfg.te_softclip_kmer_uniqueness_min = 0.10;
        Pipeline pipeline(cfg, nullptr);

        EventReadEvidence evidence;
        evidence.bp_left = 1000;
        evidence.bp_right = 1000;
        evidence.support_qnames = {
            "insert_read",
            "left_clip_a",
            "left_clip_low_complexity",
            "right_clip_a",
            "right_clip_b",
        };

        EventSegmentation segmentation;
        segmentation.pass = true;
        segmentation.insert_seq =
            "ACGTCAGTTCGATCGATCGGATCCGATCGTACGATCGATGCTAGCTAGGATCCGATCG";

        std::vector<InsertionFragment> fragments;

        InsertionFragment insert_fragment;
        insert_fragment.read_id = "insert_read";
        insert_fragment.source = InsertionFragmentSource::kCigarInsertion;
        insert_fragment.ref_junc_pos = 1000;
        insert_fragment.sequence = segmentation.insert_seq;
        fragments.push_back(insert_fragment);

        InsertionFragment left_clip;
        left_clip.read_id = "left_clip_a";
        left_clip.source = InsertionFragmentSource::kClipRefLeft;
        left_clip.ref_junc_pos = 1000;
        left_clip.anchor_len = 80;
        left_clip.nm = 0;
        left_clip.sequence = segmentation.insert_seq.substr(0, 52);
        fragments.push_back(left_clip);

        InsertionFragment low_complexity_left_clip = left_clip;
        low_complexity_left_clip.read_id = "left_clip_low_complexity";
        low_complexity_left_clip.sequence = std::string(60, 'A');
        fragments.push_back(low_complexity_left_clip);

        InsertionFragment right_clip;
        right_clip.read_id = "right_clip_a";
        right_clip.source = InsertionFragmentSource::kClipRefRight;
        right_clip.ref_junc_pos = 1000;
        right_clip.anchor_len = 80;
        right_clip.nm = 0;
        right_clip.sequence = segmentation.insert_seq.substr(
            segmentation.insert_seq.size() - 54);
        fragments.push_back(right_clip);

        InsertionFragment right_clip_with_one_mismatch = right_clip;
        right_clip_with_one_mismatch.read_id = "right_clip_b";
        right_clip_with_one_mismatch.sequence[3] =
            right_clip_with_one_mismatch.sequence[3] == 'A' ? 'C' : 'A';
        fragments.push_back(right_clip_with_one_mismatch);

        const ClipInsertConcordanceEvidence concordance =
            pipeline.analyze_clip_insert_concordance(
                evidence,
                segmentation,
                fragments);

        assert(concordance.pass);
        assert(concordance.full_insert_reads == 1);
        assert(concordance.left_clip_reads == 1);
        assert(concordance.right_clip_reads == 2);
        assert(concordance.max_left_identity >= 0.99);
        assert(concordance.max_right_identity >= 0.98);
        assert(concordance.qc == "PASS_CLIP_INSERT_CONCORDANCE");
    }

    std::remove(fasta_path.c_str());
    std::remove((fasta_path + ".fai").c_str());
    return 0;
}
