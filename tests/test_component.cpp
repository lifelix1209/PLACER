#include <iostream>
#include <cassert>
#include <vector>
#include "bam_reader.h"
#include "component_builder.h"
#include "task_queue.h"
#include "test_path_utils.h"

using namespace placer;

void test_component_builder_basic() {
    std::cout << "Testing ComponentBuilder (basic)..." << std::endl;

    ComponentBuilder builder;
    ComponentBuilderConfig config;
    config.cluster_gap = 50;
    config.max_cluster_span = 200;
    builder = ComponentBuilder(config);

    // 创建测试 reads
    std::vector<ReadSketch> reads;

    // Read 1: 左断点 ~1000
    ReadSketch r1;
    r1.qname = "read1";
    r1.tid = 0;
    r1.pos = 1000;
    r1.end_pos = 3000;
    r1.cigar_ops = {{'S', 50}, {'M', 2000}};  // 左 clip 50bp
    r1.has_sa = false;
    reads.push_back(r1);

    // Read 2: 靠近 Read 1 的断点 ~1005
    ReadSketch r2;
    r2.qname = "read2";
    r2.tid = 0;
    r2.pos = 1005;
    r2.end_pos = 3500;
    r2.cigar_ops = {{'S', 30}, {'M', 2500}};
    r2.has_sa = true;
    r2.sa_targets = {{0, 8000}};  // SA 指向远端
    reads.push_back(r2);

    // Read 3: 不同位置的断点 ~5000
    ReadSketch r3;
    r3.qname = "read3";
    r3.tid = 0;
    r3.pos = 5000;
    r3.end_pos = 7000;
    r3.cigar_ops = {{'S', 40}, {'M', 2000}};
    r3.has_sa = false;
    reads.push_back(r3);

    auto components = builder.build(reads);

    std::cout << "  Total anchors: " << builder.get_total_anchors() << std::endl;
    std::cout << "  Components: " << builder.get_total_components() << std::endl;

    // 应该形成至少 2 个 components
    assert(components.size() >= 2);
    assert(builder.get_total_anchors() > 0);

    std::cout << "  Components created:" << std::endl;
    for (const auto& comp : components) {
        std::cout << "    Comp " << comp.id << ": [" << comp.start << "-" << comp.end
                  << "], reads=" << comp.read_count << ", anchors=" << comp.anchor_count << std::endl;
    }

    std::cout << "  Basic ComponentBuilder tests passed!" << std::endl;
}

void test_component_builder_multi_chrom() {
    std::cout << "Testing ComponentBuilder (multi-chromosome)..." << std::endl;

    ComponentBuilder builder;

    std::vector<ReadSketch> reads;

    // Read 1: chr1 (tid=0)
    ReadSketch r1;
    r1.qname = "chr1_read";
    r1.tid = 0;
    r1.pos = 1000;
    r1.end_pos = 3000;
    r1.cigar_ops = {{'S', 50}, {'M', 2000}};
    r1.has_sa = false;
    reads.push_back(r1);

    // Read 2: chr2 (tid=1)
    ReadSketch r2;
    r2.qname = "chr2_read";
    r2.tid = 1;
    r2.pos = 1000;
    r2.end_pos = 3000;
    r2.cigar_ops = {{'S', 40}, {'M', 2000}};
    r2.has_sa = false;
    reads.push_back(r2);

    auto components = builder.build(reads);

    std::cout << "  Components: " << components.size() << std::endl;

    // 应该是 2 个 components（不同染色体）
    assert(components.size() == 2);

    // 验证染色体
    assert(components[0].chrom_tid != components[1].chrom_tid);

    std::cout << "  Multi-chromosome tests passed!" << std::endl;
}

void test_component_builder_span_threshold() {
    std::cout << "Testing ComponentBuilder (span threshold)..." << std::endl;

    ComponentBuilder builder;
    ComponentBuilderConfig config;
    config.cluster_gap = 50;
    config.max_cluster_span = 200;  // 限制单断点簇的最大跨度
    builder = ComponentBuilder(config);

    std::vector<ReadSketch> reads;

    // 创建跨度 500bp 的散点（应该触发 breaker）
    for (int i = 0; i < 10; i++) {
        ReadSketch r;
        r.qname = "span_read_" + std::to_string(i);
        r.tid = 0;
        r.pos = 1000 + i * 50;  // 跨度 500bp
        r.end_pos = r.pos + 2000;
        r.cigar_ops = {{'S', 30}, {'M', 2000}};
        r.has_sa = false;
        reads.push_back(r);
    }

    auto components = builder.build(reads);

    std::cout << "  Breakers triggered: " << builder.get_breakers_triggered() << std::endl;
    std::cout << "  Components: " << components.size() << std::endl;

    // 应该触发 breaker
    assert(builder.get_breakers_triggered() >= 1);

    std::cout << "  Span threshold tests passed!" << std::endl;
}

void test_component_builder_ins_event() {
    std::cout << "Testing ComponentBuilder (insertion events)..." << std::endl;

    ComponentBuilder builder;
    ComponentBuilderConfig config;
    config.min_ins_len = 50;
    config.cluster_gap = 50;
    config.max_cluster_span = 200;
    builder = ComponentBuilder(config);

    std::vector<ReadSketch> reads;

    // Read 带有 Large Insertion
    ReadSketch r;
    r.qname = "ins_read";
    r.tid = 0;
    r.pos = 1000;
    r.end_pos = 4000;
    r.cigar_ops = {{'M', 1000}, {'I', 100}, {'M', 2900}};  // 100bp insertion
    r.has_sa = false;
    reads.push_back(r);

    auto components = builder.build(reads);

    std::cout << "  Components: " << components.size() << std::endl;
    std::cout << "  Anchors: " << builder.get_total_anchors() << std::endl;

    // 应该能提取到 insertion anchor
    assert(builder.get_total_anchors() >= 1);

    std::cout << "  Insertion event tests passed!" << std::endl;
}

void test_component_builder_sa_split() {
    std::cout << "Testing ComponentBuilder (SA splits)..." << std::endl;

    ComponentBuilder builder;

    std::vector<ReadSketch> reads;

    // Read 带有多个 SA splits
    ReadSketch r;
    r.qname = "sa_read";
    r.tid = 0;
    r.pos = 1000;
    r.end_pos = 3000;
    r.cigar_ops = {{'M', 2000}};
    r.has_sa = true;
    r.sa_targets = {{1, 5000}, {2, 10000}};  // 两个 SA split
    reads.push_back(r);

    auto components = builder.build(reads);

    std::cout << "  Anchors: " << builder.get_total_anchors() << std::endl;

    // 应该为每个 SA split 生成一个 anchor
    assert(builder.get_total_anchors() == 2);

    std::cout << "  SA split tests passed!" << std::endl;
}

void test_component_builder_real_bam() {
    std::cout << "Testing ComponentBuilder (real BAM data)..." << std::endl;

    ComponentBuilder builder;
    ComponentBuilderConfig config;
    config.min_clip_len = 20;
    config.cluster_gap = 50;
    config.max_cluster_span = 200;
    builder = ComponentBuilder(config);

    const std::string bam_path = placer_test::resolve_test_file(
        "PLACER_TEST_BAM",
        {"test_data/test.bam", "tests/data/test.bam"});
    if (!placer_test::require_path_or_skip(
            bam_path, "BAM fixture", "PLACER_TEST_BAM")) {
        return;
    }

    BamReader reader(bam_path);
    if (!reader.is_valid()) {
        std::cout << "  [SKIP] Failed to open BAM fixture: " << bam_path
                  << std::endl;
        return;
    }

    std::vector<ReadSketch> reads;
    int64_t count = 0;

    auto callback = [&reads, &count](const ReadSketch& read) {
        reads.push_back(read);
        count++;
    };

    reader.stream(callback);
    std::cout << "  Loaded " << count << " reads from BAM" << std::endl;

    auto components = builder.build(reads);

    std::cout << "  Total anchors: " << builder.get_total_anchors() << std::endl;
    std::cout << "  Components: " << builder.get_total_components() << std::endl;
    std::cout << "  Breakers triggered: " << builder.get_breakers_triggered() << std::endl;

    assert(components.size() > 0);

    for (const auto& comp : components) {
        std::cout << "    Component " << comp.id << ": [" << comp.start << "-" << comp.end
                  << "], reads=" << comp.read_count << ", anchors=" << comp.anchor_count << std::endl;
    }

    std::cout << "  Real BAM tests passed!" << std::endl;
}

void test_component_builder_empty() {
    std::cout << "Testing ComponentBuilder (empty input)..." << std::endl;

    ComponentBuilder builder;

    std::vector<ReadSketch> reads;
    auto components = builder.build(reads);

    assert(components.empty());
    assert(builder.get_total_anchors() == 0);
    assert(builder.get_total_components() == 0);

    std::cout << "  Empty input tests passed!" << std::endl;
}

void test_component_builder_no_signals() {
    std::cout << "Testing ComponentBuilder (reads without signals)..." << std::endl;

    ComponentBuilder builder;

    std::vector<ReadSketch> reads;

    // Read 只有 M operation，无 clip/SA/ins
    ReadSketch r;
    r.qname = "clean_read";
    r.tid = 0;
    r.pos = 1000;
    r.end_pos = 3000;
    r.cigar_ops = {{'M', 2000}};  // 纯匹配
    r.has_sa = false;
    reads.push_back(r);

    auto components = builder.build(reads);

    std::cout << "  Components: " << components.size() << std::endl;
    std::cout << "  Anchors: " << builder.get_total_anchors() << std::endl;

    // 无信号则无 anchor
    assert(builder.get_total_anchors() == 0);
    assert(components.empty());

    std::cout << "  No signals tests passed!" << std::endl;
}

int main() {
    std::cout << "=== PLACER Phase 3 Tests ===" << std::endl;

    try {
        test_component_builder_empty();
        test_component_builder_no_signals();
        test_component_builder_basic();
        test_component_builder_multi_chrom();
        test_component_builder_span_threshold();
        test_component_builder_ins_event();
        test_component_builder_sa_split();
        test_component_builder_real_bam();

        std::cout << "\n=== All Phase 3 tests passed! ===" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Test failed: " << e.what() << std::endl;
        return 1;
    }
}
