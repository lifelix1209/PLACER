#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <atomic>
#include <cstdio>
#include <fstream>
#include <functional>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <thread>
#include <unistd.h>
#include <vector>

#include "pipeline.h"

namespace {

std::string make_temp_path(const std::string& stem, const std::string& suffix) {
    const long long pid = static_cast<long long>(::getpid());
    return "/tmp/" + stem + "_" + std::to_string(pid) + suffix;
}

void write_text_file(const std::string& path, const std::string& text) {
    std::ofstream out(path);
    assert(out.is_open());
    out << text;
}

void write_executable_script(const std::string& path, const std::string& text) {
    write_text_file(path, text);
    assert(::chmod(path.c_str(), 0700) == 0);
}

std::string build_sequence(size_t len, uint32_t seed) {
    static const char kBases[] = {'A', 'C', 'G', 'T'};
    uint32_t state = seed;
    std::string out;
    out.reserve(len);
    for (size_t i = 0; i < len; ++i) {
        state = (state * 1664525u) + 1013904223u;
        out.push_back(kBases[(state >> 24) & 3u]);
    }
    return out;
}

void write_te_fasta(
    const std::string& path,
    const std::vector<std::pair<std::string, std::string>>& entries) {
    std::ofstream out(path);
    assert(out.is_open());
    for (const auto& entry : entries) {
        out << ">" << entry.first << "\n";
        out << entry.second << "\n";
    }
}

std::string read_text_file(const std::string& path) {
    std::ifstream in(path);
    assert(in.is_open());
    std::ostringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

}  // namespace

int main() {
    using namespace placer;

    const std::string te_fasta = make_temp_path("placer_blast_te_library", ".fa");
    const std::string makeblastdb = make_temp_path("placer_fake_makeblastdb", ".sh");
    const std::string blastn = make_temp_path("placer_fake_blastn", ".sh");
    const std::string query = build_sequence(120, 17u);
    write_te_fasta(
        te_fasta,
        {
            {"GypsyTop#LTR/Gypsy", query},
            {"CopiaRunner#LTR/Copia", build_sequence(120, 97u)},
        });

    write_executable_script(
        makeblastdb,
        "#!/bin/sh\n"
        "out=\"\"\n"
        "while [ \"$#\" -gt 0 ]; do\n"
        "  if [ \"$1\" = \"-out\" ]; then shift; out=\"$1\"; fi\n"
        "  shift\n"
        "done\n"
        "printf 'fake-db\\n' > \"${out}.nhr\"\n"
        "printf 'fake-db\\n' > \"${out}.nin\"\n"
        "printf 'fake-db\\n' > \"${out}.nsq\"\n"
        "exit 0\n");
    write_executable_script(
        blastn,
        "#!/bin/sh\n"
        "query=\"\"\n"
        "out=\"\"\n"
        "db=\"\"\n"
        "while [ \"$#\" -gt 0 ]; do\n"
        "  if [ \"$1\" = \"-query\" ]; then shift; query=\"$1\"; fi\n"
        "  if [ \"$1\" = \"-out\" ]; then shift; out=\"$1\"; fi\n"
        "  if [ \"$1\" = \"-db\" ]; then shift; db=\"$1\"; fi\n"
        "  shift\n"
        "done\n"
        "id=$(grep '^>' \"$query\" | sed 's/^>//' | head -n 1)\n"
        "if [ -f \"${db}.short\" ]; then\n"
        "  printf '%s\\tShortSine#SINE/Short\\t100.000\\t60\\t60\\t1\\t60\\t1\\t60\\t120.0\\t1e-30\\n' \"$id\" > \"$out\"\n"
        "  exit 0\n"
        "fi\n"
        "printf '%s\\tGypsyTop#LTR/Gypsy\\t99.000\\t110\\t120\\t6\\t115\\t1\\t110\\t240.0\\t1e-60\\n' \"$id\" > \"$out\"\n"
        "printf '%s\\tCopiaRunner#LTR/Copia\\t88.000\\t84\\t120\\t19\\t102\\t3\\t86\\t160.0\\t1e-30\\n' \"$id\" >> \"$out\"\n");

    PipelineConfig cfg;
    cfg.te_fasta_path = te_fasta;
    cfg.te_blastn_path = blastn;
    cfg.te_makeblastdb_path = makeblastdb;
    TEKmerQuickClassifierModule classifier(cfg);

    const TEAlignmentEvidence evidence = classifier.align_insert_sequence(query);
    assert(evidence.pass);
    assert(evidence.qc_reason == "PASS_INSERT_TE_ALIGNMENT");
    assert(evidence.best_family == "Gypsy");
    assert(evidence.best_subfamily == "GypsyTop");
    assert(evidence.best_identity > 0.98);
    assert(evidence.best_query_coverage > 0.91);
    assert(evidence.best_query_coverage < 0.92);
    assert(evidence.best_score > evidence.second_score);
    assert(evidence.second_family == "Copia");
    assert(evidence.cross_family_margin > 0.25);
    assert(evidence.sequence_model_label == "TE_MODEL_IN_DISTRIBUTION");
    assert(evidence.sequence_model_score == evidence.best_score);
    assert(evidence.sequence_model_entropy == 0.0);
    assert(evidence.annotation_intervals.find("q=5-115") != std::string::npos);

    const std::string short_te_fasta = make_temp_path("placer_short_blast_te_library", ".fa");
    write_te_fasta(
        short_te_fasta,
        {
            {"ShortSine#SINE/Short", std::string(60, 'A')},
        });
    write_executable_script(
        makeblastdb,
        "#!/bin/sh\n"
        "out=\"\"\n"
        "in=\"\"\n"
        "while [ \"$#\" -gt 0 ]; do\n"
        "  if [ \"$1\" = \"-out\" ]; then shift; out=\"$1\"; fi\n"
        "  if [ \"$1\" = \"-in\" ]; then shift; in=\"$1\"; fi\n"
        "  shift\n"
        "done\n"
        "printf 'fake-db\\n' > \"${out}.nhr\"\n"
        "printf 'fake-db\\n' > \"${out}.nin\"\n"
        "printf 'fake-db\\n' > \"${out}.nsq\"\n"
        "case \"$in\" in *short_blast_te_library*) printf 'short\\n' > \"${out}.short\" ;; esac\n"
        "exit 0\n");

    PipelineConfig short_cfg;
    short_cfg.te_fasta_path = short_te_fasta;
    short_cfg.te_blastn_path = blastn;
    short_cfg.te_makeblastdb_path = makeblastdb;
    TEKmerQuickClassifierModule short_classifier(short_cfg);

    const TEAlignmentEvidence short_evidence =
        short_classifier.align_insert_sequence(std::string(60, 'A'));
    assert(short_evidence.pass);
    assert(short_evidence.best_family == "Short");
    assert(short_evidence.best_subfamily == "ShortSine");

    const std::string batch_te_fasta = make_temp_path("placer_batch_blast_te_library", ".fa");
    const std::string batch_makeblastdb = make_temp_path("placer_batch_fake_makeblastdb", ".sh");
    const std::string batch_blastn = make_temp_path("placer_batch_fake_blastn", ".sh");
    const std::string batch_counter = make_temp_path("placer_batch_blast_counter", ".txt");
    const std::string query_a = build_sequence(120, 301u);
    const std::string query_b = build_sequence(120, 709u);
    write_te_fasta(
        batch_te_fasta,
        {
            {"GypsyBatch#LTR/Gypsy", query_a},
            {"CopiaBatch#LTR/Copia", query_b},
        });
    write_executable_script(
        batch_makeblastdb,
        "#!/bin/sh\n"
        "out=\"\"\n"
        "while [ \"$#\" -gt 0 ]; do\n"
        "  if [ \"$1\" = \"-out\" ]; then shift; out=\"$1\"; fi\n"
        "  shift\n"
        "done\n"
        "printf 'fake-db\\n' > \"${out}.nhr\"\n"
        "printf 'fake-db\\n' > \"${out}.nin\"\n"
        "printf 'fake-db\\n' > \"${out}.nsq\"\n"
        "exit 0\n");
    write_executable_script(
        batch_blastn,
        "#!/bin/sh\n"
        "query=\"\"\n"
        "out=\"\"\n"
        "counter='" + batch_counter + "'\n"
        "while [ \"$#\" -gt 0 ]; do\n"
        "  if [ \"$1\" = \"-query\" ]; then shift; query=\"$1\"; fi\n"
        "  if [ \"$1\" = \"-out\" ]; then shift; out=\"$1\"; fi\n"
        "  shift\n"
        "done\n"
        "count=0\n"
        "if [ -f \"$counter\" ]; then count=$(cat \"$counter\"); fi\n"
        "count=$((count + 1))\n"
        "printf '%s\\n' \"$count\" > \"$counter\"\n"
        "ids=$(grep '^>' \"$query\" | sed 's/^>//')\n"
        "record_count=$(printf '%s\\n' \"$ids\" | sed '/^$/d' | wc -l | tr -d ' ')\n"
        "if [ \"$record_count\" -ne 2 ]; then exit 17; fi\n"
        ": > \"$out\"\n"
        "idx=0\n"
        "for id in $ids; do\n"
        "  idx=$((idx + 1))\n"
        "  if [ \"$idx\" -eq 1 ]; then\n"
        "    printf '%s\\tGypsyBatch#LTR/Gypsy\\t99.000\\t110\\t120\\t6\\t115\\t1\\t110\\t240.0\\t1e-60\\n' \"$id\" >> \"$out\"\n"
        "  else\n"
        "    printf '%s\\tCopiaBatch#LTR/Copia\\t98.000\\t108\\t120\\t7\\t114\\t2\\t109\\t230.0\\t1e-55\\n' \"$id\" >> \"$out\"\n"
        "  fi\n"
        "done\n"
        "exit 0\n");

    PipelineConfig batch_cfg;
    batch_cfg.te_fasta_path = batch_te_fasta;
    batch_cfg.te_blastn_path = batch_blastn;
    batch_cfg.te_makeblastdb_path = batch_makeblastdb;
    TEKmerQuickClassifierModule batch_classifier(batch_cfg);

    const std::vector<TEAlignmentEvidence> batch_evidence =
        batch_classifier.align_insert_sequences({query_a, query_b, query_a});
    assert(batch_evidence.size() == 3);
    assert(batch_evidence[0].pass);
    assert(batch_evidence[1].pass);
    assert(batch_evidence[2].pass);
    assert(batch_evidence[0].best_family == "Gypsy");
    assert(batch_evidence[1].best_family == "Copia");
    assert(batch_evidence[2].best_family == "Gypsy");
    assert(batch_evidence[0].best_score == batch_evidence[2].best_score);
    assert(read_text_file(batch_counter) == "1\n");
    const std::vector<TEAlignmentEvidence> cached_batch_evidence =
        batch_classifier.align_insert_sequences({query_a});
    assert(cached_batch_evidence.size() == 1);
    assert(cached_batch_evidence[0].best_family == "Gypsy");
    assert(read_text_file(batch_counter) == "1\n");

    const std::string inflight_te_fasta =
        make_temp_path("placer_inflight_blast_te_library", ".fa");
    const std::string inflight_makeblastdb =
        make_temp_path("placer_inflight_fake_makeblastdb", ".sh");
    const std::string inflight_blastn =
        make_temp_path("placer_inflight_fake_blastn", ".sh");
    const std::string inflight_counter =
        make_temp_path("placer_inflight_blast_counter", ".txt");
    const std::string query_c = build_sequence(128, 1709u);
    write_te_fasta(
        inflight_te_fasta,
        {
            {"InFlightGypsy#LTR/Gypsy", query_c},
        });
    write_executable_script(
        inflight_makeblastdb,
        "#!/bin/sh\n"
        "out=\"\"\n"
        "while [ \"$#\" -gt 0 ]; do\n"
        "  if [ \"$1\" = \"-out\" ]; then shift; out=\"$1\"; fi\n"
        "  shift\n"
        "done\n"
        "printf 'fake-db\\n' > \"${out}.nhr\"\n"
        "printf 'fake-db\\n' > \"${out}.nin\"\n"
        "printf 'fake-db\\n' > \"${out}.nsq\"\n"
        "exit 0\n");
    write_executable_script(
        inflight_blastn,
        "#!/bin/sh\n"
        "query=\"\"\n"
        "out=\"\"\n"
        "counter='" + inflight_counter + "'\n"
        "while [ \"$#\" -gt 0 ]; do\n"
        "  if [ \"$1\" = \"-query\" ]; then shift; query=\"$1\"; fi\n"
        "  if [ \"$1\" = \"-out\" ]; then shift; out=\"$1\"; fi\n"
        "  shift\n"
        "done\n"
        "printf 'x\\n' >> \"$counter\"\n"
        "sleep 1\n"
        "id=$(grep '^>' \"$query\" | sed 's/^>//' | head -n 1)\n"
        "printf '%s\\tInFlightGypsy#LTR/Gypsy\\t99.000\\t120\\t128\\t5\\t124\\t1\\t120\\t250.0\\t1e-65\\n' \"$id\" > \"$out\"\n"
        "exit 0\n");

    PipelineConfig inflight_cfg;
    inflight_cfg.te_fasta_path = inflight_te_fasta;
    inflight_cfg.te_blastn_path = inflight_blastn;
    inflight_cfg.te_makeblastdb_path = inflight_makeblastdb;
    TEKmerQuickClassifierModule inflight_classifier_a(inflight_cfg);
    TEKmerQuickClassifierModule inflight_classifier_b(inflight_cfg);

    std::atomic<int> ready_threads{0};
    std::atomic<bool> start_threads{false};
    TEAlignmentEvidence inflight_a;
    TEAlignmentEvidence inflight_b;
    auto align_inflight = [&](TEKmerQuickClassifierModule& classifier,
                              TEAlignmentEvidence& evidence) {
        ready_threads.fetch_add(1);
        while (!start_threads.load()) {
            std::this_thread::yield();
        }
        evidence = classifier.align_insert_sequence(query_c);
    };
    std::thread thread_a(
        align_inflight,
        std::ref(inflight_classifier_a),
        std::ref(inflight_a));
    std::thread thread_b(
        align_inflight,
        std::ref(inflight_classifier_b),
        std::ref(inflight_b));
    while (ready_threads.load() < 2) {
        std::this_thread::yield();
    }
    start_threads.store(true);
    thread_a.join();
    thread_b.join();
    assert(inflight_a.pass);
    assert(inflight_b.pass);
    assert(inflight_a.best_family == "Gypsy");
    assert(inflight_b.best_family == "Gypsy");
    assert(read_text_file(inflight_counter) == "x\n");

    std::remove(te_fasta.c_str());
    std::remove(short_te_fasta.c_str());
    std::remove(makeblastdb.c_str());
    std::remove(blastn.c_str());
    std::remove(batch_te_fasta.c_str());
    std::remove(batch_makeblastdb.c_str());
    std::remove(batch_blastn.c_str());
    std::remove(batch_counter.c_str());
    std::remove(inflight_te_fasta.c_str());
    std::remove(inflight_makeblastdb.c_str());
    std::remove(inflight_blastn.c_str());
    std::remove(inflight_counter.c_str());
    return 0;
}
