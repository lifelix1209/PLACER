#include "pipeline.h"
#include "gate1_module.h"

#include <algorithm>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <thread>
#include <utility>

#include <htslib/sam.h>

namespace placer {
namespace {

int32_t classify_tier(double delta_score) {
    if (delta_score >= 30.0) {
        return 1;
    }
    if (delta_score >= 10.0) {
        return 2;
    }
    return 3;
}

template <typename T>
class SafeQueue {
public:
    void push(T&& value) {
        {
            std::lock_guard<std::mutex> lock(mu_);
            queue_.push(std::move(value));
        }
        cv_.notify_one();
    }

    bool pop(T& out) {
        std::unique_lock<std::mutex> lock(mu_);
        cv_.wait(lock, [this]() { return finished_ || !queue_.empty(); });
        if (queue_.empty()) {
            return false;
        }
        out = std::move(queue_.front());
        queue_.pop();
        return true;
    }

    void close() {
        {
            std::lock_guard<std::mutex> lock(mu_);
            finished_ = true;
        }
        cv_.notify_all();
    }

private:
    std::queue<T> queue_;
    std::mutex mu_;
    std::condition_variable cv_;
    bool finished_ = false;
};

}  // namespace

std::vector<ComponentCall> LinearBinComponentModule::build(
    const std::vector<BamRecordPtr>& bin_records,
    const std::string& chrom,
    int32_t tid,
    int32_t bin_start,
    int32_t bin_end) const {
    if (bin_records.empty()) {
        return {};
    }

    ComponentCall call;
    call.chrom = chrom;
    call.tid = tid;
    call.bin_start = bin_start;
    call.bin_end = bin_end;
    call.anchor_pos = (bin_start + bin_end) / 2;

    call.read_indices.reserve(bin_records.size());
    for (size_t i = 0; i < bin_records.size(); ++i) {
        call.read_indices.push_back(i);

        if (!bin_records[i]) {
            continue;
        }

        ReadView view(bin_records[i].get());
        const uint16_t flag = view.flag();

        if (view.has_sa_tag() || ((flag & BAM_FSUPPLEMENTARY) != 0)) {
            call.split_sa_read_indices.push_back(i);
        }

        const uint32_t* cigar = view.cigar();
        const int32_t n_cigar = view.n_cigar();
        int32_t max_soft = 0;
        int32_t max_ins = 0;
        for (int32_t ci = 0; ci < n_cigar; ++ci) {
            const int op = bam_cigar_op(cigar[ci]);
            const int32_t len = static_cast<int32_t>(bam_cigar_oplen(cigar[ci]));
            if (op == BAM_CSOFT_CLIP) {
                max_soft = std::max(max_soft, len);
            } else if (op == BAM_CINS) {
                max_ins = std::max(max_ins, len);
            }
        }

        if (max_soft >= 20) {
            call.soft_clip_read_indices.push_back(i);
        }
        if (max_ins >= 50) {
            call.insertion_read_indices.push_back(i);
        }
    }

    return {std::move(call)};
}

std::vector<LocusEvidence> SimpleLocalRealignModule::collect(
    const ComponentCall& component,
    const std::vector<BamRecordPtr>& bin_records) const {
    std::vector<LocusEvidence> evidence;
    evidence.reserve(component.read_indices.size());

    const int32_t center = (component.bin_start + component.bin_end) / 2;

    for (size_t idx : component.read_indices) {
        if (idx >= bin_records.size() || !bin_records[idx]) {
            continue;
        }

        ReadView view(bin_records[idx].get());
        LocusEvidence e;
        e.tid = view.tid();
        e.pos = view.pos();
        e.read_index = idx;

        const int32_t distance = std::abs(view.pos() - center);
        e.normalized_score = std::max(0.0, 1.0 - static_cast<double>(distance) / 10000.0);
        evidence.push_back(e);
    }

    return evidence;
}

AssemblyCall LazyDecodeAssemblyModule::assemble(
    const ComponentCall& component,
    const std::vector<LocusEvidence>& evidence,
    const std::vector<BamRecordPtr>& bin_records) const {
    (void)evidence;

    AssemblyCall call;
    call.tid = component.tid;
    call.pos = component.anchor_pos;

    if (!component.read_indices.empty()) {
        const size_t idx = component.read_indices.front();
        if (idx < bin_records.size() && bin_records[idx]) {
            ReadView view(bin_records[idx].get());
            std::string seq = view.decode_sequence();  // Lazy decode: only here.
            if (seq.size() > 300) {
                seq.resize(300);
            }
            call.consensus = std::move(seq);
        }
    }

    return call;
}

PlaceabilityReport SimplePlaceabilityModule::score(
    const AssemblyCall& assembly,
    const std::vector<LocusEvidence>& evidence) const {
    PlaceabilityReport report;
    report.tid = assembly.tid;
    report.pos = assembly.pos;
    report.support_reads = static_cast<int32_t>(evidence.size());

    double score_sum = 0.0;
    for (const auto& row : evidence) {
        score_sum += row.normalized_score;
    }

    report.delta_score = report.support_reads > 0
        ? (score_sum / static_cast<double>(report.support_reads)) * 100.0
        : 0.0;
    report.tier = classify_tier(report.delta_score);

    return report;
}

GenotypeCall SimpleGenotypingModule::genotype(
    const AssemblyCall& assembly,
    const PlaceabilityReport& placeability,
    const std::vector<LocusEvidence>& evidence) const {
    (void)evidence;

    GenotypeCall call;
    call.tid = assembly.tid;
    call.pos = assembly.pos;

    const double af = std::clamp(static_cast<double>(placeability.support_reads) / 20.0, 0.0, 1.0);
    call.af = af;
    call.gq = static_cast<int32_t>(std::round(placeability.delta_score));

    if (af >= 0.8) {
        call.genotype = "1/1";
    } else if (af >= 0.2) {
        call.genotype = "0/1";
    } else {
        call.genotype = "0/0";
    }

    return call;
}

Pipeline::Pipeline(
    PipelineConfig config,
    std::unique_ptr<BamStreamReader> bam_reader,
    std::unique_ptr<IGate1Module> gate1_module,
    std::unique_ptr<IComponentModule> component_module,
    std::unique_ptr<ILocalRealignModule> local_realign_module,
    std::unique_ptr<IAssemblyModule> assembly_module,
    std::unique_ptr<IPlaceabilityModule> placeability_module,
    std::unique_ptr<IGenotypingModule> genotyping_module,
    std::unique_ptr<IInsertionFragmentModule> ins_fragment_module,
    std::unique_ptr<ITEQuickClassifierModule> te_classifier_module)
    : config_(std::move(config)),
      bam_reader_(std::move(bam_reader)),
      gate1_module_(std::move(gate1_module)),
      component_module_(std::move(component_module)),
      local_realign_module_(std::move(local_realign_module)),
      assembly_module_(std::move(assembly_module)),
      placeability_module_(std::move(placeability_module)),
      genotyping_module_(std::move(genotyping_module)),
      ins_fragment_module_(std::move(ins_fragment_module)),
      te_classifier_module_(std::move(te_classifier_module)) {}

PipelineResult Pipeline::run() const {
    if (!bam_reader_ || !bam_reader_->is_valid()) {
        throw std::runtime_error("BAM reader is not valid");
    }
    if (!gate1_module_ || !component_module_ || !local_realign_module_ ||
        !assembly_module_ || !placeability_module_ || !genotyping_module_ ||
        !ins_fragment_module_ || !te_classifier_module_) {
        throw std::runtime_error("Pipeline modules are not fully configured");
    }

    if (config_.enable_parallel) {
        return run_parallel();
    }
    return run_streaming();
}

PipelineResult Pipeline::run_streaming() const {
    PipelineResult result;
    StreamingState state;

    const auto progress_cb = [this](int64_t processed, int32_t tid) {
        std::cerr << "[Pipeline] processed=" << processed << " current_tid=" << tid << '\n';
        return true;
    };

    result.total_reads = bam_reader_->stream(
        [this, &state, &result](BamRecordPtr&& record) {
            consume_record(std::move(record), state, result);
        },
        progress_cb,
        config_.progress_interval);

    flush_current_bin(state, result);
    return result;
}

PipelineResult Pipeline::run_parallel() const {
    PipelineResult result;
    SafeQueue<std::vector<BamRecordPtr>> batch_queue;

    PipelineResult worker_result;
    std::thread worker([this, &batch_queue, &worker_result]() {
        StreamingState worker_state;
        std::vector<BamRecordPtr> batch;
        while (batch_queue.pop(batch)) {
            for (auto& record : batch) {
                consume_record(std::move(record), worker_state, worker_result);
            }
            batch.clear();
        }
        flush_current_bin(worker_state, worker_result);
    });

    std::vector<BamRecordPtr> current_batch;
    current_batch.reserve(config_.batch_size);

    const auto progress_cb = [this](int64_t processed, int32_t tid) {
        std::cerr << "[Pipeline] processed=" << processed << " current_tid=" << tid << '\n';
        return true;
    };

    const int64_t total_reads = bam_reader_->stream(
        [this, &current_batch, &batch_queue](BamRecordPtr&& record) {
            current_batch.push_back(std::move(record));
            if (current_batch.size() >= config_.batch_size) {
                batch_queue.push(std::move(current_batch));
                current_batch.clear();
                current_batch.reserve(config_.batch_size);
            }
        },
        progress_cb,
        config_.progress_interval);

    if (!current_batch.empty()) {
        batch_queue.push(std::move(current_batch));
    }

    batch_queue.close();
    worker.join();

    result = std::move(worker_result);
    result.total_reads = total_reads;
    return result;
}

void Pipeline::consume_record(
    BamRecordPtr&& record,
    StreamingState& state,
    PipelineResult& result) const {
    if (!record) {
        return;
    }

    ReadView view(record.get());
    if (!gate1_module_->pass_preliminary(view)) {
        return;
    }

    result.gate1_passed += 1;

    while (!state.active_window.empty()) {
        const auto& front = state.active_window.front();
        const bool different_tid = (view.tid() != front.tid);
        const bool beyond_window = (!different_tid && (view.pos() - front.pos > config_.window_size));
        if (!different_tid && !beyond_window) {
            break;
        }
        state.active_window.pop_front();
    }
    state.active_window.push_back({view.tid(), view.pos()});

    const int32_t bin_index = (config_.bin_size > 0) ? (view.pos() / config_.bin_size) : 0;

    if (state.current_tid < 0) {
        state.current_tid = view.tid();
        state.current_bin_index = bin_index;
    }

    if (view.tid() != state.current_tid || bin_index != state.current_bin_index) {
        flush_current_bin(state, result);
        state.current_tid = view.tid();
        state.current_bin_index = bin_index;
    }

    state.current_bin_records.push_back(std::move(record));
}

void Pipeline::flush_current_bin(
    StreamingState& state,
    PipelineResult& result) const {
    if (state.current_bin_records.empty()) {
        return;
    }

    process_bin_records(
        std::move(state.current_bin_records),
        state.current_tid,
        state.current_bin_index,
        result);

    state.current_bin_records.clear();
}

void Pipeline::process_bin_records(
    std::vector<BamRecordPtr>&& bin_records,
    int32_t tid,
    int32_t bin_index,
    PipelineResult& result) const {
    if (bin_records.empty()) {
        return;
    }

    result.processed_bins += 1;

    const int32_t bin_start = bin_index * config_.bin_size;
    const int32_t bin_end = bin_start + config_.bin_size;
    const std::string chrom = bam_reader_->chromosome_name(tid);

    auto components = component_module_->build(bin_records, chrom, tid, bin_start, bin_end);
    result.built_components += static_cast<int64_t>(components.size());

    for (const auto& component : components) {
        ClusterTECall te_call;
        if (ins_fragment_module_ && te_classifier_module_) {
            auto fragments = ins_fragment_module_->extract(component, bin_records);
            if (te_classifier_module_->is_enabled()) {
                auto hits = te_classifier_module_->classify(fragments);
                te_call = te_classifier_module_->vote_cluster(hits);
            }
        }

        auto evidence = local_realign_module_->collect(component, bin_records);
        result.evidence_rows += static_cast<int64_t>(evidence.size());

        AssemblyCall assembly = assembly_module_->assemble(component, evidence, bin_records);
        result.assembled_calls += 1;

        PlaceabilityReport placeability = placeability_module_->score(assembly, evidence);
        result.placeability_calls += 1;

        GenotypeCall genotype = genotyping_module_->genotype(assembly, placeability, evidence);
        result.genotype_calls += 1;

        FinalCall call;
        call.chrom = component.chrom;
        call.tid = assembly.tid;
        call.pos = assembly.pos;
        call.window_start = component.bin_start;
        call.window_end = component.bin_end;

        if (te_call.passed) {
            call.te_name = te_call.te_name;
        }
        call.te_vote_fraction = te_call.vote_fraction;
        call.te_median_identity = te_call.median_identity;
        call.te_fragment_count = te_call.fragment_count;

        call.tier = placeability.tier;
        call.support_reads = placeability.support_reads;
        call.genotype = genotype.genotype;
        call.af = genotype.af;
        call.gq = genotype.gq;
        result.final_calls.push_back(std::move(call));
    }
}

PipelineBuilder::PipelineBuilder(PipelineConfig config)
    : config_(std::move(config)) {}

PipelineBuilder& PipelineBuilder::with_bam_reader(std::unique_ptr<BamStreamReader> reader) {
    bam_reader_ = std::move(reader);
    return *this;
}

PipelineBuilder& PipelineBuilder::with_gate1_module(std::unique_ptr<IGate1Module> module) {
    gate1_module_ = std::move(module);
    return *this;
}

PipelineBuilder& PipelineBuilder::with_component_module(std::unique_ptr<IComponentModule> module) {
    component_module_ = std::move(module);
    return *this;
}

PipelineBuilder& PipelineBuilder::with_local_realign_module(std::unique_ptr<ILocalRealignModule> module) {
    local_realign_module_ = std::move(module);
    return *this;
}

PipelineBuilder& PipelineBuilder::with_assembly_module(std::unique_ptr<IAssemblyModule> module) {
    assembly_module_ = std::move(module);
    return *this;
}

PipelineBuilder& PipelineBuilder::with_placeability_module(std::unique_ptr<IPlaceabilityModule> module) {
    placeability_module_ = std::move(module);
    return *this;
}

PipelineBuilder& PipelineBuilder::with_genotyping_module(std::unique_ptr<IGenotypingModule> module) {
    genotyping_module_ = std::move(module);
    return *this;
}

PipelineBuilder& PipelineBuilder::with_insertion_fragment_module(std::unique_ptr<IInsertionFragmentModule> module) {
    ins_fragment_module_ = std::move(module);
    return *this;
}

PipelineBuilder& PipelineBuilder::with_te_classifier_module(std::unique_ptr<ITEQuickClassifierModule> module) {
    te_classifier_module_ = std::move(module);
    return *this;
}

std::unique_ptr<Pipeline> PipelineBuilder::build() {
    if (!bam_reader_) {
        bam_reader_ = make_bam_reader(config_.bam_path, config_.bam_threads);
    }
    if (!gate1_module_) {
        gate1_module_ = std::make_unique<SignalFirstGate1Module>();
    }
    if (!component_module_) {
        component_module_ = std::make_unique<LinearBinComponentModule>();
    }
    if (!local_realign_module_) {
        local_realign_module_ = std::make_unique<SimpleLocalRealignModule>();
    }
    if (!assembly_module_) {
        assembly_module_ = std::make_unique<LazyDecodeAssemblyModule>();
    }
    if (!placeability_module_) {
        placeability_module_ = std::make_unique<SimplePlaceabilityModule>();
    }
    if (!genotyping_module_) {
        genotyping_module_ = std::make_unique<SimpleGenotypingModule>();
    }
    if (!ins_fragment_module_) {
        ins_fragment_module_ = std::make_unique<CigarInsertionFragmentModule>(config_);
    }
    if (!te_classifier_module_) {
        te_classifier_module_ = std::make_unique<TEKmerQuickClassifierModule>(config_);
    }

    return std::make_unique<Pipeline>(
        config_,
        std::move(bam_reader_),
        std::move(gate1_module_),
        std::move(component_module_),
        std::move(local_realign_module_),
        std::move(assembly_module_),
        std::move(placeability_module_),
        std::move(genotyping_module_),
        std::move(ins_fragment_module_),
        std::move(te_classifier_module_));
}

std::unique_ptr<Pipeline> build_default_pipeline(const PipelineConfig& config) {
    PipelineBuilder builder(config);
    return builder.build();
}

}  // namespace placer
