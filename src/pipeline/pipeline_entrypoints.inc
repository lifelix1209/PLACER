std::vector<ComponentCall> LinearBinComponentModule::build(
    const std::vector<const bam1_t*>& bin_records,
    const std::string& chrom,
    int32_t tid,
    int32_t bin_start,
    int32_t bin_end) const {
    (void)bin_start;
    (void)bin_end;
    DbscanComponentModule module;
    return module.build(bin_records, chrom, tid);
}

Pipeline::Pipeline(PipelineConfig config, std::unique_ptr<BamStreamReader> bam_reader)
    : config_(std::move(config)),
      bam_reader_(std::move(bam_reader)),
      gate1_module_(),
      component_module_(),
      ins_fragment_module_(config_),
      te_classifier_module_(config_),
      tsd_detector_(config_) {}

PipelineResult Pipeline::run() const {
    if (!bam_reader_ || !bam_reader_->is_valid()) {
        throw std::runtime_error("BAM reader is not valid");
    }
    if (!bam_reader_->can_fetch()) {
        throw std::runtime_error(
            "BAM reader does not support indexed local fetch; a BAM index is required for event-level recollection");
    }
    if (!tsd_detector_.can_fetch_reference()) {
        throw std::runtime_error(
            "Reference FASTA cannot be fetched; event consensus segmentation requires a readable indexed reference");
    }
    if (!te_classifier_module_.is_enabled()) {
        throw std::runtime_error(
            "TE library alignment is unavailable; Stage 5 requires a readable TE FASTA and a usable shortlist index");
    }

    {
        std::ostringstream oss;
        oss << "[Pipeline] starting mode="
            << (config_.enable_parallel ? "parallel" : "streaming")
            << " bam_threads=" << config_.bam_threads
            << " progress_interval=" << config_.progress_interval
            << " log_stage_bins=" << bool_as_int(config_.log_stage_bins)
            << " log_stage_components=" << bool_as_int(config_.log_stage_components);
        emit_pipeline_log_line(oss.str());
    }

    if (config_.enable_parallel) {
        return run_parallel();
    }
    return run_streaming();
}

PipelineResult Pipeline::run_streaming() const {
    const auto pipeline_started = std::chrono::steady_clock::now();
    PipelineResult result;
    const int32_t bin_size = std::max(1, config_.bin_size);
    ExactBinScanState state;

    const auto process_task = [&](const ExactBinTask& task) {
        std::vector<const bam1_t*> bin_records = collect_task_record_ptrs(task);
        if (bin_records.empty()) {
            return;
        }
        process_bin_records(
            std::move(bin_records),
            task.tid,
            task.bin_index,
            task.bin_start,
            task.bin_end,
            result);
    };

    const auto process_ready = [&](std::vector<ExactBinTask> ready) {
        for (const auto& task : ready) {
            process_task(task);
        }
    };

    const auto progress_cb = [this](int64_t processed, int32_t tid) {
        std::cerr << "[Pipeline] processed=" << processed << " current_tid=" << tid << '\n';
        return true;
    };

    const auto stream_started = std::chrono::steady_clock::now();
    result.total_reads = bam_reader_->stream(
        [this, &state, &result, bin_size, &process_ready](BamRecordPtr&& record) {
            if (!record) {
                return;
            }

            ReadView view(record.get());
            const auto gate1_started = std::chrono::steady_clock::now();
            const bool gate1_pass = gate1_module_.pass_preliminary(view);
            result.performance.gate1_preliminary_seconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - gate1_started).count();
            if (!gate1_pass) {
                return;
            }

            result.gate1_passed += 1;
            const int32_t tid = view.tid();
            const int32_t bin_index = exact_bin_index_for_record(record.get(), bin_size);

            BufferedRecord buffered;
            buffered.ref_end = compute_ref_end(view);
            buffered.record = std::move(record);

            const auto scan_started = std::chrono::steady_clock::now();
            std::vector<ExactBinTask> ready = append_record_to_exact_bin_scan(
                state,
                std::move(buffered),
                tid,
                bin_index,
                bin_size,
                kCrossBinContextBins,
                kWindowBinSlackBp);
            result.performance.exact_bin_scan_seconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - scan_started).count();
            process_ready(std::move(ready));
        },
        progress_cb,
        config_.progress_interval);
    result.performance.bam_stream_wall_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - stream_started).count();

    const auto flush_started = std::chrono::steady_clock::now();
    std::vector<ExactBinTask> flush_ready = flush_exact_bin_scan(
        state,
        bin_size,
        kCrossBinContextBins,
        kWindowBinSlackBp);
    result.performance.exact_bin_scan_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - flush_started).count();
    process_ready(std::move(flush_ready));
    FinalCallFilterConfig final_filter_config;
    final_filter_config.min_raw_cigar_insert_len_bp =
        config_.min_final_raw_cigar_insert_len_bp;
    final_filter_config.report_mode = config_.final_report_mode;
    const auto finalization_started = std::chrono::steady_clock::now();
    finalize_final_calls(
        result,
        tsd_detector_,
        config_.final_fdr_q,
        final_filter_config);
    result.performance.finalization_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - finalization_started).count();
    result.performance.pipeline_wall_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - pipeline_started).count();
    return result;
}

PipelineResult Pipeline::run_parallel() const {
    const auto pipeline_started = std::chrono::steady_clock::now();
    PipelineResult result;
    const int32_t bin_size = std::max(1, config_.bin_size);
    int32_t worker_count = config_.parallel_workers;
    if (worker_count <= 0) {
        worker_count = static_cast<int32_t>(std::thread::hardware_concurrency());
    }
    worker_count = std::max(1, worker_count);

    const size_t ready_capacity = static_cast<size_t>(
        config_.parallel_queue_max_tasks > 0
            ? config_.parallel_queue_max_tasks
            : std::max(4, worker_count * 4));
    const size_t result_capacity = static_cast<size_t>(
        config_.parallel_result_buffer_max > 0
            ? config_.parallel_result_buffer_max
            : std::max(4, worker_count * 4));

    ExactBinScanState scan_state;
    CostAwareTaskQueue ready_queue(ready_capacity);
    BoundedTaskQueue<ExactBinTaskResult> result_queue(result_capacity);
    DeterministicBinReducer reducer;
    std::mutex reducer_mutex;
    std::mutex error_mutex;
    std::exception_ptr first_error;
    std::atomic<int64_t> tasks_materialized{0};
    std::atomic<int64_t> tasks_running{0};
    std::atomic<int64_t> tasks_completed{0};
    std::atomic<int64_t> workers_waiting{0};
    std::atomic<double> max_task_elapsed_s{0.0};
    std::atomic<int64_t> reducer_committed{0};
    std::atomic<int64_t> raw_bins_discovered{0};
    std::atomic<size_t> snapshot_window_size{0};
    double gate1_preliminary_seconds = 0.0;
    double exact_bin_scan_seconds = 0.0;

    const auto update_max_elapsed = [&](double elapsed_s) {
        double cur = max_task_elapsed_s.load();
        while (elapsed_s > cur &&
               !max_task_elapsed_s.compare_exchange_weak(cur, elapsed_s)) {
        }
    };

    const auto set_error = [&](std::exception_ptr error) {
        {
            std::lock_guard<std::mutex> lock(error_mutex);
            if (!first_error) {
                first_error = error;
            }
        }
        ready_queue.close();
        result_queue.close();
    };

    const auto rethrow_if_error = [&]() {
        std::exception_ptr error;
        {
            std::lock_guard<std::mutex> lock(error_mutex);
            error = first_error;
        }
        if (error) {
            std::rethrow_exception(error);
        }
    };

    std::thread reducer_thread([&]() {
        try {
            ExactBinTaskResult task_result;
            while (result_queue.pop(task_result)) {
                std::lock_guard<std::mutex> lock(reducer_mutex);
                reducer.push_ready(
                    task_result.tid,
                    task_result.bin_index,
                    std::move(task_result));
                reducer.drain([&](
                                  int32_t tid,
                                  int32_t bin_index,
                                  ExactBinTaskResult payload) {
                    if (payload.tid != tid || payload.bin_index != bin_index) {
                        throw std::runtime_error(
                            "parallel reducer received mismatched task result");
                    }
                    append_partial_pipeline_result(result, std::move(payload.partial));
                    reducer_committed.fetch_add(1);
                });
            }
        } catch (...) {
            set_error(std::current_exception());
        }
    });

    std::vector<std::thread> workers;
    workers.reserve(static_cast<size_t>(worker_count));
    for (int32_t worker_id = 0; worker_id < worker_count; ++worker_id) {
        workers.emplace_back([&, worker_id]() {
            (void)worker_id;
            try {
                PipelineConfig worker_config = config_;
                worker_config.enable_parallel = false;
                worker_config.bam_threads = 1;
                auto worker_reader = bam_reader_->clone(1);
                Pipeline worker_pipeline(worker_config, std::move(worker_reader));

                ExactBinTask task;
                while (true) {
                    workers_waiting.fetch_add(1);
                    const bool popped = ready_queue.pop(task);
                    workers_waiting.fetch_sub(1);
                    if (!popped) {
                        break;
                    }
                    tasks_running.fetch_add(1);
                    PipelineResult partial;
                    const auto started = std::chrono::steady_clock::now();
                    std::vector<const bam1_t*> bin_records =
                        collect_task_record_ptrs(task);
                    if (!bin_records.empty()) {
                        worker_pipeline.process_bin_records(
                            std::move(bin_records),
                            task.tid,
                            task.bin_index,
                            task.bin_start,
                            task.bin_end,
                            partial);
                    }

                    ExactBinTaskResult task_result;
                    task_result.tid = task.tid;
                    task_result.bin_index = task.bin_index;
                    task_result.partial = std::move(partial);
                    task_result.elapsed_s = std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - started).count();
                    update_max_elapsed(task_result.elapsed_s);
                    tasks_running.fetch_sub(1);
                    tasks_completed.fetch_add(1);
                    if (!result_queue.push(std::move(task_result))) {
                        break;
                    }
                }
            } catch (...) {
                tasks_running.store(0);
                set_error(std::current_exception());
            }
        });
    }

    const auto emit_ready_tasks = [&](std::vector<ExactBinTask> ready_tasks) {
        for (auto& task : ready_tasks) {
            {
                std::lock_guard<std::mutex> lock(reducer_mutex);
                reducer.expect(task.tid, task.bin_index);
            }
            raw_bins_discovered.fetch_add(1);
            tasks_materialized.fetch_add(1);
            if (!ready_queue.push(std::move(task))) {
                rethrow_if_error();
                throw std::runtime_error(
                    "parallel ready queue closed before all tasks were submitted");
            }
            rethrow_if_error();
        }
        snapshot_window_size.store(scan_state.recent_snapshots.size());
    };

    const auto progress_cb = [this](int64_t processed, int32_t tid) {
        std::cerr << "[Pipeline] processed=" << processed << " current_tid=" << tid << '\n';
        return true;
    };

    try {
        const auto stream_started = std::chrono::steady_clock::now();
        result.total_reads = bam_reader_->stream(
            [
                this,
                &result,
                &scan_state,
                bin_size,
                &emit_ready_tasks,
                &rethrow_if_error,
                &gate1_preliminary_seconds,
                &exact_bin_scan_seconds
            ](
                BamRecordPtr&& record) {
                rethrow_if_error();
                if (!record) {
                    return;
                }

                ReadView view(record.get());
                const auto gate1_started = std::chrono::steady_clock::now();
                const bool gate1_pass = gate1_module_.pass_preliminary(view);
                gate1_preliminary_seconds += std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - gate1_started).count();
                if (!gate1_pass) {
                    return;
                }

                result.gate1_passed += 1;
                const int32_t tid = view.tid();
                const int32_t bin_index = exact_bin_index_for_record(record.get(), bin_size);

                BufferedRecord buffered;
                buffered.ref_end = compute_ref_end(view);
                buffered.record = std::move(record);
                const auto scan_started = std::chrono::steady_clock::now();
                std::vector<ExactBinTask> ready = append_record_to_exact_bin_scan(
                    scan_state,
                    std::move(buffered),
                    tid,
                    bin_index,
                    bin_size,
                    kCrossBinContextBins,
                    kWindowBinSlackBp);
                exact_bin_scan_seconds += std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - scan_started).count();
                emit_ready_tasks(std::move(ready));
            },
            progress_cb,
            config_.progress_interval);
        result.performance.bam_stream_wall_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - stream_started).count();

        const auto flush_started = std::chrono::steady_clock::now();
        std::vector<ExactBinTask> flush_ready = flush_exact_bin_scan(
            scan_state,
            bin_size,
            kCrossBinContextBins,
            kWindowBinSlackBp);
        exact_bin_scan_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - flush_started).count();
        emit_ready_tasks(std::move(flush_ready));
    } catch (...) {
        set_error(std::current_exception());
    }

    ready_queue.close();
    for (auto& worker : workers) {
        if (worker.joinable()) {
            worker.join();
        }
    }
    result_queue.close();
    if (reducer_thread.joinable()) {
        reducer_thread.join();
    }
    rethrow_if_error();

    {
        std::lock_guard<std::mutex> lock(reducer_mutex);
        if (reducer.committed_count() != reducer.expected_count() ||
            reducer.has_pending_ready()) {
            throw std::runtime_error(
                "parallel reducer finished with undrained task results");
        }
    }

    if (config_.log_parallel_progress) {
        ParallelProgressSnapshot snapshot;
        snapshot.reads_scanned = result.total_reads;
        snapshot.raw_bins_discovered = raw_bins_discovered.load();
        snapshot.tasks_materialized = tasks_materialized.load();
        snapshot.tasks_running = tasks_running.load();
        snapshot.tasks_completed = tasks_completed.load();
        snapshot.worker_count = worker_count;
        snapshot.workers_waiting = workers_waiting.load();
        snapshot.max_task_elapsed_s = max_task_elapsed_s.load();
        snapshot.reducer_committed = reducer_committed.load();
        snapshot.snapshot_window_size = snapshot_window_size.load();
        snapshot.ready_queue_depth = ready_queue.size();
        snapshot.result_queue_depth = result_queue.size();
        emit_pipeline_log_line(render_parallel_progress(snapshot));
    }

    FinalCallFilterConfig final_filter_config;
    final_filter_config.min_raw_cigar_insert_len_bp =
        config_.min_final_raw_cigar_insert_len_bp;
    final_filter_config.report_mode = config_.final_report_mode;
    const auto finalization_started = std::chrono::steady_clock::now();
    finalize_final_calls(
        result,
        tsd_detector_,
        config_.final_fdr_q,
        final_filter_config);
    result.performance.finalization_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - finalization_started).count();
    result.performance.gate1_preliminary_seconds += gate1_preliminary_seconds;
    result.performance.exact_bin_scan_seconds += exact_bin_scan_seconds;
    result.performance.pipeline_wall_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - pipeline_started).count();
    return result;
}
