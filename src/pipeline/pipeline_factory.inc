std::unique_ptr<Pipeline> build_default_pipeline(const PipelineConfig& config) {
    auto reader = make_bam_reader(config.bam_path, config.bam_threads, config.bam_region_scope);
    return std::make_unique<Pipeline>(config, std::move(reader));
}
