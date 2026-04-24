find_package(Python3 REQUIRED COMPONENTS Interpreter)

add_executable(test_decision_policy
    test_decision_policy.cpp
)
target_include_directories(test_decision_policy PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_decision_policy PRIVATE placer_pipeline)

add_test(NAME test_decision_policy COMMAND test_decision_policy)

add_executable(test_te_alignment_decision_policy
    test_te_alignment_decision_policy.cpp
)
target_include_directories(test_te_alignment_decision_policy PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_te_alignment_decision_policy PRIVATE placer_pipeline)

add_test(NAME test_te_alignment_decision_policy COMMAND test_te_alignment_decision_policy)

add_executable(test_joint_te_emission_alt_only
    test_joint_te_emission_alt_only.cpp
)
target_include_directories(test_joint_te_emission_alt_only PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_joint_te_emission_alt_only PRIVATE placer_pipeline)

add_test(NAME test_joint_te_emission_alt_only COMMAND test_joint_te_emission_alt_only)

add_executable(test_bayesian_sv_genotyping
    test_bayesian_sv_genotyping.cpp
)
target_include_directories(test_bayesian_sv_genotyping PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_bayesian_sv_genotyping PRIVATE placer_pipeline)

add_test(NAME test_bayesian_sv_genotyping COMMAND test_bayesian_sv_genotyping)

add_executable(test_gate1
    test_gate1.cpp
)
target_include_directories(test_gate1 PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_gate1 PRIVATE placer_pipeline)

add_test(NAME test_gate1 COMMAND test_gate1)

add_executable(test_te_quick_classifier_multik
    test_te_quick_classifier_multik.cpp
)
target_include_directories(test_te_quick_classifier_multik PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_te_quick_classifier_multik PRIVATE placer_pipeline)

add_test(NAME test_te_quick_classifier_multik COMMAND test_te_quick_classifier_multik)

add_executable(test_tsd_detector
    test_tsd_detector.cpp
)
target_include_directories(test_tsd_detector PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_tsd_detector PRIVATE placer_pipeline)

add_test(NAME test_tsd_detector COMMAND test_tsd_detector)

add_executable(test_event_segmentation
    test_event_segmentation.cpp
)
target_include_directories(test_event_segmentation PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_event_segmentation PRIVATE placer_pipeline)

add_test(NAME test_event_segmentation COMMAND test_event_segmentation)

add_executable(test_breakpoint_hypothesis_enumeration
    test_breakpoint_hypothesis_enumeration.cpp
)
target_include_directories(test_breakpoint_hypothesis_enumeration PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_breakpoint_hypothesis_enumeration PRIVATE placer_pipeline)

add_test(NAME test_breakpoint_hypothesis_enumeration COMMAND test_breakpoint_hypothesis_enumeration)

add_executable(test_anchor_seed_bins
    test_anchor_seed_bins.cpp
)
target_include_directories(test_anchor_seed_bins PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_anchor_seed_bins PRIVATE placer_pipeline)

add_test(NAME test_anchor_seed_bins COMMAND test_anchor_seed_bins)

add_executable(test_event_read_evidence
    test_event_read_evidence.cpp
)
target_include_directories(test_event_read_evidence PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_event_read_evidence PRIVATE placer_pipeline)

add_test(NAME test_event_read_evidence COMMAND test_event_read_evidence)

add_executable(test_pre_segmentation_gate
    test_pre_segmentation_gate.cpp
)
target_include_directories(test_pre_segmentation_gate PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_pre_segmentation_gate PRIVATE placer_pipeline)

add_test(NAME test_pre_segmentation_gate COMMAND test_pre_segmentation_gate)

add_executable(test_insert_te_alignment
    test_insert_te_alignment.cpp
)
target_include_directories(test_insert_te_alignment PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_insert_te_alignment PRIVATE placer_pipeline)

add_test(NAME test_insert_te_alignment COMMAND test_insert_te_alignment)

add_executable(test_te_classifier_shared_cache
    test_te_classifier_shared_cache.cpp
)
target_include_directories(test_te_classifier_shared_cache PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_te_classifier_shared_cache PRIVATE placer_pipeline)

add_test(NAME test_te_classifier_shared_cache COMMAND test_te_classifier_shared_cache)

add_executable(test_te_family_prefilter
    test_te_family_prefilter.cpp
)
target_include_directories(test_te_family_prefilter PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src/component
)
target_link_libraries(test_te_family_prefilter PRIVATE placer_pipeline)

add_test(NAME test_te_family_prefilter COMMAND test_te_family_prefilter)

add_executable(test_te_repetitive_seed_downweight
    test_te_repetitive_seed_downweight.cpp
)
target_include_directories(test_te_repetitive_seed_downweight PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src/component
)
target_link_libraries(test_te_repetitive_seed_downweight PRIVATE placer_pipeline)

add_test(NAME test_te_repetitive_seed_downweight COMMAND test_te_repetitive_seed_downweight)

add_executable(test_te_exact_alignment_budget
    test_te_exact_alignment_budget.cpp
)
target_include_directories(test_te_exact_alignment_budget PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_te_exact_alignment_budget PRIVATE placer_pipeline)

add_test(NAME test_te_exact_alignment_budget COMMAND test_te_exact_alignment_budget)

add_executable(bench_te_library_alignment
    bench_te_library_alignment.cpp
)
target_include_directories(bench_te_library_alignment PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(bench_te_library_alignment PRIVATE placer_pipeline)

add_executable(test_parallel_exact_bin_task_equivalence
    test_parallel_exact_bin_task_equivalence.cpp
)
target_include_directories(test_parallel_exact_bin_task_equivalence PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src/pipeline
)
target_link_libraries(test_parallel_exact_bin_task_equivalence PRIVATE placer_pipeline)

add_test(
    NAME test_parallel_exact_bin_task_equivalence
    COMMAND test_parallel_exact_bin_task_equivalence
)

add_executable(test_parallel_bounded_scheduler
    test_parallel_bounded_scheduler.cpp
)
target_include_directories(test_parallel_bounded_scheduler PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src/pipeline
)
target_link_libraries(test_parallel_bounded_scheduler PRIVATE placer_pipeline)

add_test(
    NAME test_parallel_bounded_scheduler
    COMMAND test_parallel_bounded_scheduler
)

add_executable(test_parallel_local_interval_cache
    test_parallel_local_interval_cache.cpp
)
target_include_directories(test_parallel_local_interval_cache PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src/pipeline
)
target_link_libraries(test_parallel_local_interval_cache PRIVATE placer_pipeline)

add_test(
    NAME test_parallel_local_interval_cache
    COMMAND test_parallel_local_interval_cache
)

add_executable(test_validator_input_counts
    test_validator_input_counts.cpp
)
target_include_directories(test_validator_input_counts PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_validator_input_counts PRIVATE placer_pipeline)

add_test(
    NAME test_validator_input_counts
    COMMAND test_validator_input_counts
)

add_executable(test_component_signal_cache
    test_component_signal_cache.cpp
)
target_include_directories(test_component_signal_cache PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_component_signal_cache PRIVATE placer_pipeline)

add_test(
    NAME test_component_signal_cache
    COMMAND test_component_signal_cache
)

add_executable(test_component_locality
    test_component_locality.cpp
)
target_include_directories(test_component_locality PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_component_locality PRIVATE placer_pipeline)

add_test(NAME test_component_locality COMMAND test_component_locality)

add_executable(test_component_multi_emit_selection
    test_component_multi_emit_selection.cpp
)
target_include_directories(test_component_multi_emit_selection PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_component_multi_emit_selection PRIVATE placer_pipeline)

add_test(NAME test_component_multi_emit_selection COMMAND test_component_multi_emit_selection)

add_executable(test_hypothesis_pruning
    test_hypothesis_pruning.cpp
)
target_include_directories(test_hypothesis_pruning PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_hypothesis_pruning PRIVATE placer_pipeline)

add_test(NAME test_hypothesis_pruning COMMAND test_hypothesis_pruning)

add_executable(test_dbscan_component_module
    test_dbscan_component_module.cpp
)
target_include_directories(test_dbscan_component_module PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_dbscan_component_module PRIVATE placer_pipeline)

add_test(NAME test_dbscan_component_module COMMAND test_dbscan_component_module)

add_executable(test_denovo_cli
    test_denovo_cli.cpp
)
target_include_directories(test_denovo_cli PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_denovo_cli PRIVATE placer_denovo)

add_test(NAME test_denovo_cli COMMAND test_denovo_cli)

add_test(
    NAME test_run_sharded_placer
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_run_sharded_placer.py
)

add_test(
    NAME test_release_readiness_gate
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_release_readiness_gate.py
)
