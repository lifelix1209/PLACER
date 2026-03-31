find_package(Python3 REQUIRED COMPONENTS Interpreter)

add_executable(test_decision_policy
    test_decision_policy.cpp
)
target_include_directories(test_decision_policy PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_decision_policy PRIVATE placer_pipeline)

add_test(NAME test_decision_policy COMMAND test_decision_policy)

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

add_executable(test_event_read_evidence
    test_event_read_evidence.cpp
)
target_include_directories(test_event_read_evidence PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_event_read_evidence PRIVATE placer_pipeline)

add_test(NAME test_event_read_evidence COMMAND test_event_read_evidence)

add_executable(test_insert_te_alignment
    test_insert_te_alignment.cpp
)
target_include_directories(test_insert_te_alignment PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_insert_te_alignment PRIVATE placer_pipeline)

add_test(NAME test_insert_te_alignment COMMAND test_insert_te_alignment)

add_executable(test_local_recollection_fragment_refresh
    test_local_recollection_fragment_refresh.cpp
)
target_include_directories(test_local_recollection_fragment_refresh PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_local_recollection_fragment_refresh PRIVATE placer_pipeline)

add_test(NAME test_local_recollection_fragment_refresh COMMAND test_local_recollection_fragment_refresh)

add_executable(test_ref_span_dominant_rejection
    test_ref_span_dominant_rejection.cpp
)
target_include_directories(test_ref_span_dominant_rejection PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_ref_span_dominant_rejection PRIVATE placer_pipeline)

add_test(NAME test_ref_span_dominant_rejection COMMAND test_ref_span_dominant_rejection)

add_executable(test_component_locality
    test_component_locality.cpp
)
target_include_directories(test_component_locality PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_component_locality PRIVATE placer_pipeline)

add_test(NAME test_component_locality COMMAND test_component_locality)

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
    NAME test_random_truth_interval_eval
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_random_truth_interval_eval.py
)

add_test(
    NAME test_submit_placer_slurm
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_submit_placer_slurm.py
)
