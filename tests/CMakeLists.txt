add_executable(test_decision_policy
    test_decision_policy.cpp
)
target_include_directories(test_decision_policy PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_decision_policy PRIVATE placer_pipeline)

add_test(NAME test_decision_policy COMMAND test_decision_policy)

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

add_executable(test_no_softclip_rescue
    test_no_softclip_rescue.cpp
)
target_include_directories(test_no_softclip_rescue PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_no_softclip_rescue PRIVATE placer_pipeline)

add_test(NAME test_no_softclip_rescue COMMAND test_no_softclip_rescue)

add_executable(test_short_insertion_rescue
    test_short_insertion_rescue.cpp
)
target_include_directories(test_short_insertion_rescue PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_short_insertion_rescue PRIVATE placer_pipeline)

add_test(NAME test_short_insertion_rescue COMMAND test_short_insertion_rescue)
