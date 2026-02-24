add_executable(test_decision_policy
    test_decision_policy.cpp
)
target_include_directories(test_decision_policy PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_decision_policy PRIVATE placer_pipeline)

add_test(NAME test_decision_policy COMMAND test_decision_policy)
