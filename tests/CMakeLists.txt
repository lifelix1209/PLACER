# Test configuration for Phase 1

add_executable(test_stream
    test_stream.cpp
)

target_link_libraries(test_stream PRIVATE
    stream
    hts::htslib
)

target_include_directories(test_stream PRIVATE
    ${HTSLIB_INCLUDE_DIR}
)

add_test(NAME StreamTest COMMAND test_stream)
