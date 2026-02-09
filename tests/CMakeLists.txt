# Test configuration for Phase 1, 2 & 3

add_executable(test_stream
    test_stream.cpp
)

target_link_libraries(test_stream PRIVATE
    stream
    gate1
    hts::htslib
)

target_include_directories(test_stream PRIVATE
    ${HTSLIB_INCLUDE_DIR}
)

add_test(NAME StreamTest COMMAND test_stream)

# Phase 3 Component Builder Tests
add_executable(test_component
    test_component.cpp
)

target_link_libraries(test_component PRIVATE
    stream
    gate1
    component
    hts::htslib
)

target_include_directories(test_component PRIVATE
    ${HTSLIB_INCLUDE_DIR}
)

add_test(NAME ComponentTest COMMAND test_component)

# Phase 4 Local Realigner Tests
add_executable(test_local_realign
    test_local_realign.cpp
)

target_link_libraries(test_local_realign PRIVATE
    stream
    gate1
    component
    local_realign
    hts::htslib
)

target_include_directories(test_local_realign PRIVATE
    ${HTSLIB_INCLUDE_DIR}
)

add_test(NAME LocalRealignTest COMMAND test_local_realign)

# Phase 5 Assembly + Collapsing Tests
add_executable(test_assembly
    test_assembly.cpp
)

target_link_libraries(test_assembly PRIVATE
    stream
    gate1
    component
    local_realign
    assembly
    hts::htslib
)

target_include_directories(test_assembly PRIVATE
    ${HTSLIB_INCLUDE_DIR}
)

add_test(NAME AssemblyTest COMMAND test_assembly)
