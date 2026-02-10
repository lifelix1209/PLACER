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

# Phase 6 Placeability + Tier Tests
add_executable(test_placeability
    test_placeability.cpp
)

target_link_libraries(test_placeability PRIVATE
    stream
    gate1
    component
    local_realign
    assembly
    placeability
    hts::htslib
)

target_include_directories(test_placeability PRIVATE
    ${HTSLIB_INCLUDE_DIR}
)

add_test(NAME PlaceabilityTest COMMAND test_placeability)

# Placeability Integration Tests (using tldr_optimized/test data)
add_executable(test_placeability_integration
    test_placeability_integration.cpp
)

target_link_libraries(test_placeability_integration PRIVATE
    stream
    gate1
    component
    local_realign
    assembly
    placeability
    hts::htslib
)

target_include_directories(test_placeability_integration PRIVATE
    ${HTSLIB_INCLUDE_DIR}
)

add_test(NAME PlaceabilityIntegrationTest COMMAND test_placeability_integration)

# Phase 7 Genotyping Tests
add_executable(test_genotyping
    test_genotyping.cpp
)

target_link_libraries(test_genotyping PRIVATE
    stream
    gate1
    component
    local_realign
    assembly
    placeability
    genotyping
    hts::htslib
)

target_include_directories(test_genotyping PRIVATE
    ${HTSLIB_INCLUDE_DIR}
)

add_test(NAME GenotypingTest COMMAND test_genotyping)

# Phase 8 TE Reverse Index Tests
add_executable(test_te_reverse_index
    test_te_reverse_index.cpp
)

target_link_libraries(test_te_reverse_index PRIVATE
    stream
    gate1
    component
    local_realign
    assembly
    placeability
    genotyping
    te_reverse_index
    hts::htslib
)

target_include_directories(test_te_reverse_index PRIVATE
    ${HTSLIB_INCLUDE_DIR}
)

add_test(NAME TEReverseIndexTest COMMAND test_te_reverse_index)
