# Gate 1 Module

add_library(gate1 STATIC
    gate1.cpp
    gate1_filter.cpp
)

target_include_directories(gate1 PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${HTSLIB_INCLUDE_DIR}
)

target_link_libraries(gate1 PRIVATE
    hts::htslib
    stream
)
