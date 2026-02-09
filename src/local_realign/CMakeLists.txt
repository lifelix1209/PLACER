# Local Realigner Module

add_library(local_realign STATIC
    local_realign.cpp
)

target_include_directories(local_realign PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${HTSLIB_INCLUDE_DIR}
)

target_link_libraries(local_realign PRIVATE
    stream
    hts::htslib
)
