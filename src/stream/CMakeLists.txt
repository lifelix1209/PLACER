# Stream Layer Module

add_library(stream STATIC
    bam_reader.cpp
    window_buffer.cpp
    trigger.cpp
    task_queue.cpp
)

target_include_directories(stream PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)

# Link with system htslib
target_include_directories(stream PRIVATE ${HTSLIB_INCLUDE_DIR})
target_link_libraries(stream PRIVATE hts::htslib)
