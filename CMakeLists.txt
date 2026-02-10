cmake_minimum_required(VERSION 3.16)
project(PLACER VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Use pysam's bundled htslib headers
execute_process(
    COMMAND python3 -c "import pysam; inc = pysam.get_include(); print(inc[1] if len(inc) > 1 else inc[0])"
    OUTPUT_VARIABLE HTSLIB_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

message(STATUS "Found htslib headers: ${HTSLIB_INCLUDE_DIR}")

# System htslib library path
set(HTSLIB_LIBRARY "/usr/lib/x86_64-linux-gnu/libhts.so.3")

if(NOT EXISTS ${HTSLIB_LIBRARY})
    message(FATAL_ERROR "Cannot find htslib library at ${HTSLIB_LIBRARY}")
endif()

message(STATUS "Found htslib library: ${HTSLIB_LIBRARY}")

# Create imported target for htslib
add_library(hts::htslib STATIC IMPORTED)
set_target_properties(hts::htslib PROPERTIES
    IMPORTED_LOCATION "${HTSLIB_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORY "${HTSLIB_INCLUDE_DIR}"
)

# Main executable
add_executable(placer src/main.cpp)

target_include_directories(placer PRIVATE ${HTSLIB_INCLUDE_DIR})
target_link_libraries(placer PRIVATE hts::htslib stream gate1 component local_realign assembly placeability genotyping te_reverse_index)

# Add subdirectories
add_subdirectory(src/stream)
add_subdirectory(src/gate1)
add_subdirectory(src/component)
add_subdirectory(src/local_realign)
add_subdirectory(src/assembly)
add_subdirectory(src/placeability)
add_subdirectory(src/genotyping)

# Test configuration
enable_testing()
add_subdirectory(tests)

# Add TE Reverse Index module (Phase 8)
add_subdirectory(src/te_reverse_index)
