cmake_minimum_required(VERSION 3.16)
project(PLACER VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# ============================================================================
# htslib discovery
# ============================================================================
set(HTSLIB_ROOT "" CACHE PATH "Root directory of htslib installation")
set(HTSLIB_INCLUDE_DIR "" CACHE PATH "Directory containing htslib/sam.h")
set(HTSLIB_LIBRARY "" CACHE FILEPATH "Path to libhts library file")

set(_HTSLIB_INCLUDE_HINTS)
if(NOT HTSLIB_INCLUDE_DIR)
    execute_process(
        COMMAND python3 -c "import pysam; print(';'.join(pysam.get_include()))"
        OUTPUT_VARIABLE _PYSAM_INCLUDE_HINTS
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE _PYSAM_INCLUDE_STATUS
        ERROR_QUIET
    )
    if(_PYSAM_INCLUDE_STATUS EQUAL 0 AND _PYSAM_INCLUDE_HINTS)
        list(APPEND _HTSLIB_INCLUDE_HINTS ${_PYSAM_INCLUDE_HINTS})
    endif()
endif()

set(_HTSLIB_HINT_ROOTS
    ${HTSLIB_ROOT}
    $ENV{HTSLIB_ROOT}
    $ENV{CONDA_PREFIX}
    /opt/homebrew
    /usr/local
    /usr
    /opt/anaconda3
)

if(NOT HTSLIB_INCLUDE_DIR)
    find_path(_HTSLIB_INCLUDE_DIR_FOUND
        NAMES htslib/sam.h
        HINTS ${_HTSLIB_INCLUDE_HINTS}
        PATHS ${_HTSLIB_HINT_ROOTS}
        PATH_SUFFIXES include include/htslib htslib
    )
    if(_HTSLIB_INCLUDE_DIR_FOUND)
        set(HTSLIB_INCLUDE_DIR "${_HTSLIB_INCLUDE_DIR_FOUND}" CACHE PATH
            "Directory containing htslib/sam.h" FORCE)
    endif()
endif()

if(NOT HTSLIB_LIBRARY)
    find_library(_HTSLIB_LIBRARY_FOUND
        NAMES hts libhts
        PATHS ${_HTSLIB_HINT_ROOTS}
        PATH_SUFFIXES lib lib64
    )
    if(_HTSLIB_LIBRARY_FOUND)
        set(HTSLIB_LIBRARY "${_HTSLIB_LIBRARY_FOUND}" CACHE FILEPATH
            "Path to libhts library file" FORCE)
    endif()
endif()

if(NOT HTSLIB_INCLUDE_DIR OR NOT HTSLIB_LIBRARY)
    message(FATAL_ERROR
        "Cannot find htslib.\n"
        "Detected HTSLIB_INCLUDE_DIR='${HTSLIB_INCLUDE_DIR}'\n"
        "Detected HTSLIB_LIBRARY='${HTSLIB_LIBRARY}'\n"
        "Install htslib or pass explicit paths:\n"
        "  -DHTSLIB_ROOT=<prefix>\n"
        "  -DHTSLIB_INCLUDE_DIR=<path containing htslib/sam.h>\n"
        "  -DHTSLIB_LIBRARY=<full path to libhts>"
    )
endif()

message(STATUS "Found htslib headers: ${HTSLIB_INCLUDE_DIR}")
message(STATUS "Found htslib library: ${HTSLIB_LIBRARY}")

add_library(hts::htslib UNKNOWN IMPORTED)
set_target_properties(hts::htslib PROPERTIES
    IMPORTED_LOCATION "${HTSLIB_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${HTSLIB_INCLUDE_DIR}"
)

# ============================================================================
# PLACER modules
# ============================================================================
add_library(placer_bam_io
    src/stream/bam_io.cpp
)
target_include_directories(placer_bam_io
    PUBLIC
        include
        ${HTSLIB_INCLUDE_DIR}
)
target_link_libraries(placer_bam_io PUBLIC hts::htslib)

add_library(placer_pipeline
    src/pipeline/pipeline.cpp
    src/gate1/gate1_module.cpp
    src/component/insert_fragment_module.cpp
    src/component/te_quick_classifier.cpp
)
target_include_directories(placer_pipeline PUBLIC include)
target_link_libraries(placer_pipeline PUBLIC placer_bam_io)

add_executable(placer
    src/main.cpp
)
target_include_directories(placer PRIVATE include)
target_link_libraries(placer PRIVATE placer_pipeline)
