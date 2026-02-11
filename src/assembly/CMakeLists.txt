# Assembly module (Phase 5: Assembly + Collapsing)

# Headers
set(ASSEMBLY_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/../../include/assembly.h
)

# Sources
set(ASSEMBLY_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/assembly.cpp
)

# Create static library
add_library(assembly STATIC
    ${ASSEMBLY_SOURCES}
    ${ASSEMBLY_HEADERS}
)

# Include directories
target_include_directories(assembly PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/../../include
    ${CMAKE_SOURCE_DIR}/third_party/abPOA/src
)

# Link dependencies
target_link_libraries(assembly
    PUBLIC
        component
        local_realign
        abpoa_wrapper
)

# Compile options
target_compile_options(assembly PRIVATE
    ${PLACER_CXX_FLAGS}
)

# Installation
install(TARGETS assembly
    ARCHIVE DESTINATION lib
)
