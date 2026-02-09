# Component Builder Module

add_library(component STATIC
    component_builder.cpp
)

target_include_directories(component PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${HTSLIB_INCLUDE_DIR}
)

target_link_libraries(component PRIVATE
    stream
    hts::htslib
)
