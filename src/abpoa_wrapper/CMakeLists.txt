# abPOA Wrapper Module

add_library(abpoa_wrapper STATIC
    abpoa_wrapper.cpp
)

target_include_directories(abpoa_wrapper PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/third_party/abPOA/src
)

target_link_libraries(abpoa_wrapper PUBLIC abpoa)
