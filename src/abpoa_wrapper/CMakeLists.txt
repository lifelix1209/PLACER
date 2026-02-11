# abPOA Wrapper Module

add_library(abpoa_wrapper STATIC
    abpoa_wrapper.cpp
)

target_include_directories(abpoa_wrapper PRIVATE
    ${CMAKE_SOURCE_DIR}/include
)

if(TARGET abpoa)
    target_include_directories(abpoa_wrapper PRIVATE
        ${CMAKE_SOURCE_DIR}/third_party/abPOA/src
        ${CMAKE_SOURCE_DIR}/third_party/abPOA/include
    )
    target_compile_definitions(abpoa_wrapper PUBLIC PLACER_HAS_ABPOA=1)
    target_link_libraries(abpoa_wrapper PUBLIC abpoa)
else()
    target_compile_definitions(abpoa_wrapper PUBLIC PLACER_HAS_ABPOA=0)
endif()
