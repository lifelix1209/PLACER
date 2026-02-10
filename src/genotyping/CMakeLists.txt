# Genotyping module CMakeLists.txt

add_library(genotyping STATIC
    genotyping.cpp
)

target_include_directories(genotyping PRIVATE
    ${HTSLIB_INCLUDE_DIR}
)

target_link_libraries(genotyping PRIVATE
    placeability
    assembly
    local_realign
    hts::htslib
)
