# TE Reverse Index library
add_library(te_reverse_index STATIC te_reverse_index.cpp)
target_include_directories(te_reverse_index PRIVATE ${HTSLIB_INCLUDE_DIR})
target_link_libraries(te_reverse_index PRIVATE hts::htslib gate1 component)
