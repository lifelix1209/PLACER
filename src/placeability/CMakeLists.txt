# Placeability Module CMakeLists.txt

add_library(placeability
    placeability.cpp
)

target_link_libraries(placeability PUBLIC
    stream
    gate1
    component
    local_realign
    assembly
    hts::htslib
)

target_include_directories(placeability PUBLIC
    ${PROJECT_SOURCE_DIR}/include
    ${HTSLIB_INCLUDE_DIR}
)

# Precompiled header (optional)
if(USE_PRECOMPILED_HEADER)
    target_precompile_headers(placeability
        <vector>
        <string>
        <cmath>
        <algorithm>
        <numeric>
        <sstream>
        <iomanip>
        <unordered_map>
        <memory>
    )
endif()
