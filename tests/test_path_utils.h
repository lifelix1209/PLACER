#ifndef PLACER_TEST_PATH_UTILS_H
#define PLACER_TEST_PATH_UTILS_H

#include <cstdlib>
#include <filesystem>
#include <initializer_list>
#include <iostream>
#include <string>

namespace placer_test {

inline std::string absolute_if_exists(const std::filesystem::path& path) {
    if (path.empty() || !std::filesystem::exists(path)) {
        return {};
    }
    return std::filesystem::absolute(path).string();
}

inline std::string find_from_cwd_upwards(const std::filesystem::path& relative_path,
                                         int max_levels = 8) {
    if (relative_path.empty()) return {};

    if (relative_path.is_absolute()) {
        return absolute_if_exists(relative_path);
    }

    std::filesystem::path current = std::filesystem::current_path();
    for (int level = 0; level <= max_levels; ++level) {
        const auto candidate = current / relative_path;
        auto found = absolute_if_exists(candidate);
        if (!found.empty()) {
            return found;
        }
        if (!current.has_parent_path()) break;
        current = current.parent_path();
    }
    return {};
}

inline std::string resolve_test_file(const char* env_name,
                                     std::initializer_list<const char*> repo_candidates) {
    if (env_name != nullptr) {
        const char* env_value = std::getenv(env_name);
        if (env_value != nullptr && env_value[0] != '\0') {
            auto env_path = absolute_if_exists(std::filesystem::path(env_value));
            if (!env_path.empty()) return env_path;
            std::cout << "  [WARN] " << env_name << " is set but file not found: "
                      << env_value << std::endl;
        }
    }

    for (const auto* candidate : repo_candidates) {
        if (candidate == nullptr || candidate[0] == '\0') continue;
        auto resolved = find_from_cwd_upwards(candidate);
        if (!resolved.empty()) {
            return resolved;
        }
    }

    return {};
}

inline bool require_path_or_skip(const std::string& path,
                                 const char* label,
                                 const char* env_name) {
    if (!path.empty()) return true;
    std::cout << "  [SKIP] Missing " << label;
    if (env_name != nullptr) {
        std::cout << " (set " << env_name << " to enable this test)";
    }
    std::cout << std::endl;
    return false;
}

inline std::string make_temp_dir(const std::string& dirname) {
    auto path = std::filesystem::temp_directory_path() / dirname;
    std::filesystem::create_directories(path);
    return path.string();
}

}  // namespace placer_test

#endif  // PLACER_TEST_PATH_UTILS_H
