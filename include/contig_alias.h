#ifndef PLACER_CONTIG_ALIAS_H
#define PLACER_CONTIG_ALIAS_H

#include <algorithm>
#include <cctype>
#include <string>
#include <vector>

namespace placer {

namespace contig_alias_detail {

inline bool starts_with_chr_prefix(const std::string& chrom) {
    return chrom.size() > 3 &&
           std::tolower(static_cast<unsigned char>(chrom[0])) == 'c' &&
           std::tolower(static_cast<unsigned char>(chrom[1])) == 'h' &&
           std::tolower(static_cast<unsigned char>(chrom[2])) == 'r';
}

inline std::string upper_ascii(std::string s) {
    for (char& c : s) {
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    }
    return s;
}

inline void append_unique(std::vector<std::string>& names, const std::string& name) {
    if (name.empty()) {
        return;
    }
    if (std::find(names.begin(), names.end(), name) == names.end()) {
        names.push_back(name);
    }
}

}  // namespace contig_alias_detail

inline std::vector<std::string> contig_name_aliases(const std::string& chrom) {
    std::vector<std::string> aliases;
    contig_alias_detail::append_unique(aliases, chrom);
    if (chrom.empty()) {
        return aliases;
    }

    const bool has_chr = contig_alias_detail::starts_with_chr_prefix(chrom);
    const std::string stripped = has_chr ? chrom.substr(3) : chrom;
    const std::string stripped_upper = contig_alias_detail::upper_ascii(stripped);

    if (has_chr) {
        contig_alias_detail::append_unique(aliases, stripped);
    } else {
        contig_alias_detail::append_unique(aliases, "chr" + chrom);
    }

    if (stripped_upper == "M") {
        contig_alias_detail::append_unique(aliases, "MT");
        contig_alias_detail::append_unique(aliases, "chrM");
        contig_alias_detail::append_unique(aliases, "chrMT");
    } else if (stripped_upper == "MT") {
        contig_alias_detail::append_unique(aliases, "M");
        contig_alias_detail::append_unique(aliases, "chrMT");
        contig_alias_detail::append_unique(aliases, "chrM");
    }

    return aliases;
}

}  // namespace placer

#endif  // PLACER_CONTIG_ALIAS_H
