#ifndef PLACER_BOOST_FORMAT_SHIM_HPP
#define PLACER_BOOST_FORMAT_SHIM_HPP

#include <string>

namespace boost {

class format {
public:
    explicit format(const char* fmt) : fmt_(fmt ? fmt : "") {}
    explicit format(const std::string& fmt) : fmt_(fmt) {}

    template <typename T>
    format& operator%(const T&) {
        return *this;
    }

    std::string str() const {
        return fmt_;
    }

private:
    std::string fmt_;
};

}  // namespace boost

#endif  // PLACER_BOOST_FORMAT_SHIM_HPP
