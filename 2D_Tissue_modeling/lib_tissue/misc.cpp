#include  "util.h"

// template <typename T>
std::string get_str_from_float(double data) {

    std::ostringstream strs;
    strs << data;
    std::string str = strs.str();
    return str;
}


