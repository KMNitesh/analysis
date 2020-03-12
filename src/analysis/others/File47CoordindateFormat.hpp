#ifndef TINKER_FILE47COORDINDATEFORMAT_HPP
#define TINKER_FILE47COORDINDATEFORMAT_HPP

#include <string_view>

class File47CoordindateFormat {
public:
    [[nodiscard]] static std::string_view title() { return "File .47 Coordinate Format"; }

    static void process();
};

#endif  // TINKER_FILE47COORDINDATEFORMAT_HPP
