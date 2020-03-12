#ifndef TINKER_GRORENUMBER_HPP
#define TINKER_GRORENUMBER_HPP

#include <string_view>

class GroRenumber {
public:
    [[nodiscard]] static std::string_view title() { return "Renumber GRO"; }

    static void process();
};

#endif // TINKER_GRORENUMBER_HPP
