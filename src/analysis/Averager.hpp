

#ifndef TINKER_AVERAGER_HPP
#define TINKER_AVERAGER_HPP


#include <string_view>

class Averager {
public:

    static void process();

    [[nodiscard]] static std::string_view title() { return "Averager"; }
};


#endif //TINKER_AVERAGER_HPP
