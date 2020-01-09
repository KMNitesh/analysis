#ifndef TINKER_TRAJCONVERTER_HPP
#define TINKER_TRAJCONVERTER_HPP

#include <string_view>

class TrajConverter {
public:
    [[nodiscard]] static std::string_view title() { return "Traj Converter (trj -> pdb)"; }

    static void process();
};


#endif //TINKER_TRAJCONVERTER_HPP
