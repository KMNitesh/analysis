#ifndef TINKER_TRAJCONVERTER_HPP
#define TINKER_TRAJCONVERTER_HPP

#include <string_view>

class TrajConverter {
public:
    [[nodiscard]] static std::string_view title() { return "Traj Converter (trj -> pdb)"; }

    static void process();

    static int parse_header(const std::string &line);

    static boost::fusion::vector<std::string, double, double, double> parse_atom(const std::string &line);
};


#endif //TINKER_TRAJCONVERTER_HPP
