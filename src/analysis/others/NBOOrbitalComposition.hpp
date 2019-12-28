#ifndef TINKER_NBOORBITALCOMPOSITION_HPP
#define TINKER_NBOORBITALCOMPOSITION_HPP

#include <string_view>
#include <optional>
#include <boost/fusion/sequence.hpp>

class NBOOrbitalComposition {
public:
    [[nodiscard]] static std::string_view title() { return "NBO Oribital Composition Analysis"; }

    static void process();

    static int findOccupancy(std::istream &is);

    static bool match(const std::string &line);

    static void driveMultiwfn(const std::string &file, int alpha_orbitals);

    static std::optional<boost::fusion::vector<int, std::string, std::string, double>>
    parseLine(const std::string &line);

};


#endif //TINKER_NBOORBITALCOMPOSITION_HPP
