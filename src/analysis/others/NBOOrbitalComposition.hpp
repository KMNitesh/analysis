#ifndef TINKER_NBOORBITALCOMPOSITION_HPP
#define TINKER_NBOORBITALCOMPOSITION_HPP

#include <string_view>
#include <optional>
#include <boost/fusion/sequence.hpp>

class NBOOrbitalComposition {
public:
    [[nodiscard]] static std::string_view title() { return "NBO Orbital Composition Analysis"; }

    static void process();

    [[nodiscard]] static int findOccupancy(std::istream &is);

    [[nodiscard]] static bool match(const std::string &line);

    static void driveMultiwfn(const std::string &file, int alpha_orbitals, int beta_orbitals);

    [[nodiscard]] static std::optional<boost::fusion::vector<int, std::string, std::string, double>>
    parseLine(const std::string &line);

    [[nodiscard]] static std::map<int, std::map<std::pair<int, std::string>, double>, std::greater<>>
    read_contributions(const std::vector<boost::fusion::vector<int, std::string>> &attrs, std::istream &is,
                       int orbital_number);

    static void print_contributions(
            std::string_view descriptions,
            std::map<int, std::map<std::pair<int, std::string>, double>, std::greater<>> &contributions,
            const std::vector<boost::fusion::vector<int, std::string>> &attrs);

};


#endif //TINKER_NBOORBITALCOMPOSITION_HPP
