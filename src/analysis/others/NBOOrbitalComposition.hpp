#ifndef TINKER_NBOORBITALCOMPOSITION_HPP
#define TINKER_NBOORBITALCOMPOSITION_HPP

#include <boost/fusion/sequence.hpp>
#include <boost/optional.hpp>
#include <boost/variant.hpp>

#include "utils/std.hpp"

class NBOOrbitalComposition {
public:
    [[nodiscard]] static std::string_view title() { return "NBO Orbital Composition Analysis"; }

    static void process();

    [[nodiscard]] static int findOccupancy(std::istream &is);

    [[nodiscard]] static bool match(const std::string &line);

    struct AtomComposition {
        int atom_no;
        std::string symbol;
        std::string type_name;
        std::string orbial_name;
        double contribution;
    };

    static void driveMultiwfn(const std::string &file, int alpha_orbitals, int beta_orbitals);

    [[nodiscard]] static std::optional<boost::fusion::vector<int, std::string, std::string, std::string, double>>
    parseLine(const std::string &line);

    [[nodiscard]] static std::map<int, std::vector<AtomComposition>, std::greater<>> read_contributions(
        std::istream &is, int orbital_number);

    static void print_contributions(std::string_view descriptions,
                                    std::map<int, std::vector<double>, std::greater<>> &contributions,
                                    const std::vector<std::string> &column_names);

    [[nodiscard]] static std::pair<std::map<int, std::vector<double>, std::greater<>>, std::vector<std::string>> filter(
        const std::vector<boost::fusion::vector<std::vector<boost::fusion::vector<int, boost::optional<int>>>,
                                                boost::optional<std::string>,
                                                boost::variant<std::vector<std::string>, std::string>>> &attrs,
        const std::map<int, std::vector<AtomComposition>, std::greater<>> &contributions);
};

#endif  // TINKER_NBOORBITALCOMPOSITION_HPP
