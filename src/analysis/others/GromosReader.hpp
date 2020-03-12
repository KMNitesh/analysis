#ifndef TINKER_GROMOSREADER_HPP
#define TINKER_GROMOSREADER_HPP

#include "utils/std.hpp"

class GromosReader {
public:
    [[nodiscard]] static std::string_view title() { return "GROMOS omd reader"; }

    static void process();

    static void printMenu(const std::vector<std::string> &menuStrings, std::size_t width);

    struct Energy {
        uint step;
        double time;
        std::vector<double> energies;

        std::array<std::vector<double>, 3> intergroups;
    };

    static void save2Xml(const std::vector<Energy> &energies, const std::vector<std::string> &menuStrings,
                         const std::array<std::string, 3> &energy_names, const std::vector<std::string> &group_names,
                         std::ostream &os = std::cout);

    struct Bundle {
        std::vector<Energy> energies;
        std::vector<std::string> menuStrings;
        std::array<std::string, 3> energy_names;
        std::vector<std::string> group_names;
    };

    static void printEnergies(const std::vector<Energy> &energies, const std::vector<std::string> &menuStrings,
                              const std::array<std::string, 3> &energy_names,
                              const std::vector<std::string> &group_names);

    [[nodiscard]] static Bundle readOmd(const std::string &filename);

    [[nodiscard]] static Bundle readXml(std::istream &is);

    [[nodiscard]] static Bundle readXml(const std::string &filename);
};

#endif  // TINKER_GROMOSREADER_HPP
