#ifndef TINKER_GROMOSREADER_HPP
#define TINKER_GROMOSREADER_HPP

#include <string_view>

class GromosReader {
public:
    [[nodiscard]] static std::string_view title() { return "GROMOS omd reader"; }

    static void process();

    static void printMenu(std::vector<std::string> &menuStrings, std::size_t width);

    struct Energy {
        uint step;
        double time;
        std::vector<double> energies;

        std::array<std::vector<double>, 3> intergroups;
    };

    static void save2Xml(const std::vector<Energy> &energies,
                         const std::vector<std::string> &menuStrings,
                         const std::array<std::string, 3> &energy_names,
                         const std::vector<std::string> &group_names,
                         std::ostream &os = std::cout);
};


#endif //TINKER_GROMOSREADER_HPP
