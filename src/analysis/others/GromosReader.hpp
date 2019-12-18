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
    };
};


#endif //TINKER_GROMOSREADER_HPP
