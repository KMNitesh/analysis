

#ifndef TINKER_AVERAGER_HPP
#define TINKER_AVERAGER_HPP

#include <iosfwd>
#include <string_view>
#include <vector>

class Averager {
public:
    static void process();

    [[nodiscard]] static std::string_view title() { return "Averager"; }

    static void process_file(std::ifstream &ifs, std::ofstream &ofs, int range);

    [[nodiscard]] static std::vector<std::vector<double>> readfile(std::ifstream &ifs);

    static void output_average(std::ostream &os, int range, const std::vector<std::vector<double>> &lines);
};

#endif  // TINKER_AVERAGER_HPP
