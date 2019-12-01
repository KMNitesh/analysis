
#include <boost/algorithm/string.hpp>
#include <boost/range/algorithm.hpp>
#include "Averager.hpp"
#include "utils/common.hpp"

void Averager::process() {

    std::string input_file = choose_file("Input xvg file > ").isExist(true);
    std::string output_file = choose_file("Output file > ").isExist(false);

    std::ifstream ifs(input_file);
    if (!ifs) {
        std::cerr << "ERROR !! Cannot open input file " << input_file << '\n';
        exit(EXIT_FAILURE);
    }

    std::ofstream ofs(output_file);
    if (!ofs) {
        std::cerr << "ERROR !! Cannot open output file " << output_file << '\n';
        exit(EXIT_FAILURE);
    }


    std::vector<std::vector<double>> lines;

    std::string line;

    while (!ifs.eof()) {
        std::getline(ifs, line);
        if (line.empty()) continue;
        auto field = split(line);

        std::vector<double> item;
        item.reserve(field.size());
        for (const auto &s : field) {
            item.emplace_back(std::stod(s));
        }
        lines.emplace_back(std::move(item));
    }

    int range = choose(2, "Enter range (line) > ");

    std::vector<double> avg;
    avg.resize(lines.front().size());
    for (std::size_t i = 0; i < lines.size(); ++i) {
        const auto &v = lines[i];
        for (std::size_t j = 0; j < v.size(); ++j) {
            avg[j] += v[j];
        }

        if ((i + 1) % range == 0) {
            for (auto value : avg) {
                ofs << (value / range) << "  ";
            }
            ofs << '\n';
            boost::fill(avg, 0.0);
        }
    }
}
