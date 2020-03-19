
#include "Averager.hpp"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/range/algorithm.hpp>

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

    int range = choose(2, "Enter range (line) > ");

    process_file(ifs, ofs, range);
}

void Averager::process_file(std::ifstream &ifs, std::ofstream &ofs, const int range) {
    std::vector<std::vector<double>> lines = readfile(ifs);

    output_average(ofs, range, lines);
}

void Averager::output_average(std::ostream &os, int range, const std::vector<std::vector<double>> &lines) {
    using namespace boost::accumulators;

    std::vector<accumulator_set<double, features<tag::mean, tag::variance>>> statics(lines.front().size());

    std::vector<double> avg(lines.front().size());

    for (std::size_t i = 0; i < lines.size(); ++i) {
        const auto &v = lines[i];
        for (std::size_t j = 0; j < v.size(); ++j) {
            avg[j] += v[j];
            statics[j](v[j]);
        }

        if ((i + 1) % range == 0) {
            for (auto value : avg) {
                os << (value / range) << "  ";
            }
            os << '\n';
            boost::fill(avg, 0.0);
        }
    }
    std::cout << "Mean : ";
    for (auto &s : statics)
        std::cout << std::setw(15) << mean(s) << "  ";
    std::cout << '\n';

    std::cout << "STD  : ";
    for (auto &s : statics)
        std::cout << std::setw(15) << std::sqrt(variance(s)) << "  ";
    std::cout << '\n';
}

std::vector<std::vector<double>> Averager::readfile(std::ifstream &ifs) {
    std::vector<std::vector<double>> lines;
    std::string line;

    while (!ifs.eof()) {
        std::getline(ifs, line);
        if (line.empty())
            continue;
        if (boost::starts_with(line, "#") or boost::starts_with(line, "@"))
            continue;

        auto field = split(line);

        std::vector<double> item;
        item.reserve(field.size());
        for (const auto &s : field) {
            item.emplace_back(std::stod(s));
        }
        lines.emplace_back(std::move(item));
    }
    return lines;
}
