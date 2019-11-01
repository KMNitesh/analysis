//
// Created by xiamr on 11/1/19.
//

#include "NBOSpin.hpp"
#include "common.hpp"
#include <boost/xpressive/xpressive.hpp>
#include <boost/algorithm/string.hpp>

double NBOSpin::total_spin(std::string line) {
    using namespace boost::xpressive;

    static sregex rex = +~(set = '(') >> '(' >> _s >> (s1 = digit >> '.' >> digit >> digit) >> ')';
    smatch what;
    double total = 0.0;
    while (regex_search(line, what, rex)) {
        total += std::stod(what[1]);
        line = what.suffix().str();
    }
    return total;
}

void NBOSpin::do_process(const std::string &filename) {
    std::ifstream ifs(filename);
    if (!ifs) {
        std::cerr << "ERROR!! Cannot open file " << filename << '\n';
        exit(EXIT_FAILURE);
    }

    std::map<int, std::pair<std::string, std::array<double, 3>>> table = getElectronSpin(ifs);

    std::cout << boost::format("%6s %4s %6s %6s %6s %6s\n")
                 % "Atom" % "No" % "Total" % "Alpha" % "Beta" % "Alpha-Beta";
    for (auto &item : table) {
        auto diff = item.second.second[1] - item.second.second[2];
        std::cout << boost::format("%6s %4d %6.2f %6.2f %6.2f %6.2f  %s\n")
                     % item.second.first
                     % item.first
                     % item.second.second[0]
                     % item.second.second[1]
                     % item.second.second[2]
                     % diff
                     % (std::abs(diff) >= 0.1 ? '*' : ' ');
    }
}

std::map<int, std::pair<std::string, std::array<double, 3>>> NBOSpin::getElectronSpin(std::istream &ifs) {
    std::map<int, std::pair<std::string, std::array<double, 3>>> table;
    int current_pos = 0;
    std::string line;
    while (std::getline(ifs, line)) {
        if (boost::ends_with(line, "Atom No         Natural Electron Configuration")) {
            std::getline(ifs, line);
            for (;;) {
                std::getline(ifs, line);
                boost::trim(line);
                if (line.empty()) {
                    ++current_pos;
                    break;
                }
                auto field = split(line);
                auto &item = table[std::stoi(field[1])];
                item.first = field[0];
                item.second[current_pos] = total_spin(line);
            }
        }
    }
    return table;
}

void NBOSpin::process() {
    std::string filename = choose_file("Enter NBO file : ");
    do_process(filename);
}
