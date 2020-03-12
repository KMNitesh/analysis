//
// Created by xiamr on 11/1/19.
//

#include "NBOSpin.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/xpressive/regex_actions.hpp>
#include <boost/xpressive/xpressive.hpp>

#include "utils/common.hpp"

double NBOSpin::total_spin(std::string_view line) {
    using namespace boost::xpressive;
    placeholder<double> _result;
    double total = 0.0;
    static cregex rex =
        bos >> +(*~(set = '(') >> '(' >> *blank >> (+~(set = ')'))[_result += as<double>(_)] >> ')') >> eos;
    cmatch what;
    what.let(_result = total);

    if (!regex_match(line, what, rex)) throw std::runtime_error("Regex not match");
    return total;
}

void NBOSpin::do_process(const std::string &filename) {
    std::ifstream ifs(filename);
    if (!ifs) {
        std::cerr << "ERROR!! Cannot open file " << filename << '\n';
        exit(EXIT_FAILURE);
    }

    std::map<int, std::pair<std::string, std::array<double, 3>>> table = getElectronSpin(ifs);

    std::cout << boost::format("%6s %4s %6s %6s %6s %6s\n") % "Atom" % "No" % "Total" % "Alpha" % "Beta" % "Alpha-Beta";

    const auto fmt = boost::format("%6s %4d %6.2f %6.2f %6.2f %6.2f  %s\n");
    for (auto &[atomNo, content] : table) {
        auto &[name, spins] = content;
        auto diff = spins[1] - spins[2];
        std::cout << boost::format(fmt) % name % atomNo % spins[0] % spins[1] % spins[2] % diff %
                         (std::abs(diff) >= 0.1 ? '*' : ' ');
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
                    if (++current_pos == 3) return table;
                    break;
                }
                auto field = split(line);
                auto &[name, spins] = table[std::stoi(field[1])];
                if (current_pos == 0) name = field[0];
                spins[current_pos] = total_spin(line);
            }
        }
    }
    return table;
}

void NBOSpin::process() {
    std::string filename = choose_file("Enter NBO file : ").isExist(true);
    do_process(filename);
}
