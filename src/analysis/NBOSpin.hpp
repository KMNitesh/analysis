//
// Created by xiamr on 11/1/19.
//

#ifndef TINKER_NBOSPIN_HPP
#define TINKER_NBOSPIN_HPP

#include "std.hpp"

class NBOSpin {
public:
    static std::string title() { return "NBO Spin"; }

    static void process();

    static void do_process(const std::string &filename);

    static double total_spin(std::string line);

    static std::map<int, std::pair<std::string, std::array<double, 3>>> getElectronSpin(std::istream &ifs);
};


#endif //TINKER_NBOSPIN_HPP
