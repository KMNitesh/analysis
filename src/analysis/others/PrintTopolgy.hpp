//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_PRINTTOPOLGY_HPP
#define TINKER_PRINTTOPOLGY_HPP

#include <string>

class PrintTopolgy {
public:
    static void action(const std::string &topology_filename);

    [[nodiscard]] static std::string_view title() { return "Print Selected Atoms in Topolgoy File"; }
};

#endif  // TINKER_PRINTTOPOLGY_HPP
