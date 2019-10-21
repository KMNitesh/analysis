//
// Created by xiamr on 10/18/19.
//

#ifndef TINKER_GQUADRUPLEXPDB2GMX_HPP
#define TINKER_GQUADRUPLEXPDB2GMX_HPP

#include "std.hpp"

class Frame;

class GQuadruplexPdb2gmx {
public:
    static std::string title() { return "Convert 3dnus G-Quadruplex pdb for gmx pdb2gmx"; }

    static void convert();

    static void superpose_and_move();
};


#endif //TINKER_GQUADRUPLEXPDB2GMX_HPP
