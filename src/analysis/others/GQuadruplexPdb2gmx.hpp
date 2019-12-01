//
// Created by xiamr on 10/18/19.
//

#ifndef TINKER_GQUADRUPLEXPDB2GMX_HPP
#define TINKER_GQUADRUPLEXPDB2GMX_HPP

#include "utils/std.hpp"

class Frame;

class GQuadruplexPdb2gmx {
public:

    [[nodiscard]] static std::string_view title() { return "Convert 3dnus G-Quadruplex pdb for gmx pdb2gmx"; }

    static void convert();

    static void superpose_and_move();

    static void renumberAtomAndResidueNum();
};


#endif //TINKER_GQUADRUPLEXPDB2GMX_HPP
