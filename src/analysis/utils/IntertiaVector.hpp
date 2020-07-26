
#ifndef TINKER_INTERTIAVECTOR_HPP
#define TINKER_INTERTIAVECTOR_HPP

#include "ana_module/Angle.hpp"
#include "dsl/AmberMask.hpp"

class IntertiaVector {
public:
    // IntertiaVector(const AmberMask &mask, Angle::AxisType axis) : mask(mask), axis(axis){}

    IntertiaVector(const AmberMask &mask, const std::string &axis);

    AmberMask mask;
    Angle::AxisType axis;
};

#endif // TINKER_INTERTIAVECTOR_HPP