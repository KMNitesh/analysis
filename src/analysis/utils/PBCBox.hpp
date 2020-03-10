#ifndef TINKER_PBCBOX_HPP
#define TINKER_PBCBOX_HPP

#include <array>

namespace gmx {

#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/pbc.h"

}

class PBCUtils;

class PBCBox {
public:
 enum class Type { orthogonal, octahedron, other };

    PBCBox() = default;

    explicit PBCBox(double xbox, double ybox, double zbox, double alpha, double beta, double gamma);

    explicit PBCBox(gmx::matrix box);

    const std::array<double, 3> &getAxis() const { return axis; }

    std::array<double, 6> getBoxParameter() const {
        return {axis[0], axis[1], axis[2], angles[0], angles[1], angles[2]};
    }

    void getBoxParameter(gmx::matrix box) const;

    void image(double &xr, double &yr, double &zr) const;

    double volume() const;

    Type get_box_type() { return box_type; }

private:

    bool numeric_near(double value, double target, double difference = 0.01) const {
        return std::abs(value - target) <= difference;
    }

    void check_box_type();

    Type box_type;
    std::array<double, 3> axis;
    std::array<double, 3> angles;

    friend class PBCUtils;
    mutable gmx::matrix box;
    mutable gmx::t_pbc pbc;
    mutable bool dirty = true;
};


#endif //TINKER_PBCBOX_HPP
