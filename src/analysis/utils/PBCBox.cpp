#include <cmath>
#include <boost/format.hpp>
#include "PBCBox.hpp"
#include "common.hpp"

void PBCBox::check_box_type() {
    if (numeric_near(angles[0], 90.0) and numeric_near(angles[1], 90.0) and numeric_near(angles[2], 90.0)) {
        box_type = Type::orthogonal;
    } else if (numeric_near(axis[0], axis[1], 1E-2) and numeric_near(axis[0], axis[2], 1E-2) and
            numeric_near(angles[0], 70.53, 0.1) and numeric_near(angles[1], 109.47, 0.1) and
            numeric_near(angles[2], 70.53, 0.1)) {
        box_type = Type::octahedron;
    } else {
        throw std::runtime_error((boost::format("unknown box type : %f %f %f %f %f %f\n")
                    % axis[0] % axis[1] % axis[2] % angles[0] % angles[1] % angles[2]).str());
    }
}

PBCBox::PBCBox(double xbox, double ybox, double zbox, double alpha, double beta, double gamma) :
    axis{xbox, ybox, zbox},
    angles{alpha, beta, gamma} {
        check_box_type();
    }

PBCBox::PBCBox(gmx::matrix box) {
    const auto &[v1x, v1y, v1z] = box[0];
    const auto &[v2x, v2y, v2z] = box[1];
    const auto &[v3x, v3y, v3z] = box[2];

    auto d1 = std::sqrt(v1x * v1x + v1y * v1y + v1z * v1z);
    auto d2 = std::sqrt(v2x * v2x + v2y * v2y + v2z * v2z);
    auto d3 = std::sqrt(v3x * v3x + v3y * v3y + v3z * v3z);

    axis[0] = 10 * d1;
    axis[1] = 10 * d2;
    axis[2] = 10 * d3;

    angles[0] = radian * std::acos((v2x * v3x + v2y * v3y + v2z * v3z) / (d2 * d3));
    angles[1] = radian * std::acos((v1x * v3x + v1y * v3y + v1z * v3z) / (d1 * d3));
    angles[2] = radian * std::acos((v1x * v2x + v1y * v2y + v1z * v2z) / (d1 * d2));

    check_box_type();

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j ) {
            this->box[i][j] = 10 * box[i][j];
        }
    }
    gmx::set_pbc(&pbc,-1,this->box);
    dirty = false;
}

void PBCBox::getBoxParameter(gmx::matrix box) const {

    if (dirty) {
        box[0][0] = axis[0];
        box[0][1] = 0.0;
        box[0][2] = 0.0;

        box[1][0] = axis[1] * std::cos(angles[2] / radian);
        box[1][1] = axis[1] * std::sin(angles[2] / radian);
        box[1][2] = 0.0;

        box[2][0] = axis[2] * std::cos(angles[1] / radian);
        box[2][1] = (axis[1] * axis[2] * std::cos(angles[0] / radian) - box[1][0] * box[2][0]) / box[1][1];
        box[2][2] = std::sqrt(axis[2] * axis[2] - box[2][0] * box[2][0] - box[2][1] * box[2][1]);

        box[0][0] *= 0.1;
        box[1][0] *= 0.1;
        box[1][1] *= 0.1;
        box[2][0] *= 0.1;
        box[2][1] *= 0.1;
        box[2][2] *= 0.1;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j ) {
                this->box[i][j] = 10 * box[i][j];
            }
        }
        gmx::set_pbc(&pbc,-1,this->box);
        dirty = false;
    } else {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j ) {
                box[i][j] = 0.1 * this->box[i][j];
            }
        }
    }
}

void PBCBox::image(double &xr, double &yr, double &zr) const {

    if (dirty) {
        gmx::matrix box;
        getBoxParameter(box);
    }

    if (box_type == Type::orthogonal) {
        while (std::abs(xr) > 0.5 * axis[0]) xr -= sign(axis[0], xr);
        while (std::abs(yr) > 0.5 * axis[1]) yr -= sign(axis[1], yr);
        while (std::abs(zr) > 0.5 * axis[2]) zr -= sign(axis[2], zr);
    } if (box_type == Type::octahedron) {
        gmx::rvec x1 = { xr , yr, zr };
        gmx::rvec x2 = { 0,0,0}; 
        gmx::rvec dx;
        gmx::pbc_dx(&pbc,x1,x2,dx);
        xr = dx[0];
        yr = dx[1];
        zr = dx[2];
    }

}

double PBCBox::volume() const {
    if (box_type == Type::orthogonal) return axis[0] * axis[1] * axis[2];

    const auto factor = 4.0 * std::sqrt(3) / 9.0;
    return factor * std::pow(axis[0], 3);
}
