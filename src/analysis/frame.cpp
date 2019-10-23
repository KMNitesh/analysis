//
// Created by xiamr on 3/17/19.
//
#include "config.h"
#include "frame.hpp"
#include "common.hpp"
#include "molecule.hpp"
#include "atom.hpp"


void Frame::image(double &xr, double &yr, double &zr) const {
    if (!enable_bound) return;
    while (std::abs(xr) > a_axis_half) xr -= sign(a_axis, xr);
    while (std::abs(yr) > b_axis_half) yr -= sign(b_axis, yr);
    while (std::abs(zr) > c_axis_half) zr -= sign(c_axis, zr);
}

void Frame::image(Eigen::Array3d &r) const {
    image(r[0], r[1], r[2]);
}

void Frame::image(std::tuple<double, double, double> &r) const {
    auto &[xr, yr, zr] = r;
    image(xr, yr, zr);
}

std::tuple<double, double, double> Frame::getDipole() {
    auto frame = shared_from_this();
    std::tuple<double, double, double> system_dipole{};
    for (auto &mol : this->molecule_list) {
        system_dipole += mol->calc_dipole(frame);
    }
    return system_dipole;
}


