//
// Created by xiamr on 3/17/19.
//
#include "frame.hpp"

#include "config.h"
#include "data_structure/atom.hpp"
#include "molecule.hpp"
#include "utils/common.hpp"

void Frame::image(double &xr, double &yr, double &zr) const {
    if (!enable_bound) return;
    box.image(xr, yr, zr);
}

void Frame::image(std::array<double, 3> &r) const { image(r[0], r[1], r[2]); }

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

void Frame::build_graph() {
    for (auto &mol : molecule_list) mol->build_graph(shared_from_this());
}
