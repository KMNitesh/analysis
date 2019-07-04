//
// Created by xiamr on 3/17/19.
//
#include "config.h"
#include "frame.hpp"
#include "common.hpp"


void Frame::image(double &xr, double &yr, double &zr) {
    if (!enable_bound) return;
    while (std::abs(xr) > a_axis_half) xr -= sign(a_axis, xr);
    while (std::abs(yr) > b_axis_half) yr -= sign(b_axis, yr);
    while (std::abs(zr) > c_axis_half) zr -= sign(c_axis, zr);
}

void Frame::image(std::tuple<double, double, double> &r) {
    auto &[xr, yr, zr] = r;
    image(xr, yr, zr);
}
