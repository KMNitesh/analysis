//
// Created by xiamr on 3/17/19.
//

#ifndef TINKER_FRAME_HPP
#define TINKER_FRAME_HPP

#include "config.h"
#include <string>
#include <list>
#include <unordered_map>
#include <memory>
#include <cassert>

class Atom;
class Molecule;

class Frame {
public:
    double a_axis = 0.0;
    double b_axis = 0.0;
    double c_axis = 0.0;

    double a_axis_half = 0.0;
    double b_axis_half = 0.0;
    double c_axis_half = 0.0;

    double alpha = 0.0;
    double beta = 0.0;
    double gamma = 0.0;

    std::string title;

    std::list<std::shared_ptr<Atom>> atom_list;
    std::unordered_map<int, std::shared_ptr<Atom>> atom_map;

    std::list<std::shared_ptr<Molecule>> molecule_list;

    bool enable_bound = false;

    void image(double &xr, double &yr, double &zr);

    double volume() const {
        assert(enable_bound);
        return a_axis * b_axis * c_axis;
    }
};

#endif //TINKER_FRAME_HPP
