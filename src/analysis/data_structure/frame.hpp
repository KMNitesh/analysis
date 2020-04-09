//
// Created by xiamr on 3/17/19.
//

#ifndef TINKER_FRAME_HPP
#define TINKER_FRAME_HPP

#include <Eigen/Eigen>
#include <cassert>
#include <list>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utils/PBCBox.hpp>

#include "config.h"
#include "utils/common.hpp"

namespace gmx {
#include <gromacs/utility/real.h>
}

class Atom;

class Molecule;

class Frame : public std::enable_shared_from_this<Frame> {
    std::optional<float> current_time;

public:
    const std::optional<float> &getCurrentTime() const { return current_time; }

    void setCurrentTime(const std::optional<float> &currentTime) { current_time = currentTime; }

    PBCBox box;
    std::string title;

    std::vector<std::shared_ptr<Atom>> atom_list;
    std::unordered_map<int, std::shared_ptr<Atom>> atom_map;

    std::vector<std::shared_ptr<Molecule>> molecule_list;

    bool enable_bound = false;

    bool has_velocity = false;

    void image(double &xr, double &yr, double &zr) const;

    void image(std::array<double, 3> &r) const;

    void image(std::tuple<double, double, double> &r) const;

    double volume() const {
        assert(enable_bound);
        return box.volume();
    }

    std::tuple<double, double, double> getDipole();

    void build_graph();

    struct harmonic {
        double krA, rA;
    };

    struct pdihs {
        double phiA, cpA;
        int mult;
    };

    std::map<std::array<std::shared_ptr<Atom>, 2>, harmonic> f_bond_params;
    std::map<std::array<std::shared_ptr<Atom>, 3>, harmonic> f_angle_params;
    std::multimap<std::array<std::shared_ptr<Atom>, 4>, pdihs> f_dihedral_params;
    std::multimap<std::array<std::shared_ptr<Atom>, 4>, pdihs> f_improper_dihedral_params;

    std::map<uint, std::shared_ptr<Atom>> real_residue_num_map;
};

#endif // TINKER_FRAME_HPP
