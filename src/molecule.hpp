//
// Created by xiamr on 3/17/19.
//

#ifndef TINKER_MOLECULE_HPP
#define TINKER_MOLECULE_HPP

#include "config.h"
#include <list>
#include <memory>
#include <tuple>
#include <boost/optional.hpp>

class Frame;

class Atom;

class Molecule {
public:
    double mass;  // molecular mass
    std::list<std::shared_ptr<Atom>> atom_list;

    boost::optional<int> sequence;

    double center_x, center_y, center_z;

    bool bExculde = false;

    void calc_geom_center(const std::shared_ptr<Frame> &frame);

    void calc_mass();

    std::tuple<double, double, double> calc_weigh_center(const std::shared_ptr<Frame> &frame);

    std::tuple<double, double, double> calc_charge_center(const std::shared_ptr<Frame> &frame);

    std::tuple<double, double, double> calc_dipole(const std::shared_ptr<Frame> &frame);

    /**
     * @return return seq of first atom for temporary use
     *         else return 0
     */
    int seq();
};


double min_distance(std::shared_ptr<Molecule> &mol1, std::shared_ptr<Molecule> &mol2, std::shared_ptr<Frame> &frame);

#endif //TINKER_MOLECULE_HPP
