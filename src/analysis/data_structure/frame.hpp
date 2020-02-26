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
#include <tuple>
#include <cassert>
#include <Eigen/Eigen>
#include <utils/PBCBox.hpp>
#include "utils/common.hpp"

class Atom;

class Molecule;

class Frame : public std::enable_shared_from_this<Frame> {

    std::optional<float> current_time;
public:
    const std::optional<float> &getCurrentTime() const {
        return current_time;
    }

    void setCurrentTime(const std::optional<float> &currentTime) {
        current_time = currentTime;
    }

    PBCBox box;
    std::string title;

    std::list<std::shared_ptr<Atom>> atom_list;
    std::unordered_map<int, std::shared_ptr<Atom>> atom_map;
    std::unordered_map<boost::graph_traits<graph_t>::vertex_descriptor, std::shared_ptr<Atom>> vertex_descriptor_map;

    std::list<std::shared_ptr<Molecule>> molecule_list;

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
};

#endif //TINKER_FRAME_HPP
