//
// Created by xiamr on 3/17/19.
//

#ifndef TINKER_MOLECULE_HPP
#define TINKER_MOLECULE_HPP

#include <boost/optional.hpp>

#include "utils/common.hpp"
#include "utils/std.hpp"

class Frame;

class Atom;

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                              boost::property<boost::vertex_distance_t, double>,
                              boost::property<boost::edge_weight_t, double>>
    Graph;

class Molecule {
public:
    double mass;  // molecular mass
    std::list<std::shared_ptr<Atom>> atom_list;
    boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;

    unsigned int sequence;

    bool bExculde = false;

    void calc_mass();

    std::tuple<double, double, double> calc_geom_center(const std::shared_ptr<Frame> &frame);

    std::tuple<double, double, double> calc_geom_center_inplace();

    std::tuple<double, double, double> calc_weigh_center(const std::shared_ptr<Frame> &frame,
                                                         bool includeHydrogen = true);

    std::tuple<double, double, double> calc_charge_center(const std::shared_ptr<Frame> &frame);

    std::tuple<double, double, double> calc_dipole(const std::shared_ptr<Frame> &frame);

    void do_aggregate(const std::shared_ptr<Frame> &frame);
    /**
     * @return return seq of first atom for temporary use
     *         else return 0
     */
    int seq();

    graph_t g;
    void build_graph(const std::shared_ptr<Frame> &frame);

    std::tuple<double, double, double> inplace_geometry_center;
};

std::pair<double, std::array<std::shared_ptr<Atom>, 2>> min_distance(const std::shared_ptr<Molecule> &mol1,
                                                                     const std::shared_ptr<Molecule> &mol2,
                                                                     const std::shared_ptr<Frame> &frame);

#endif  // TINKER_MOLECULE_HPP
