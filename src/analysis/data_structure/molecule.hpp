//
// Created by xiamr on 3/17/19.
//

#ifndef TINKER_MOLECULE_HPP
#define TINKER_MOLECULE_HPP

#include "utils/std.hpp"
#include <boost/optional.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>

class Frame;

class Atom;

class Molecule {
public:
    double mass;  // molecular mass
    std::list<std::shared_ptr<Atom>> atom_list;

    unsigned int sequence;

    bool bExculde = false;

    void calc_mass();

    std::tuple<double, double, double> calc_geom_center(const std::shared_ptr<Frame> &frame);

    std::tuple<double, double, double> calc_weigh_center(const std::shared_ptr<Frame> &frame, bool includeHydrogen = true);

    std::tuple<double, double, double> calc_charge_center(const std::shared_ptr<Frame> &frame);

    std::tuple<double, double, double> calc_dipole(const std::shared_ptr<Frame> &frame);

    void do_aggregate(std::shared_ptr<Frame> &frame);
    /**
     * @return return seq of first atom for temporary use
     *         else return 0
     */
    int seq();

    void build_graph(std::shared_ptr<Frame> &frame);

private:
    graph_t g;

    boost::graph_traits<graph_t>::vertex_descriptor vertex_descriptor;

    class AggregateVisitor : boost::default_bfs_visitor {
    public:
        explicit AggregateVisitor(std::shared_ptr<Frame> &frame) : frame(frame) {};

        template<typename Edge, typename Graph>
        void tree_edge(Edge e, const Graph &g) const;

    private:
        std::shared_ptr<Frame> &frame;
    };
};


double min_distance(std::shared_ptr<Molecule> &mol1, std::shared_ptr<Molecule> &mol2, std::shared_ptr<Frame> &frame);

#endif //TINKER_MOLECULE_HPP
