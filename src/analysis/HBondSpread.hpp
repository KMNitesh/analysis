//
// Created by xiamr on 8/11/19.
//

#ifndef TINKER_HBONDSPREAD_HPP
#define TINKER_HBONDSPREAD_HPP

#include "std.hpp"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>

#include "BasicAnalysis.hpp"
#include "atom.hpp"


class Frame;

struct VertexNode {
    typedef boost::vertex_property_tag kind;
};

class MyVisitor : public boost::default_dfs_visitor {

public:
    MyVisitor(std::unordered_set<int> &all_connected_atoms) : all_connected_atoms(all_connected_atoms) {};

    template<typename Vertex, typename Graph>
    void discover_vertex(Vertex v, const Graph &g) const {
        all_connected_atoms.emplace(boost::get(VertexNode(), g, v));
    }

private:

    std::unordered_set<int> &all_connected_atoms;
};

class HBondSpread : public BasicAnalysis {
public:
    HBondSpread();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    static const std::string title() { return "Hydrogen Bond Spread Analysis"; }

protected:

    AmberMask center_Metal_atom_mask;
    AmberMask Ow_atom_mask;

    double dist_R_cutoff;
    double angle_HOO_cutoff;

    double cutoff2;


    typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, boost::property<VertexNode, int>> graph_t;
    graph_t g;

    std::unordered_set<std::shared_ptr<Atom>> metal;
    std::vector<std::shared_ptr<Atom>> Ow;
    std::vector<std::shared_ptr<Atom>> hydrogens;

    std::unordered_map<int, boost::graph_traits<graph_t>::vertex_descriptor> Ow_mapping;


    std::deque<std::pair<int, double>> hbond_connected_num_and_max_distance;
};


#endif //TINKER_HBONDSPREAD_HPP
