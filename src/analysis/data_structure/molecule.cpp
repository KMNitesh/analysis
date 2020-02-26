//
// Created by xiamr on 3/17/19.
//

#include <boost/graph/breadth_first_search.hpp>
#include "config.h"
#include "utils/common.hpp"
#include "molecule.hpp"
#include "data_structure/atom.hpp"
#include "data_structure/frame.hpp"
#include "forcefield.hpp"
#include "utils/ThrowAssert.hpp"
#include "ana_module/HBond.hpp"

int Molecule::seq() {
    if (atom_list.empty()) return 0;
    return atom_list.front()->seq;
}

void Molecule::calc_mass() {
    double mol_mass = 0.0;
    for (auto &atom : atom_list) {
        assert(atom->mass);
        mol_mass += atom->mass.get();
    }
    mass = mol_mass;
}


namespace {
    class AggregateVisitor : public boost::default_bfs_visitor {
        public:
            explicit AggregateVisitor(const std::shared_ptr<Frame> &frame) : frame(frame) {};

            template<typename Edge, typename Graph>
                void tree_edge(Edge e, const Graph &g) const {
                    auto &source = g[boost::source(e, g)];
                    auto &target = g[boost::target(e, g)];
                    auto r = target->getCoordinate() - source->getCoordinate();
                    frame->image(r);
                    std::tie(target->x, target->y, target->z) = r + source->getCoordinate();
                }

        private:
            const std::shared_ptr<Frame> &frame;
    };
    
    template<typename Func>
    class CenterVisitor : public boost::default_bfs_visitor {
        public:
            explicit CenterVisitor(const std::shared_ptr<Frame> &frame, Func f) : frame(frame), func(std::move(f)) {};

            template<typename Edge, typename Graph>
                void tree_edge(Edge e, const Graph &g) {
                    const auto &source_loc = map[boost::source(e, g)];
                    const auto &target = g[boost::target(e, g)];
                    auto r = target->getCoordinate() - source_loc;
                    frame->image(r);
                    r += source_loc;
                    func(r, target);
                    map[boost::target(e, g)] = r;
                }
        private:
            const std::shared_ptr<Frame> &frame;
            std::map<boost::graph_traits<graph_t>::vertex_descriptor, std::tuple<double,double,double>> map;
            Func func;

    };

}

std::tuple<double, double, double> Molecule::calc_weigh_center(const std::shared_ptr<Frame> &frame, bool includeHydrogen) {
    if (boost::num_vertices(g) == 0) build_graph(frame);
    double mol_mass = atom_list.front()->mass.value();
    auto sum = atom_list.front()->getCoordinate()*mol_mass;

    CenterVisitor vis(frame, [&sum, &mol_mass, includeHydrogen](auto &r, const std::shared_ptr<Atom> &atom){
            if (!includeHydrogen and which(atom) == Symbol::Hydrogen) return;
            auto weight = atom->mass.value();
            mol_mass += weight;
            sum += weight * r;
    });
    boost::breadth_first_search(g, atom_list.front()->vertex_descriptor, boost::visitor(vis));
    return sum / mol_mass;
}

std::tuple<double, double, double> Molecule::calc_charge_center(const std::shared_ptr<Frame> &frame) {
    if (boost::num_vertices(g) == 0) build_graph(frame);
    double mol_charge = atom_list.front()->charge.value();
    std::tuple<double,double,double> sum = atom_list.front()->getCoordinate() * mol_charge;

    CenterVisitor vis(frame, [&sum, &mol_charge](auto &r, const std::shared_ptr<Atom> &atom){
            auto charge = atom->charge.get();
            mol_charge += charge;
            sum += charge * r;
    });
    boost::breadth_first_search(g, atom_list.front()->vertex_descriptor, boost::visitor(vis));
    
    if (std::abs(mol_charge) < 1E-3) {
        return calc_geom_center(frame);
    }
    return sum / mol_charge;
}

std::tuple<double, double, double> Molecule::calc_dipole(const std::shared_ptr<Frame> &frame) {
    if (boost::num_vertices(g) == 0) build_graph(frame);
    auto dipole = atom_list.front()->getCoordinate() * atom_list.front()->charge.value();

    CenterVisitor vis(frame, [&dipole](auto &r, const std::shared_ptr<Atom> &atom){
            auto charge = atom->charge.get();
            dipole += charge * r;
    });
    boost::breadth_first_search(g, atom_list.front()->vertex_descriptor, boost::visitor(vis));
    return dipole;
}

std::tuple<double,double,double> Molecule::calc_geom_center(const std::shared_ptr<Frame> &frame) {
    if (boost::num_vertices(g) == 0) build_graph(frame);
    auto sum = atom_list.front()->getCoordinate();
    CenterVisitor vis(frame, [&sum](auto &r, const auto&){ sum += r;});
    boost::breadth_first_search(g, atom_list.front()->vertex_descriptor, boost::visitor(vis));
    return sum / atom_list.size();
}


void Molecule::do_aggregate(std::shared_ptr<Frame> &frame){
    if (boost::num_vertices(g) == 0) build_graph(frame);
    boost::breadth_first_search(g, atom_list.front()->vertex_descriptor, boost::visitor(AggregateVisitor(frame)));

}

void Molecule::build_graph(const std::shared_ptr<Frame> &frame) {
    for (auto &atom : atom_list) {
        atom->vertex_descriptor = boost::add_vertex(atom, g);
    }
    for (auto &atom : atom_list) {
        for ( auto i : atom->con_list) {
            boost::add_edge(atom->vertex_descriptor, frame->atom_map[i]->vertex_descriptor, g);
        }
    }
}


double min_distance(std::shared_ptr<Molecule> &mol1, std::shared_ptr<Molecule> &mol2, std::shared_ptr<Frame> &frame) {
    double mindistance2 = std::numeric_limits<double>::max();
    for (auto &atom1 : mol1->atom_list) {
        for (auto &atom2 : mol2->atom_list) {
            mindistance2 = std::min(mindistance2, atom_distance2(atom1, atom2, frame));
        }
    }
    return std::sqrt(mindistance2);
}


