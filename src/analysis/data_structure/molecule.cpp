//
// Created by xiamr on 3/17/19.
//

#include "molecule.hpp"

#include <boost/graph/breadth_first_search.hpp>
#include <boost/range/numeric.hpp>

#include "ana_module/HBond.hpp"
#include "config.h"
#include "data_structure/atom.hpp"
#include "data_structure/forcefield.hpp"
#include "data_structure/frame.hpp"
#include "utils/PBCUtils.hpp"
#include "utils/ThrowAssert.hpp"
#include "utils/common.hpp"

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

template <typename Func>
class CenterVisitor : public boost::default_bfs_visitor {
   public:
    explicit CenterVisitor(const std::shared_ptr<Frame> &frame, std::shared_ptr<Atom> &atom, Func f)
        : frame(frame), func(std::move(f)) {
        map[atom->vertex_descriptor] = atom->getCoordinate();
    };

    template <typename Edge, typename Graph>
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
    std::map<boost::graph_traits<graph_t>::vertex_descriptor, std::tuple<double, double, double>> map;
    Func func;
};

}  // namespace

std::tuple<double, double, double> Molecule::calc_weigh_center(const std::shared_ptr<Frame> &frame,
                                                               bool includeHydrogen) {
    double mol_mass = atom_list.front()->mass.value();
    auto sum = atom_list.front()->getCoordinate() * mol_mass;

    CenterVisitor vis(frame, atom_list.front(),
                      [&sum, &mol_mass, includeHydrogen](auto &r, const std::shared_ptr<Atom> &atom) {
                          if (!includeHydrogen and which(atom) == Symbol::Hydrogen) return;
                          auto weight = atom->mass.value();
                          mol_mass += weight;
                          sum += weight * r;
                      });
    boost::breadth_first_search(g, atom_list.front()->vertex_descriptor, boost::visitor(vis));
    return sum / mol_mass;
}

std::tuple<double, double, double> Molecule::calc_charge_center(const std::shared_ptr<Frame> &frame) {
    double mol_charge = atom_list.front()->charge.value();
    std::tuple<double, double, double> sum = atom_list.front()->getCoordinate() * mol_charge;

    CenterVisitor vis(frame, atom_list.front(), [&sum, &mol_charge](auto &r, const std::shared_ptr<Atom> &atom) {
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
    auto dipole = atom_list.front()->getCoordinate() * atom_list.front()->charge.value();

    CenterVisitor vis(frame, atom_list.front(), [&dipole](auto &r, const std::shared_ptr<Atom> &atom) {
        auto charge = atom->charge.get();
        dipole += charge * r;
    });
    boost::breadth_first_search(g, atom_list.front()->vertex_descriptor, boost::visitor(vis));
    return dipole;
}

std::tuple<double, double, double> Molecule::calc_geom_center(const std::shared_ptr<Frame> &frame) {
    auto sum = atom_list.front()->getCoordinate();
    CenterVisitor vis(frame, atom_list.front(), [&sum](auto &r, const auto &) { sum += r; });
    boost::breadth_first_search(g, atom_list.front()->vertex_descriptor, boost::visitor(vis));
    return sum / atom_list.size();
}

std::tuple<double, double, double> Molecule::calc_geom_center_inplace() {
    inplace_geometry_center =
        boost::accumulate(atom_list, std::tuple<double, double, double>{},
                          [](auto sum, const auto &atom) { return sum + atom->getCoordinate(); }) /
        atom_list.size();

    return inplace_geometry_center;
}

void Molecule::do_aggregate(const std::shared_ptr<Frame> &frame) {
    boost::breadth_first_search(g, atom_list.front()->vertex_descriptor, boost::visitor(AggregateVisitor(frame)));
}

std::pair<double, std::array<std::shared_ptr<Atom>, 2>> min_distance(const std::shared_ptr<Molecule> &mol1,
                                                                     const std::shared_ptr<Molecule> &mol2,
                                                                     const std::shared_ptr<Frame> &frame) {
    double mindistance2 = std::numeric_limits<double>::max();
    std::array<std::shared_ptr<Atom>, 2> atom_pair;
    for (auto &atom1 : mol1->atom_list) {
        for (auto &atom2 : mol2->atom_list) {
            auto d2 = atom_distance2(atom1, atom2, frame);
            if (d2 < mindistance2) {
                mindistance2 = d2;
                atom_pair[0] = atom1;
                atom_pair[1] = atom2;
            }
        }
    }
    return {std::sqrt(mindistance2), atom_pair};
}

void Molecule::build_graph(const std::shared_ptr<Frame> &frame) {
    for (auto &atom : atom_list) {
        atom->vertex_descriptor = boost::add_vertex(atom, g);
    }
    for (auto &atom : atom_list) {
        for (auto i : atom->con_list) {
            boost::add_edge(atom->vertex_descriptor, frame->atom_map[i]->vertex_descriptor, g);
        }
    }
}
