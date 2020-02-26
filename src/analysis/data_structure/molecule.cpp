//
// Created by xiamr on 3/17/19.
//

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


std::tuple<double, double, double> Molecule::calc_weigh_center(const std::shared_ptr<Frame> &frame, bool includeHydrogen) {
    std::tuple<double,double,double> sum;
    double mol_mass = 0.0;
    auto pre = atom_list.front()->getCoordinate();
    for (auto &atom : atom_list) {
        if (!includeHydrogen and which(atom) == Symbol::Hydrogen) continue;
        auto r = atom->getCoordinate() - pre;
        frame->image(r);
        pre += r;
        auto weight = atom->mass.value();
        mol_mass += weight;
        sum += mol_mass * r;
    }
    return sum / mol_mass;
}

std::tuple<double, double, double> Molecule::calc_charge_center(const std::shared_ptr<Frame> &frame) {
    std::tuple<double,double,double> sum{};
    double mol_charge{};
    auto pre = atom_list.front()->getCoordinate();

    for (auto &atom : atom_list) {
        auto r = atom->getCoordinate() - pre;
        frame->image(r);
        pre += r;
        auto charge = *(atom->charge);
        mol_charge += charge;
        sum += charge * pre;
    }
    if (abs(mol_charge) < 1E-3) {
        return calc_geom_center(frame);
    }
    return sum / mol_charge;
}

std::tuple<double, double, double> Molecule::calc_dipole(const std::shared_ptr<Frame> &frame) {
    std::tuple<double,double,double> dipole{};
    auto pre = atom_list.front()->getCoordinate();
    for (auto &atom : atom_list) {
        auto r = atom->getCoordinate() - pre;
        frame->image(r);
        pre += r;
        dipole += atom->charge.value() * pre;
    }
    return dipole;
}

std::tuple<double,double,double> Molecule::calc_geom_center(const std::shared_ptr<Frame> &frame) {
    std::tuple<double,double,double> sum{};
    auto pre_coord = atom_list.front()->getCoordinate();
    for (auto &atom : atom_list) {
        auto r = atom->getCoordinate() - pre_coord;
        frame->image(r);
        pre_coord += r;
        sum += pre_coord;
    }
    return sum / atom_list.size();
}


void Molecule::do_aggregate(std::shared_ptr<Frame> &frame){
    auto pre_coord = atom_list.front()->getCoordinate();
    for (auto &atom : atom_list){
        auto r = atom->getCoordinate() - pre_coord;
        frame->image(r);
        pre_coord += r;
        std::tie(atom->x, atom->y, atom->z) = pre_coord;
    }
}

void Molecule::build_graph(std::shared_ptr<Frame> &frame) {
    for(auto &atom : atom_list) {
        for ( auto i : atom->con_list) {
            boost::add_edge(atom->vertex_descriptor, frame->atom_map[i]->vertex_descriptor, g);
        }
    }
}

template<typename Edge, typename Graph>
void Molecule::AggregateVisitor::tree_edge(Edge e, const Graph &g) const {
    auto &source = frame->vertex_descriptor_map[boost::source(e, g)];
    auto &target = frame->vertex_descriptor_map[boost::source(e, g)];
    auto r = target->getCoordinate() - source->getCoordinate();
    frame->image(r);
    r += source->getCoordinate();
    std::tie(target->x, target->y, target->z) = r;
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


