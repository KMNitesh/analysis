//
// Created by xiamr on 3/17/19.
//

#include "config.h"
#include "common.hpp"
#include "molecule.hpp"
#include "atom.hpp"
#include "frame.hpp"
#include "forcefield.hpp"
#include "ThrowAssert.hpp"
#include "HBond.hpp"

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


std::tuple<double, double, double>
Molecule::calc_weigh_center(const std::shared_ptr<Frame> &frame, bool includeHydrogen) {
    double xmid = 0.0;
    double ymid = 0.0;
    double zmid = 0.0;
    bool first_atom = true;
    double first_x, first_y, first_z;
    double mol_mass = 0.0;
    for (auto &atom : atom_list) {
        if (!includeHydrogen and which(atom) == Symbol::Hydrogen) continue;

        if (first_atom) {
            first_atom = false;
            first_x = atom->x;
            first_y = atom->y;
            first_z = atom->z;
            assert(atom->mass);
            double weigh = atom->mass.get();
            mol_mass += weigh;
            xmid = first_x * weigh;
            ymid = first_y * weigh;
            zmid = first_z * weigh;

        } else {
            double xr = atom->x - first_x;
            double yr = atom->y - first_y;
            double zr = atom->z - first_z;
            frame->image(xr, yr, zr);
            assert(atom->mass);
            double weigh = atom->mass.get();
            mol_mass += weigh;
            xmid += (first_x + xr) * weigh;
            ymid += (first_y + yr) * weigh;
            zmid += (first_z + zr) * weigh;
        }
    }

    return std::make_tuple(xmid / mol_mass, ymid / mol_mass, zmid / mol_mass);
}

std::tuple<double, double, double> Molecule::calc_charge_center(const std::shared_ptr<Frame> &frame) {
    double xmid = 0.0;
    double ymid = 0.0;
    double zmid = 0.0;
    bool first_atom = true;
    double first_x, first_y, first_z;
    double mol_charge = 0.0;
    for (auto &atom : atom_list) {
        if (first_atom) {
            first_atom = false;
            first_x = atom->x;
            first_y = atom->y;
            first_z = atom->z;
            assert(atom->charge);
            double charge = atom->charge.value();
            mol_charge += charge;
            xmid = first_x * charge;
            ymid = first_y * charge;
            zmid = first_z * charge;

        } else {
            double xr = atom->x - first_x;
            double yr = atom->y - first_y;
            double zr = atom->z - first_z;
            frame->image(xr, yr, zr);
            assert(atom->charge);
            double charge = atom->charge.value();
            mol_charge += charge;
            xmid += (first_x + xr) * charge;
            ymid += (first_y + yr) * charge;
            zmid += (first_z + zr) * charge;
        }
    }
    if (abs(mol_charge) < 1E-3) {
        calc_geom_center(frame);
        return {center_x, center_y, center_z};
    }

    return {xmid / mol_charge, ymid / mol_charge, zmid / mol_charge};
}

std::tuple<double, double, double> Molecule::calc_dipole(const std::shared_ptr<Frame> &frame) {

    double dipole_x = 0.0;
    double dipole_y = 0.0;
    double dipole_z = 0.0;

    bool first_atom = true;

    double first_x, first_y, first_z;

    for (auto &atom : atom_list) {
        if (first_atom) {
            first_x = atom->x;
            first_y = atom->y;
            first_z = atom->z;
            first_atom = false;
        } else {
            double xr = atom->x - first_x;
            double yr = atom->y - first_y;
            double zr = atom->z - first_z;

            frame->image(xr, yr, zr);

            dipole_x += atom->charge.value() * xr;
            dipole_y += atom->charge.value() * yr;
            dipole_z += atom->charge.value() * zr;
        }
    }
    return {dipole_x, dipole_y, dipole_z};
}

void Molecule::calc_geom_center(const std::shared_ptr<Frame> &frame) {
    double sum_x = 0.0;
    double sum_y = 0.0;
    double sum_z = 0.0;
    bool first_atom = true;
    double first_x, first_y, first_z;
    for (auto &atom : atom_list) {
        if (first_atom) {
            first_atom = false;
            sum_x = first_x = atom->x;
            sum_y = first_y = atom->y;
            sum_z = first_z = atom->z;
        } else {
            double xr = atom->x - first_x;
            double yr = atom->y - first_y;
            double zr = atom->z - first_z;
            frame->image(xr, yr, zr);
            sum_x += first_x + xr;
            sum_y += first_y + yr;
            sum_z += first_z + zr;
        }
    }
    auto len = atom_list.size();
    center_x = sum_x / len;
    center_y = sum_y / len;
    center_z = sum_z / len;
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


