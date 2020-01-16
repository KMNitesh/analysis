//
// Created by xiamr on 6/29/19.
//

#include "PBCUtils.hpp"
#include "data_structure/frame.hpp"
#include "data_structure/atom.hpp"
#include "data_structure/molecule.hpp"


std::shared_ptr<Atom> PBCUtils::find_atom(AmberMask &mask, std::shared_ptr<Frame> &frame) {
    std::shared_ptr<Atom> ret;
    for (auto &atom : frame->atom_list) {
        if (Atom::is_match(atom, mask)) {
            if (ret) {
                throw std::runtime_error("More then one atom selected");
            }
            ret = atom;
        }
    }
    return ret;
}

std::shared_ptr<Molecule> PBCUtils::find_molecule(AmberMask &mask, std::shared_ptr<Frame> &frame) {
    std::shared_ptr<Molecule> ret;
    for (auto &atom : frame->atom_list) {
        if (Atom::is_match(atom, mask)) {
            if (ret and atom->molecule.lock()) {
                throw std::runtime_error("More then one moleule selected");
            }
            ret = atom->molecule.lock();
        }
    }
    return ret;
}


void PBCUtils::do_move_center_basedto_atom(AmberMask &mask, std::shared_ptr<Frame> &frame) const {
    auto selected_atom = find_atom(mask, frame);
    auto center_x = selected_atom->x;
    auto center_y = selected_atom->y;
    auto center_z = selected_atom->z;

    for (auto &mol : frame->molecule_list) {
        for (auto &atom : mol->atom_list) {
            atom->x -= center_x;
            atom->y -= center_y;
            atom->z -= center_z;
        }
    }
}

void PBCUtils::do_move_center_basedto_molecule(AmberMask &mask, std::shared_ptr<Frame> &frame) const {
    auto center_mol = find_molecule(mask, frame);
    center_mol->calc_geom_center(frame);
    auto center_x = center_mol->center_x;
    auto center_y = center_mol->center_y;
    auto center_z = center_mol->center_z;
    for (auto &mol : frame->molecule_list) {
        for (auto &atom : mol->atom_list) {
            atom->x -= center_x;
            atom->y -= center_y;
            atom->z -= center_z;
        }
    }
}

void PBCUtils::do_molecule_aggregate(std::shared_ptr<Frame> &frame) const {
    for (auto &mol : frame->molecule_list) {
        mol->calc_geom_center(frame);
        double x_move = 0.0;
        double y_move = 0.0;
        double z_move = 0.0;
        while (mol->center_x > frame->a_axis_half) {
            mol->center_x -= frame->a_axis;
            x_move -= frame->a_axis;
        }
        while (mol->center_x < -frame->a_axis_half) {
            mol->center_x += frame->a_axis;
            x_move += frame->a_axis;
        }
        while (mol->center_y > frame->b_axis_half) {
            mol->center_y -= frame->b_axis;
            y_move -= frame->b_axis;
        }
        while (mol->center_y < -frame->b_axis_half) {
            mol->center_y += frame->b_axis;
            y_move += frame->b_axis;
        }
        while (mol->center_z > frame->c_axis_half) {
            mol->center_z -= frame->c_axis;
            z_move -= frame->c_axis;
        }
        while (mol->center_z < -frame->c_axis_half) {
            mol->center_z += frame->c_axis;
            z_move += frame->c_axis;
        }
        for (auto &atom : mol->atom_list) {
            atom->x += x_move + frame->a_axis_half;
            atom->y += y_move + frame->b_axis_half;
            atom->z += z_move + frame->c_axis_half;
        }
    }
}

void PBCUtils::doPBC(Trajconv::PBCType pbc_mode, AmberMask &mask, std::shared_ptr<Frame> &frame) const {
    if (pbc_mode == Trajconv::PBCType::OneAtom) {
        do_move_center_basedto_atom(mask, frame);
    } else if (pbc_mode == Trajconv::PBCType::OneMol) {
        do_move_center_basedto_molecule(mask, frame);
    }
    if (pbc_mode != Trajconv::PBCType::None) {
        do_molecule_aggregate(frame);
    }
}
