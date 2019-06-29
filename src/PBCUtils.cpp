//
// Created by xiamr on 6/29/19.
//

#include "PBCUtils.hpp"
#include "frame.hpp"
#include "atom.hpp"
#include "molecule.hpp"

void PBCUtils::do_move_center_basedto_atom(int num, std::shared_ptr<Frame> &frame) const {
    auto center_x = frame->atom_map[num]->x;
    auto center_y = frame->atom_map[num]->y;
    auto center_z = frame->atom_map[num]->z;

    for (auto &mol : frame->molecule_list) {
        for (auto &atom : mol->atom_list) {
            atom->x -= center_x;
            atom->y -= center_y;
            atom->z -= center_z;
        }
    }
}

void PBCUtils::do_move_center_basedto_molecule(int num, std::shared_ptr<Frame> &frame) const {
    auto center_mol = frame->atom_map[num]->molecule.lock();
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

void PBCUtils::doPBC(Trajconv::PBCType pbc_mode, int num, std::shared_ptr<Frame> &frame) const {
    if (pbc_mode == Trajconv::PBCType::OneAtom) {
        do_move_center_basedto_atom(num, frame);
    } else if (pbc_mode == Trajconv::PBCType::OneMol) {
        do_move_center_basedto_molecule(num, frame);
    }
    if (pbc_mode != Trajconv::PBCType::None) {
        do_molecule_aggregate(frame);
    }
}
