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

    auto &[a_axis, b_axis, c_axis] = frame->box.getAxis();
    for (auto &mol : frame->molecule_list) {
        mol->calc_geom_center(frame);
        double center_x_origin = mol->center_x;
        double center_y_origin = mol->center_y;
        double center_z_origin = mol->center_z;
        while (std::abs(mol->center_x) > 0.5 * a_axis) {
            mol->center_x -= sign(a_axis, mol->center_x);
        }
        while (std::abs(mol->center_y) > 0.5 * b_axis) {
            mol->center_y -= sign(b_axis, mol->center_y);
        }
        while (std::abs(mol->center_z) > 0.5 * c_axis) {
            mol->center_z -= sign(c_axis, mol->center_z);
        }
        if (frame->box.get_box_type() == PBCBox::Type::octahedron) {
            if (std::abs(mol->center_x) + std::abs(mol->center_y) + std::abs(mol->center_z) > 0.75 * a_axis) {
                mol->center_x -= sign(0.5 * a_axis, mol->center_x);
                mol->center_y -= sign(0.5 * b_axis, mol->center_y);
                mol->center_z -= sign(0.5 * c_axis, mol->center_z);
            }
        }
        auto x_move = mol->center_x - center_x_origin;
        auto y_move = mol->center_y - center_y_origin;
        auto z_move = mol->center_z - center_z_origin;
        for (auto &atom : mol->atom_list) {
            atom->x += x_move + 0.5 * a_axis;
            atom->y += y_move + 0.5 * b_axis;
            atom->z += z_move + 0.5 * c_axis;
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
