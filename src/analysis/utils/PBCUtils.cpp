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

std::vector<std::shared_ptr<Atom>> PBCUtils::find_atoms(AmberMask &mask, std::shared_ptr<Frame> &frame) {
    std::vector<std::shared_ptr<Atom>> ret;
    for (auto &atom : frame->atom_list) {
        if (Atom::is_match(atom, mask)) {
            ret.push_back(atom);
        }
    }
    return ret;
}

std::shared_ptr<Molecule> PBCUtils::find_molecule(AmberMask &mask, std::shared_ptr<Frame> &frame) {
    std::shared_ptr<Molecule> ret;
    for (auto &atom : frame->atom_list) {
        if (Atom::is_match(atom, mask)) {
            if (ret and ret != atom->molecule.lock()) {
                throw std::runtime_error("More then one moleule selected");
            }
            ret = atom->molecule.lock();
        }
    }
    return ret;
}


void PBCUtils::do_move_center_basedto_atom(AmberMask &mask, std::shared_ptr<Frame> &frame) const {
    do_molecule_aggregate(frame);
    auto selected_atom = find_atom(mask, frame);
    auto center = selected_atom->getCoordinate();
    for (auto &mol : frame->molecule_list) {
        auto r = mol->calc_geom_center(frame) - center;
        auto old = r;
        frame->image(r);
        r -= old;
        for (auto &atom : mol->atom_list) {
            std::tie(atom->x, atom->y, atom->z) += r;
        }
    } 
    for (auto &atom : frame->atom_list){
        std::tie(atom->x, atom->y, atom->z) -= center - 
            0.5 * ( std::make_tuple(frame->box.box[0][0], frame->box.box[0][1], frame->box.box[0][2]) 
                    + std::make_tuple(frame->box.box[1][0], frame->box.box[1][1], frame->box.box[1][2])
                    + std::make_tuple(frame->box.box[2][0], frame->box.box[2][1], frame->box.box[2][2]));
    }
}

void PBCUtils::do_move_center_basedto_atom_group(AmberMask &mask, std::shared_ptr<Frame> &frame) const {
    do_molecule_aggregate(frame);
    auto selected_atoms = find_atoms(mask, frame);
    auto pre = selected_atoms.at(0)->getCoordinate();
    std::tuple<double,double,double> sum{};

    for(auto &atom : selected_atoms) {
        auto r = atom->getCoordinate() - pre;
        frame->image(r);
        pre += r;
        sum += pre;
    }

    auto center = sum / selected_atoms.size();

    for (auto &mol : frame->molecule_list) {
        auto r = mol->calc_geom_center(frame) - center;
        auto old = r;
        frame->image(r);
        r -= old;
        for (auto &atom : mol->atom_list) {
            std::tie(atom->x, atom->y, atom->z) += r;
        }
    }
    for (auto &atom : frame->atom_list){
        std::tie(atom->x, atom->y, atom->z) -= center - 
            0.5 * ( std::make_tuple(frame->box.box[0][0], frame->box.box[0][1], frame->box.box[0][2]) 
                    + std::make_tuple(frame->box.box[1][0], frame->box.box[1][1], frame->box.box[1][2])
                    + std::make_tuple(frame->box.box[2][0], frame->box.box[2][1], frame->box.box[2][2]));
    }
}

void PBCUtils::do_move_center_basedto_molecule(AmberMask &mask, std::shared_ptr<Frame> &frame) const {
    do_molecule_aggregate(frame);
    auto center_mol = find_molecule(mask, frame);
    auto center = center_mol->calc_geom_center(frame);
    for (auto &mol : frame->molecule_list) {
        auto r = mol->calc_geom_center(frame) - center;
        auto old = r;
        frame->image(r);
        r -= old;
        for (auto &atom : mol->atom_list) {
            std::tie(atom->x, atom->y, atom->z) += r;
        }
    }
    for (auto &atom : frame->atom_list){
        std::tie(atom->x, atom->y, atom->z) -= center - 
            0.5 * ( std::make_tuple(frame->box.box[0][0], frame->box.box[0][1], frame->box.box[0][2]) 
                    + std::make_tuple(frame->box.box[1][0], frame->box.box[1][1], frame->box.box[1][2])
                    + std::make_tuple(frame->box.box[2][0], frame->box.box[2][1], frame->box.box[2][2]));
    }
}

void PBCUtils::do_molecule_aggregate(std::shared_ptr<Frame> &frame) const {
    for (auto &mol : frame->molecule_list) {
        mol->do_aggregate(frame);
    }
}

void PBCUtils::do_move_all_atom_into_box(std::shared_ptr<Frame> &frame) const {
    for(auto &atom: frame->atom_list){
        frame->image(atom->x, atom->y, atom->z);
    }
}

void PBCUtils::doPBC(Trajconv::PBCType pbc_mode, AmberMask &mask, std::shared_ptr<Frame> &frame) const {
    if (pbc_mode == Trajconv::PBCType::AllIntoBox) {
        do_move_all_atom_into_box(frame);
    } else if (pbc_mode == Trajconv::PBCType::OneAtom) {
        do_move_center_basedto_atom(mask, frame);
    } else if (pbc_mode == Trajconv::PBCType::OneMol) {
        do_move_center_basedto_molecule(mask, frame);
    } else if (pbc_mode == Trajconv::PBCType::AtomGroup){
        do_move_center_basedto_atom_group(mask, frame);
    }
}
