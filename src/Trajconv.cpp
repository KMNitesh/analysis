//
// Created by xiamr on 6/14/19.
//

#include "common.hpp"
#include "Trajconv.hpp"

#include "molecule.hpp"
#include "frame.hpp"
#include "atom.hpp"


void Trajconv::process(std::shared_ptr<Frame> &frame) {
    if (pbc_type == PBCType::OneAtom) {
        auto center_x = frame->atom_map[this->num]->x;
        auto center_y = frame->atom_map[this->num]->y;
        auto center_z = frame->atom_map[this->num]->z;

        for (auto &mol : frame->molecule_list) {
            for (auto &atom : mol->atom_list) {
                atom->x -= center_x;
                atom->y -= center_y;
                atom->z -= center_z;
            }
        }
    } else if (pbc_type == PBCType::OneMol) {
        auto center_mol = frame->atom_map[this->num]->molecule.lock();
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
    if (pbc_type != PBCType::None) {
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

    if (step == 0) {
        if (enable_gro) {
            GROWriter w;
            w.open(grofilename);
            w.write(frame);
            w.close();
        }
        if (enable_xtc) xtc.open(xtcfilename);
        if (enable_trr) trr.open(trrfilename);
        if (enable_mdcrd) mdcrd.open(mdcrdfilename);
    }
    if (enable_xtc) xtc.write(frame);
    if (enable_trr) trr.write(frame);
    if (enable_mdcrd) mdcrd.write(frame);
    step++;
}

void Trajconv::print() {
    if (enable_xtc) xtc.close();
    if (enable_trr) trr.close();
    if (enable_mdcrd) mdcrd.close();
}

void Trajconv::readInfo() {
    for (;;) {
        grofilename = choose_file("output gro file : ", false, "gro", true);
        if (grofilename.empty()) {
            enable_gro = false;
        }

        xtcfilename = choose_file("output xtc file : ", false, "xtc", true);
        if (xtcfilename.empty()) {
            enable_xtc = false;
        }

        trrfilename = choose_file("output trr file : ", false, "trr", true);
        if (trrfilename.empty()) {
            enable_trr = false;
        }

        mdcrdfilename = choose_file("output amber nc file : ", false, "nc", true);
        if (mdcrdfilename.empty()) {
            enable_mdcrd = false;
        }

        if (enable_gro || enable_xtc || enable_trr || enable_mdcrd) {
            break;
        }
        std::cerr << "ERROR !! none of output selected !\n";
    }
    std::cout << "PBC transform option\n";
    std::cout << "(0) Do nothing\n";
    std::cout << "(1) Make atom i as center\n";
    std::cout << "(2) Make molecule i as center\n";

    while (true) {
        int option = choose(0, 2, "Choose : ");
        switch (option) {
            case 0:
                pbc_type = PBCType::None;
                break;
            case 1:
                pbc_type = PBCType::OneAtom;
                num = choose(1, std::numeric_limits<int>::max(), "Plese enter the atom NO. : ");
                break;
            case 2:
                pbc_type = PBCType::OneMol;
                num = choose(1, std::numeric_limits<int>::max(),
                             "Plese enter one atom NO. that the molecule include: ");
                break;
            default:
                std::cerr << "option not found !\n";
                continue;
        }
        break;
    }
}

void Trajconv::fastConvertTo(std::string target) {
    boost::trim(target);
    if (target.empty()) {
        std::cerr << "ERROR !! empty target trajectory file \n";
        exit(EXIT_FAILURE);
    }
    pbc_type = PBCType::None;
    auto ext = ext_filename(target);
    if (ext == "xtc") {
        enable_xtc = true;
        enable_trr = false;
        enable_gro = false;
        enable_mdcrd = false;
        xtcfilename = target;
    } else if (ext == "trr") {
        enable_xtc = false;
        enable_trr = true;
        enable_gro = false;
        enable_mdcrd = false;
        trrfilename = target;
    } else if (ext == "nc") {
        enable_xtc = false;
        enable_trr = false;
        enable_gro = false;
        enable_mdcrd = true;
        mdcrdfilename = target;
    } else {
        std::cerr << "ERROR !!  wrong type of target trajectory file (" << target << ")\n";
        exit(EXIT_FAILURE);
    }
}


