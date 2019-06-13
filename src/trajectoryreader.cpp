//
// Created by xiamr on 3/17/19.
//
#include "config.h"
#include <list>
#include <iostream>
#include <fstream>
#include <boost/regex.hpp>

#include "common.hpp"
#include "trajectoryreader.hpp"
#include "atom.hpp"
#include "molecule.hpp"
#include "frame.hpp"

namespace gmx {

#include "gromacs/fileio/tpxio.h"
#include "gromacs/topology/mtop_util.c"

}


void TrajectoryReader::close() {
    if (isnetcdf) {
        netcdfClose(&NC);
        isnetcdf = false;
    } else if (istrr) {
        gmx::close_trn(fio);
        fio = nullptr;
        istrr = false;
    } else if (isxtc) {
        gmx::sfree(x);
        x = nullptr;
        gmx::close_xtc(fio);
        fio = nullptr;

        isxtc = false;
    } else {
        position_file.close();
    }
    if (velocity_file.is_open()) velocity_file.close();
}

void TrajectoryReader::readOneFrameVel() {
    std::string line;
    std::vector<std::string> field;
    getline(velocity_file, line);
    for (auto &atom : frame->atom_list) {
        std::getline(velocity_file, line);
        field = split(line);
        auto len1 = field[2].length();
        auto len2 = field[3].length();
        auto len3 = field[4].length();
        atom->vx = std::stod(field[2].replace(len1 - 4, 1, "E"));
        atom->vy = std::stod(field[3].replace(len2 - 4, 1, "E"));
        atom->vz = std::stod(field[4].replace(len3 - 4, 1, "E"));
    }
}

void add_to_mol(std::shared_ptr<Atom> &atom, std::shared_ptr<Molecule> &mol, std::shared_ptr<Frame> &frame) {
    if (atom->molecule.lock()) return;
    mol->atom_list.push_back(atom);
    atom->molecule = mol;
    for (auto &i : atom->con_list)
        add_to_mol(frame->atom_map[i], mol, frame);
}


std::shared_ptr<Frame> TrajectoryReader::readOneFrameTraj() {
    char str[256];
    if (!frame) {
        frame = std::make_shared<Frame>();
        int length;
        position_file.read((char *) &length, 4);
        int atom_num;
        position_file.read((char *) &atom_num, 4);
        position_file.read(str, length - 4);
        frame->title = std::string(str, static_cast<unsigned long>(length - 4));
        for (int i = 0; i < atom_num; i++) {
            auto atom_ptr = std::make_shared<Atom>();
            position_file.read((char *) &length, 4);
            position_file.read((char *) &atom_ptr->seq, 4);
            position_file.read((char *) &atom_ptr->typ, 4);
            position_file.read((char *) &str, length - 8);
            atom_ptr->atom_name = std::string(str, static_cast<unsigned long>(length - 8));
            int num;
            position_file.read((char *) &num, 4);
            int n;
            for (int k = 0; k < num; k++) {
                position_file.read((char *) &n, 4);
                atom_ptr->con_list.push_back(n);
            }

            frame->atom_list.push_back(atom_ptr);
            frame->atom_map[atom_ptr->seq] = atom_ptr;
        }
        for (auto &atom : frame->atom_list) {
            if (!atom->molecule.lock()) {
                auto molecule = std::make_shared<Molecule>();
                add_to_mol(atom, molecule, frame);
                frame->molecule_list.push_back(molecule);
            }
        }
    }
    float box[6];
    position_file.read((char *) box, 24);
    frame->a_axis = box[0];
    frame->b_axis = box[1];
    frame->c_axis = box[2];
    frame->alpha = box[3];
    frame->beta = box[4];
    frame->gamma = box[5];
    for (auto &atom : frame->atom_list) {
        float vect[3];
        position_file.read((char *) vect, 12);
        atom->x = vect[0];
        atom->y = vect[1];
        atom->z = vect[2];
    }
    frame->a_axis_half = frame->a_axis / 2;
    frame->b_axis_half = frame->b_axis / 2;
    frame->c_axis_half = frame->c_axis / 2;

    return frame;


}

int TrajectoryReader::readOneFrameNetCDF() {
    auto coord = new double[NC.ncatom3];
    double box[6];
    int ret = netcdfGetNextFrame(&NC, coord, box);
    if (ret == 0) return ret;
    int i = 0;
    for (auto &atom : frame->atom_list) {
        atom->x = coord[3 * i];
        atom->y = coord[3 * i + 1];
        atom->z = coord[3 * i + 2];
        i++;
    }
    frame->a_axis = box[0];
    frame->b_axis = box[1];
    frame->c_axis = box[2];
    frame->alpha = box[3];
    frame->beta = box[4];
    frame->gamma = box[5];

    if (frame->enable_bound) {
        frame->a_axis_half = frame->a_axis / 2;
        frame->b_axis_half = frame->b_axis / 2;
        frame->c_axis_half = frame->c_axis / 2;
    }

    delete[] coord;
    return ret;
}

int TrajectoryReader::readOneFrameTrr() {
    gmx::t_trnheader trnheader;
    gmx::gmx_bool bOK;
    gmx::rvec box[3];
    gmx::rvec *coord = nullptr;
    if (gmx::fread_trnheader(fio, &trnheader, &bOK)) {
        if (bOK) {
            coord = new gmx::rvec[trnheader.natoms];
            if (trnheader.box_size) {
                gmx::fread_htrn(fio, &trnheader, box, coord, NULL, NULL);
                translate(box, &(frame->a_axis), &(frame->b_axis), &(frame->c_axis),
                          &(frame->alpha), &(frame->beta), &(frame->gamma));
                if (frame->enable_bound) {
                    frame->a_axis_half = frame->a_axis / 2;
                    frame->b_axis_half = frame->b_axis / 2;
                    frame->c_axis_half = frame->c_axis / 2;
                }
            } else {
                if (frame->enable_bound) {
                    std::cerr << "WARING ! then trajectory have PBC enabled" << std::endl;
                    exit(1);
                }
                frame->a_axis = 0.0;
                frame->b_axis = 0.0;
                frame->c_axis = 0.0;
                frame->alpha = 0.0;
                frame->beta = 0.0;
                frame->gamma = 0.0;
                frame->enable_bound = false;
                gmx::fread_htrn(fio, &trnheader, NULL, coord, NULL, NULL);
            }
            if (static_cast<int>(frame->atom_list.size()) != trnheader.natoms) {
                std::cerr << "ERROR! the atom number do not match" << std::endl;
                exit(1);
            }
            int i = 0;
            for (auto &atom : frame->atom_list) {
                atom->x = coord[i][0] * 10;
                atom->y = coord[i][1] * 10;
                atom->z = coord[i][2] * 10;
                i++;
            }
            delete[] coord;
            return 1;
        }
    }
    return 0;

}

int TrajectoryReader::readOneFrameXtc() {
    gmx::matrix box;
    gmx::gmx_bool bOK;

    int ret;
    if (!x) {
        ret = gmx::read_first_xtc(fio, &natoms, &step, &time, box, &x, &prec, &bOK);
    } else {
        ret = gmx::read_next_xtc(fio, natoms, &step, &time, box, x, &prec, &bOK);
    }
    if (!ret) return 0;

    if (bOK) {
        if (natoms != static_cast<int>(frame->atom_list.size())) {
            std::cerr << "ERROR! the atom number do not match" << std::endl;
            exit(1);
        }
        if (frame->enable_bound) {
            translate(box, &(frame->a_axis), &(frame->b_axis), &(frame->c_axis),
                      &(frame->alpha), &(frame->beta), &(frame->gamma));
            frame->a_axis_half = frame->a_axis / 2;
            frame->b_axis_half = frame->b_axis / 2;
            frame->c_axis_half = frame->c_axis / 2;
        }
        int i = 0;
        for (auto &atom : frame->atom_list) {
            atom->x = x[i][0] * 10;
            atom->y = x[i][1] * 10;
            atom->z = x[i][2] * 10;
            i++;
        }
        return 1;
    }
    std::cerr << "\nWARNING: Incomplete frame at time "
              << std::scientific << time << std::defaultfloat << "\n";
    return 0;
}


std::shared_ptr<Frame> TrajectoryReader::readOneFrameTpr() {


    if (!frame) {
        frame = std::make_shared<Frame>();
        if (enable_binaray_file) frame->enable_bound = true;
    }


    gmx::t_state state;
    gmx::t_tpxheader tpx;
    gmx::t_inputrec ir;
    gmx::gmx_mtop_t mtop;
    gmx::t_topology top;

    gmx::read_tpxheader(topology_filename.c_str(), &tpx, TRUE, nullptr, nullptr);

    gmx::read_tpx_state(topology_filename.c_str(),
                        tpx.bIr ? &ir : nullptr,
                        &state, nullptr,
                        tpx.bTop ? &mtop : nullptr);

    top = gmx::gmx_mtop_t_to_t_topology(&mtop);


    auto &atoms = top.atoms;

    if (!tpx.bX) {
        std::cerr << "Tpr Topolgy do not have coordinates\n";
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < atoms.nr; i++) {
        auto atom = std::make_shared<Atom>();
        atom->seq = i + 1;
        atom->atom_name = (*(atoms.atomname[i]));
        atom->type_name = (*(atoms.atomtype[i]));
        atom->charge = atoms.atom[i].q;
        atom->mass = atoms.atom[i].m;

        atom->x = 10 * state.x[i][XX];
        atom->y = 10 * state.x[i][YY];
        atom->z = 10 * state.x[i][ZZ];

        atom->residue_name = *atoms.resinfo[atoms.atom[i].resind].name;
        atom->residue_num = atoms.resinfo[atoms.atom[i].resind].nr;
        frame->atom_list.push_back(atom);
        frame->atom_map[atom->seq] = atom;
    }


    auto ilist = &top.idef.il[gmx::F_BONDS];
    gmx::t_iatom *iatoms = ilist->iatoms;
    for (int i = 0; i < ilist->nr; i += 3) {
        auto type = *(iatoms++);
        auto ftype = top.idef.functype[type];

        int atom_num1 = *(iatoms++) + 1;
        int atom_num2 = *(iatoms++) + 1;


        auto atom1 = frame->atom_map[atom_num1];
        auto atom2 = frame->atom_map[atom_num2];

        atom1->con_list.push_back(atom_num2);
        atom2->con_list.push_back(atom_num1);

    }
    ilist = &top.idef.il[gmx::F_SETTLE];
    iatoms = ilist->iatoms;
    for (int i = 0; i < ilist->nr;) {
        auto type = *(iatoms++);
        auto ftype = top.idef.functype[type];

        auto nratoms = gmx::interaction_function[ftype].nratoms;

        int atom_num1 = *(iatoms++) + 1;
        auto atom1 = frame->atom_map[atom_num1];


        for (int j = 1; j < nratoms; j++) {
            int atom_num2 = *(iatoms++) + 1;
            auto atom2 = frame->atom_map[atom_num2];

            atom1->con_list.push_back(atom_num2);
            atom2->con_list.push_back(atom_num1);
        }

        i += 1 + nratoms;
    }


    if (first_time) {
        for (auto &atom : frame->atom_list) {
            if (!atom->molecule.lock()) {
                auto molecule = std::make_shared<Molecule>();
                add_to_mol(atom, molecule, frame);
                frame->molecule_list.push_back(molecule);
            }
        }
    }

    if (atoms.nres != boost::numeric_cast<int>(frame->molecule_list.size())) {
        std::cout << "Residue numbers and Molecule numbers not match\n";
//        exit(EXIT_FAILURE);
    }

    return frame;
}

std::shared_ptr<Frame> TrajectoryReader::readOneFrameMol2() {
    std::string line;
    std::vector<std::string> fields;

    enum class State {
        NONE, MOLECULE, ATOM, BOND, SUBSTRUCTURE
    } state;
    state = State::NONE;
    auto whitespace_regex = boost::regex("\\s+");
    // watch out memery leak

    if (!frame) {
        frame = std::make_shared<Frame>();
        if (enable_binaray_file) frame->enable_bound = true;
    }

    auto iter = frame->atom_list.begin();

    for (;;) {
        std::getline(position_file, line);
        boost::trim(line);
        if (line.empty()) continue;
        fields.clear();

        bool break_loop = false;

        if (boost::starts_with(line, "@<TRIPOS>MOLECULE")) {
            state = State::MOLECULE;
        } else if (boost::starts_with(line, "@<TRIPOS>ATOM")) {
            state = State::ATOM;
        } else if (boost::starts_with(line, "@<TRIPOS>BOND")) {
            state = State::BOND;
        } else if (boost::starts_with(line, "@<TRIPOS>SUBSTRUCTURE")) {
            state = State::SUBSTRUCTURE;
        } else {
            switch (state) {
                case State::MOLECULE:
                    continue;
                case State::ATOM: {
                    boost::regex_split(std::back_inserter(fields), line, whitespace_regex);
                    std::shared_ptr<Atom> atom;
                    if (first_time) {
                        atom = std::make_shared<Atom>();
                        atom->seq = boost::lexical_cast<int>(fields[0]);
                        atom->atom_name = fields[1];
                        atom->type_name = fields[5];
                        atom->residue_name = fields[7];
                        atom->residue_num = boost::lexical_cast<uint>(fields[6]);
                        atom->charge = boost::lexical_cast<double>(fields[8]);
                        frame->atom_list.push_back(atom);
                        frame->atom_map[atom->seq] = atom;
                    } else {
                        atom = *iter;
                        ++iter;
                    }
                    atom->x = boost::lexical_cast<double>(fields[2]);
                    atom->y = boost::lexical_cast<double>(fields[3]);
                    atom->z = boost::lexical_cast<double>(fields[4]);
                }
                    break;
                case State::BOND: {
                    if (!first_time) {
                        break_loop = true;
                        break;
                    }
                    boost::regex_split(std::back_inserter(fields), line, whitespace_regex);
                    int atom_num1 = boost::lexical_cast<int>(fields[1]);
                    int atom_num2 = boost::lexical_cast<int>(fields[2]);

                    auto atom1 = frame->atom_map[atom_num1];
                    auto atom2 = frame->atom_map[atom_num2];

                    atom1->con_list.push_back(atom_num2);
                    atom2->con_list.push_back(atom_num1);
                }
                    break;
                case State::SUBSTRUCTURE:
                    break_loop = true;
                    break;
                default:
                    break;
            }
        }
        if (break_loop) {
            break;
        }
    }
    if (first_time) {
        for (auto &atom : frame->atom_list) {
            if (!atom->molecule.lock()) {
                auto molecule = std::make_shared<Molecule>();
                add_to_mol(atom, molecule, frame);
                frame->molecule_list.push_back(molecule);
            }
        }
    }
    return frame;

}

std::shared_ptr<Frame> TrajectoryReader::readOneFrameArc() {
    int atom_num = 0;
    std::string line;
    std::vector<std::string> field;
    if (!frame) {
        frame = std::make_shared<Frame>();
        std::getline(position_file, line);
        field = split(line);
        atom_num = std::stoi(field[0]);
        frame->title = line.substr(line.rfind(field[0]) + field[0].size());
        boost::trim(frame->title);
    } else {
        std::getline(position_file, line);
        if (line.empty()) throw std::exception();
    }

    if (frame->enable_bound) {
        std::getline(position_file, line);
        field = split(line);
        if (field.empty()){
            throw std::exception();
        }
        frame->a_axis = std::stod(field[0]);
        frame->b_axis = std::stod(field[1]);
        frame->c_axis = std::stod(field[2]);
        frame->alpha = std::stod(field[3]);
        frame->beta = std::stod(field[4]);
        frame->gamma = std::stod(field[5]);
    }

    if (first_time) {
        bool first_loop = true;
        for (int i = 0; i < atom_num; i++) {
            std::getline(position_file, line);
            field = split(line);
            if (first_loop) {
                first_loop = false;
                if (field.empty()) throw std::exception();
                try {
                    std::stod(field[1]);
                    frame->a_axis = std::stod(field[0]);
                    frame->b_axis = std::stod(field[1]);
                    frame->c_axis = std::stod(field[2]);
                    frame->alpha = std::stod(field[3]);
                    frame->beta = std::stod(field[4]);
                    frame->gamma = std::stod(field[5]);
                    frame->enable_bound = true;
                    i--;
                    continue;
                } catch (std::exception &e) {

                }
            }

            auto atom = std::make_shared<Atom>();
            atom->seq = std::stoi(field[0]);
            atom->atom_name = field[1];
            atom->x = std::stod(field[2]);
            atom->y = std::stod(field[3]);
            atom->z = std::stod(field[4]);
            atom->typ = std::stoi(field[5]);
            for (size_t j = 6; j < field.size(); j++) {
                atom->con_list.push_back(std::stoi(field[j]));
            }
            frame->atom_list.push_back(atom);
            frame->atom_map[atom->seq] = atom;
        }

        for (auto &atom : frame->atom_list) {
            if (!atom->molecule.lock()) {
                auto molecule = std::make_shared<Molecule>();
                add_to_mol(atom, molecule, frame);
                frame->molecule_list.push_back(molecule);
            }
        }

    } else {
        for (auto &atom : frame->atom_list) {
            std::getline(position_file, line);
            field = split(line);
            atom->x = std::stod(field[2]);
            atom->y = std::stod(field[3]);
            atom->z = std::stod(field[4]);
        }
    }
    if (frame->enable_bound) {
        frame->a_axis_half = frame->a_axis / 2;
        frame->b_axis_half = frame->b_axis / 2;
        frame->c_axis_half = frame->c_axis / 2;
    }
    first_time = false;
    return frame;
}

void TrajectoryReader::open(const std::string &filename) {
    auto field = split(filename, ".");
    auto ext = ext_filename(filename);

    if (enable_binaray_file) {
        if (ext == "nc") {
            if (netcdfLoad(&NC, filename.c_str())) {
                std::cerr << "error open NETCDF file: " << filename << std::endl;
                exit(2);
            }
            isnetcdf = true;
            istrr = false;
            isxtc = false;
        } else if (ext == "trr") {
            fio = gmx::open_trn(filename.c_str(), "r");
            isnetcdf = false;
            istrr = true;
            isxtc = false;

        } else if (ext == "xtc") {
            fio = gmx::open_xtc(filename.c_str(), "r");
            isnetcdf = false;
            istrr = false;
            isxtc = true;
        }
    } else if (ext == "traj") {
        isbin = true;
        position_file.exceptions(std::ios::eofbit | std::ios::failbit | std::ios::badbit);
        position_file.open(filename, std::ios::in | std::ios::binary);
    } else {
        isbin = false;
        position_file.open(filename);
    }
    if (enable_read_velocity) {
        velocity_file.open(field[0] + ".vel");
        velocity_file.exceptions(std::ios::eofbit | std::ios::failbit | std::ios::badbit);
        if (velocity_file.good()) this->openvel = true;
        else {
            std::cerr << "Error open velocity file" << std::endl;
            exit(3);
        }
    }
}

void TrajectoryReader::add_filename(const std::string &filename) {
    arc_filename_list.push_back(filename);
}

void TrajectoryReader::add_topology(const std::string &filename) {
    enable_binaray_file = true;
    topology_filename = filename;
    std::string ext_name = ext_filename(filename);
    boost::to_lower(ext_name);

    if (ext_name == "arc" or ext_name == "xyz") {
        topology_type = TOPOLOGY_TYPE::ARC;
    } else if (ext_name == "mol2") {
        topology_type = TOPOLOGY_TYPE::MOL2;
    } else if (ext_name == "tpr") {
        topology_type = TOPOLOGY_TYPE::TPR;
    } else {
        std::cerr << " Error file type of topology file [ " << filename << "] " << std::endl;
        exit(1);
    }

}

std::shared_ptr<Frame> TrajectoryReader::readOneFrame() {
    if (first_time) {
        readTopology();
        std::string filename = arc_filename_list.front();
        open(filename);
        arc_filename_list.pop_front();
    }
    loop:
    try {
        if (isnetcdf) {
            if (!readOneFrameNetCDF()) {
                throw std::exception();
            }
        } else if (istrr) {
            if (!readOneFrameTrr()) {
                throw std::exception();
            }
        } else if (isxtc) {
            if (!readOneFrameXtc()) {
                throw std::exception();
            }
        } else if (isbin) {
            frame = readOneFrameTraj();
            frame->enable_bound = frame->a_axis != 0.0;
        } else {
            frame = readOneFrameArc();
        }

        if (openvel && frame) {
            readOneFrameVel();
        }
    } catch (std::exception &e) {
        if (arc_filename_list.empty()) {
            frame.reset();
            close();
        } else {
            close();
            std::string filename = arc_filename_list.front();
            open(filename);
            arc_filename_list.pop_front();
            goto loop;
        }
    }

    first_time = false;

    return frame;
}

std::shared_ptr<Frame> TrajectoryReader::readTopology() {
    if (enable_binaray_file) {
        position_file.open(topology_filename, std::ios::in);
        switch (topology_type) {
            case TOPOLOGY_TYPE::ARC:
                frame = readOneFrameArc();
                break;
            case TOPOLOGY_TYPE::MOL2:
                frame = readOneFrameMol2();
                break;
            case TOPOLOGY_TYPE::TPR:
                frame = readOneFrameTpr();
                break;
        }
        position_file.close();
    }
    return frame;
}
