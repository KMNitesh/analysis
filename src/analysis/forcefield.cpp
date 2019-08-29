//
// Created by xiamr on 3/17/19.
//
#include "config.h"
#include <memory>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include "forcefield.hpp"
#include "common.hpp"
#include "frame.hpp"
#include "ThrowAssert.hpp"


void Forcefield::read(const std::string &filename) {
    if (getFileType(filename) == FileType::PRM) {
        read_tinker_prm(filename);
    } else {
        std::cerr << "Unrecognized file extension" << std::endl;
        exit(6);
    }
    isvalid = true;
}

void Forcefield::read_tinker_prm(const std::string &filename) {
    if (!mapping.empty())
        return;

    std::fstream f;
    f.open(filename, std::ios::in);
    f.exceptions(std::ios::eofbit | std::ios::failbit | std::ios::badbit);
    std::string line;
    while (true) {
        try {
            std::getline(f, line);
        } catch (std::exception &e) {
            return;
        }
        auto field = split_quoted(line);
        if (field.empty()) continue;
        if (field[0] == "atom") {
            int type = stoi(field[1]);
            int class_num = stoi(field[2]);
            auto name = field[3];
            int at_no = stod(field[field.size() - 3]);
            double mass = stod(field[field.size() - 2]);
            auto res = field[4];
            mapping.emplace(type, AtomItem(type, class_num, name, res, at_no, mass));
        } else if (field[0] == "multipole") {
            int type = stoi(field[1]);
            double charge = stod(field[field.size() - 1]);
            auto it = mapping.find(type);
            throw_assert(it != mapping.end(), "Tinker prm is illegal!");
            it->second.charge = charge;
        } else if (field[0] == "charge") {
            int type = stoi(field[1]);
            double charge = stod(field[2]);
            auto it = mapping.find(type);
            throw_assert(it != mapping.end(), "Tinker prm is illegal!");
            it->second.charge = charge;
        }
    }
}


double Forcefield::find_mass(const std::shared_ptr<Atom> &atom) {
    if (atom->mass) return atom->mass.value();
    auto it = mapping.find(atom->typ);
    if (it != mapping.end()) {
        return it->second.mass;
    }
    throw std::runtime_error("atom mass not found");
}

void Forcefield::assign_forcefield(std::shared_ptr<Frame> &frame) {
    for (auto &atom : frame->atom_list) {
        auto it = mapping.find(atom->typ);
        throw_assert(it != mapping.end(), "Atom mass not found for type = " << atom->typ << "!");
        atom->mass = it->second.mass;
        atom->charge = it->second.charge.value();
        atom->setAtNo(it->second.at_no);
        atom->atom_name = it->second.name;
        atom->residue_name = it->second.res;
    }
}



