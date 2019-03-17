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


void Forcefield::read(const std::string &filename) {
    auto ext = ext_filename(filename);
    if (ext == "prm") {
        read_tinker_prm(filename);
    } else if (ext == "map") {
        read_mass_map(filename);
    } else {
        std::cerr << "Unrecognized file extension" << std::endl;
        exit(6);
    }
}

void Forcefield::read_tinker_prm(const std::string &filename) {
    std::fstream f;
    f.open(filename, std::ios::in);
    f.exceptions(std::ios::eofbit | std::ios::failbit | std::ios::badbit);
    std::string line;
    while (true) {
        try {
            std::getline(f, line);
            auto field = split(line);
            if (field.empty()) continue;
            if (field[0] != "atom") continue;
            int type = stoi(field[1]);
            double mass = stod(field[field.size() - 2]);
            items.emplace_back(Atom::AtomIndenter(std::make_shared<Atom::atom_types>(type)), mass);
        } catch (std::exception &e) {
            f.close();
            return;
        }
    }
}

void Forcefield::read_mass_map(const std::string &filename) {
    std::fstream f;
    f.open(filename, std::ios::in);
    f.exceptions(std::ios::eofbit | std::ios::failbit | std::ios::badbit);
    std::string line;
    while (true) {
        try {
            std::getline(f, line);
            boost::trim(line);
            auto field = split(line);
            if (field.empty()) continue;
            if (field.size() != 3) {
                std::cerr << "Force field mass map syntax error : " << line << std::endl;
                exit(8);
            }
            auto type = field[0];
            double mass = boost::lexical_cast<double>(field[2]);

            if (type == "name") {
                items.emplace_back(Atom::AtomIndenter(std::make_shared<Atom::atom_name_nums>(field[1])), mass);
            } else if (type == "type") {
                items.emplace_back(Atom::AtomIndenter(std::make_shared<Atom::atom_types>(field[1])), mass);
            } else if (type == "typenumber") {
                items.emplace_back(
                        Atom::AtomIndenter(std::make_shared<Atom::atom_types>(boost::lexical_cast<int>(field[1]))), mass);
            } else {
                std::cerr << "unrecognized keyword : " << type << std::endl;
                exit(8);
            }

        } catch (boost::bad_lexical_cast &e) {
            std::cerr << "boost::bad_lexical_cast  " << e.what() << std::endl;
            exit(9);
        } catch (std::exception &e) {
            f.close();
            return;
        }
    }
}

double Forcefield::find_mass(const std::shared_ptr<Atom> &atom) {
    for (auto &id : items) {
        if (Atom::is_match(atom, id.first)) return id.second;
    }
    std::cerr << "Atom mass not found !" << std::endl;
    exit(7);
}

