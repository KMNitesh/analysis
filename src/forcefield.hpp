//
// Created by xiamr on 3/17/19.
//

#ifndef TINKER_FORCEFIELD_HPP
#define TINKER_FORCEFIELD_HPP

#include "config.h"
#include <string>
#include <memory>
#include <list>
#include <utility>
#include "atom.hpp"

class AtomItem {
public:
    int typ;
    double mass;
};

class Forcefield {
public:

    void read(const std::string &filename);

    void read_tinker_prm(const std::string &filename);

    void read_mass_map(const std::string &filename);

    double find_mass(const std::shared_ptr<Atom> &atom);

private:
    std::list<std::pair<Atom::AtomIndenter, double>> items;
};



#endif //TINKER_FORCEFIELD_HPP
