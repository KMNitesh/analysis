#include <utility>

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
#include <unordered_map>
#include <boost/optional.hpp>
#include "atom.hpp"

class Frame;

class AtomItem {
public:
    AtomItem(int typ, int classNum, std::string name, double mass) :
            typ(typ), class_num(classNum), name(std::move(name)), mass(mass) {}

    int typ;
    int class_num;

    std::string name;
    double mass;

    boost::optional<double> charge;
};

class Forcefield {
public:

    void read(const std::string &filename);

    void read_tinker_prm(const std::string &filename);

    double find_mass(const std::shared_ptr<Atom> &atom);

    void assign_forcefield(std::shared_ptr<Frame> &frame);

private:
    std::unordered_map<int, AtomItem> mapping;
};


#endif //TINKER_FORCEFIELD_HPP
