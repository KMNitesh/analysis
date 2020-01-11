#ifndef TINKER_TOPOLOGY_UTILS_HPP
#define TINKER_TOPOLOGY_UTILS_HPP

#include <memory>

class Atom;

class Molecule;

class Frame;

class topology_utils {
public:
    static void assgin_atom_to_molecule(std::shared_ptr<Frame> &frame);

    static void add_to_mol(std::shared_ptr<Atom> &atom, std::shared_ptr<Molecule> &mol, std::shared_ptr<Frame> &frame);
};


#endif //TINKER_TOPOLOGY_UTILS_HPP
