#include "topology_utils.hpp"

#include "data_structure/atom.hpp"
#include "data_structure/frame.hpp"
#include "data_structure/molecule.hpp"

void topology_utils::add_to_mol(std::shared_ptr<Atom> &atom, std::shared_ptr<Molecule> &mol,
                                std::shared_ptr<Frame> &frame) {
    if (atom->molecule.lock()) return;
    mol->atom_list.push_back(atom);
    atom->molecule = mol;
    for (auto &i : atom->con_list) add_to_mol(frame->atom_map[i], mol, frame);
}

void topology_utils::assgin_atom_to_molecule(std::shared_ptr<Frame> &frame) {
    for (auto &atom : frame->atom_list) {
        if (!atom->molecule.lock()) {
            auto molecule = std::make_shared<Molecule>();
            add_to_mol(atom, molecule, frame);
            frame->molecule_list.push_back(molecule);
        }
    }
}
