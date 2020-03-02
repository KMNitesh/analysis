#include <boost/range/irange.hpp>
#include "PrmtopReader.hpp"
#include "topology_utils.hpp"
#include "data_structure/atom.hpp"
#include "data_structure/molecule.hpp"
#include "data_structure/frame.hpp"
#include "trajectory_reader/PrmtopParser.hpp"
#include "utils/common.hpp"

std::shared_ptr<Frame> PrmtopReader::read(const std::string &filename) {
    auto frame = std::make_shared<Frame>();

    std::ifstream ifs(filename);
    auto prmtop = PrmtopParser::parse(ifs);
    if (!prmtop) {
        std::cerr << "Error read prmtop file\n";
        exit(EXIT_FAILURE);
    }

    frame->title = prmtop->title;

    constexpr double charge_unit_convert_factor = 1 / 18.222615;

    for (auto i : boost::irange(prmtop->pointers[PrmtopStruct::NATOM])) {
        auto atom = std::make_shared<Atom>();
        atom->seq = i + 1;
        atom->atom_name = prmtop->atom_name[i];
        atom->type_name = prmtop->amber_atom_type[i];
        atom->charge = prmtop->charge[i] * charge_unit_convert_factor;
        atom->mass = prmtop->mass[i];

        atom->x = atom->y = atom->z = 0.0;
        frame->atom_list.push_back(atom);
        frame->atom_map[atom->seq] = atom;
    }

    int residue_num = 1;
    for (auto it = std::begin(prmtop->residue_pointer), next_it = it + 1;;
         it = next_it, ++next_it, ++residue_num) {

        if (next_it != std::end(prmtop->residue_pointer)) {
            for (auto i : boost::irange(*it, *next_it)) {
                auto &atom = frame->atom_map[i];
                atom->residue_num = residue_num;
                atom->residue_name = prmtop->residue_label[residue_num - 1];
            }
        } else {
            for (auto i : boost::irange<int>(*it, frame->atom_list.size() + 1)) {
                auto &atom = frame->atom_map[i];
                atom->residue_num = residue_num;
                atom->residue_name = prmtop->residue_label[residue_num - 1];
            }
            break;
        }
    }

    auto all_bonds = join(prmtop->bonds_inc_hydrogen, prmtop->bonds_without_hydrogen);
    for (auto it = std::begin(all_bonds); it != std::end(all_bonds); it += 3) {
        auto atom_num1 = (*it) / 3 + 1;
        auto atom_num2 = (*(it + 1)) / 3 + 1;

        auto &atom1 = frame->atom_map[atom_num1];
        auto &atom2 = frame->atom_map[atom_num2];

        atom1->con_list.push_back(atom_num2);
        atom2->con_list.push_back(atom_num1);
    }

    topology_utils::assgin_atom_to_molecule(frame);
    frame->build_graph();
    return frame;
}
