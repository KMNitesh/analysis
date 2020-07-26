#include "PrmtopReader.hpp"

#include <boost/range/irange.hpp>

#include "data_structure/atom.hpp"
#include "data_structure/frame.hpp"
#include "data_structure/molecule.hpp"
#include "topology_utils.hpp"
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
    for (auto it = std::begin(prmtop->residue_pointer), next_it = it + 1;; it = next_it, ++next_it, ++residue_num) {
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
        auto atom_num1 = (*(it + 0)) / 3 + 1;
        auto atom_num2 = (*(it + 1)) / 3 + 1;

        auto &atom1 = frame->atom_map[atom_num1];
        auto &atom2 = frame->atom_map[atom_num2];

        atom1->con_list.push_back(atom_num2);
        atom2->con_list.push_back(atom_num1);

        auto bond_index = *(it + 2);

        frame->f_bond_params.emplace(std::array{atom1, atom2},
                                     Frame::harmonic{.krA = prmtop->bond_force_constant[bond_index],
                                                     .rA = prmtop->bond_equil_value[bond_index]});
    }

    auto all_angles = join(prmtop->angles_inc_hydrogen, prmtop->angles_without_hydrogen);
    for (auto it = std::begin(all_angles); it != std::end(all_angles); it += 4) {
        auto atom_num1 = (*(it + 0)) / 3 + 1;
        auto atom_num2 = (*(it + 1)) / 3 + 1;
        auto atom_num3 = (*(it + 2)) / 3 + 1;

        auto angle_index = *(it + 3);

        auto &atom1 = frame->atom_map[atom_num1];
        auto &atom2 = frame->atom_map[atom_num2];
        auto &atom3 = frame->atom_map[atom_num3];

        frame->f_angle_params.emplace(std::array{atom1, atom2, atom3},
                                      Frame::harmonic{.krA = prmtop->angle_force_constant[angle_index],
                                                      .rA = prmtop->angle_equil_value[angle_index]});
    }

    auto all_dihedrals = join(prmtop->dihedrals_inc_hydrogen, prmtop->dihedrals_without_hydrogen);
    for (auto it = std::begin(all_dihedrals); it != std::end(all_dihedrals); it += 5) {
        auto atom_num1 = (*(it + 0)) / 3 + 1;
        auto atom_num2 = (*(it + 1)) / 3 + 1;
        auto atom_num3 = (*(it + 2)) / 3 + 1;
        auto atom_num4 = (*(it + 3)) / 3 + 1;

        auto dihedral_index = *(it + 4);

        auto &atom1 = frame->atom_map[atom_num1];
        auto &atom2 = frame->atom_map[atom_num2];
        auto &atom3 = frame->atom_map[std::abs(atom_num3)];
        auto &atom4 = frame->atom_map[std::abs(atom_num4)];

        std::array key{atom1, atom2, atom3, atom4};

        Frame::pdihs param{.phiA = prmtop->dihedral_phase[dihedral_index],
                           .cpA = prmtop->dihedral_force_constant[dihedral_index],
                           .mult = static_cast<int>(std::abs(prmtop->dihedral_periodicity[dihedral_index]))};

        if (atom_num4 > 0) {
            frame->f_dihedral_params.emplace(std::move(key), param);
        } else {
            frame->f_improper_dihedral_params.emplace(std::move(key), param);
        }
    }

    topology_utils::assgin_atom_to_molecule(frame);
    frame->build_graph();
    return frame;
}
