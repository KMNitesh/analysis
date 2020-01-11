

#include <boost/algorithm/string_regex.hpp>
#include <boost/lexical_cast.hpp>
#include "data_structure/atom.hpp"
#include "data_structure/molecule.hpp"
#include "data_structure/frame.hpp"
#include "utils/common.hpp"
#include "Mol2Reader.hpp"
#include "topology_utils.hpp"

std::shared_ptr<Frame> Mol2Reader::read(const std::string &filename) {
    std::string line;
    std::vector<std::string> fields;

    std::ifstream position_file(filename);

    enum class State {
        NONE, MOLECULE, ATOM, BOND, SUBSTRUCTURE
    } state;
    state = State::NONE;
    auto whitespace_regex = boost::regex("\\s+");

    auto frame = std::make_shared<Frame>();

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

                    atom = std::make_shared<Atom>();
                    atom->seq = boost::lexical_cast<int>(fields[0]);
                    atom->atom_name = fields[1];
                    atom->type_name = fields[5];
                    atom->residue_name = fields[7];
                    atom->residue_num = boost::lexical_cast<uint>(fields[6]);
                    atom->charge = boost::lexical_cast<double>(fields[8]);
                    frame->atom_list.push_back(atom);
                    frame->atom_map[atom->seq] = atom;

                    atom->x = boost::lexical_cast<double>(fields[2]);
                    atom->y = boost::lexical_cast<double>(fields[3]);
                    atom->z = boost::lexical_cast<double>(fields[4]);
                }
                    break;
                case State::BOND: {
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
    topology_utils::assgin_atom_to_molecule(frame);
    return frame;
}
