//
// Created by xiamr on 6/14/19.
//

#include "PrintTopolgy.hpp"

#include "center_selection_grammar.hpp"
#include "trajectoryreader.hpp"
#include "frame.hpp"
#include "atom.hpp"
#include "forcefield.hpp"

namespace qi = boost::spirit::qi;

void PrintTopolgy::action(const std::string &topology_filename) {
    TrajectoryReader reader;
    reader.add_topology(topology_filename);
    auto frame = reader.readTopology();
    if (forcefield.isVaild()) {
        forcefield.assign_forcefield(frame);
    }

    int sequence = 1;
    for (auto &mol : frame->molecule_list) {
        mol->sequence = sequence++;
    }

    enum class Mode {
        Mass,
        Geom,
        Noop,
    } mode;

    for (;;) {
        beg:
        CenterRuleNode r;
        selectCentergroup(r, "> ");
        Atom::Node ast;

        if (auto ret = boost::get<std::shared_ptr<MassCenterRuleNode>>(&r)) {
            ast = (*ret)->SelectionMask;
            mode = Mode::Mass;
        } else if (auto ret = boost::get<std::shared_ptr<GeomCenterRuleNode>>(&r)) {
            ast = (*ret)->SelectionMask;
            mode = Mode::Geom;
        } else if (auto ret = boost::get<std::shared_ptr<NoopRuleNode>>(&r)) {
            ast = (*ret)->SelectionMask;
            mode = Mode::Noop;
        } else if (boost::get<std::shared_ptr<QuitRuleNode>>(&r)) {
            break;
        } else if (boost::get<std::shared_ptr<HelpRuleNode>>(&r)) {
            std::cout << "Help : \n"
                      << " 1. com  of ambermask\n"
                      << " 2. geom of ambermask\n"
                      << " 3. help\n"
                      << " 4. quit\n";
            continue;
        }

        std::cout << boost::format("%-6s %-7s %4s %-7s %4s %-6s  %6s %8s  %8s%8s%8s\n")
                     % "#Atom" % "Name" % "#Res" % "Name" % "#Mol" % "Type" % "Charge" % "Mass"
                     % "X(Ang)" % "Y(Ang)" % "Z(Ang)";
        double weight = 0;
        double sum_x = 0.0;
        double sum_y = 0.0;
        double sum_z = 0.0;
        for (auto &atom : frame->atom_list) {
            if (is_match_impl(atom, ast)) {
                std::cout << boost::format("%6d %-7s %4s %-7s %4s %-6s % 8f %8s %8.3f%8.3f%8.3f\n")
                             % atom->seq % atom->atom_name
                             % (atom->residue_num ? boost::lexical_cast<std::string>(atom->residue_num.get()) : "-")
                             % (atom->residue_name ? atom->residue_name.get() : "-")
                             % atom->molecule.lock()->sequence.get()
                             % atom->type_name
                             % (atom->charge ? (boost::format("%8.4f") % atom->charge.get()).str() : "-")
                             % (atom->mass ? (boost::format("%8.4f") % atom->mass.get()).str() : "-")
                             % atom->x
                             % atom->y
                             % atom->z;
                switch (mode) {
                    case Mode::Mass:
                        if (!atom->mass) {
                            std::cerr << "atom mass not available !\n";
                            goto beg;
                        }
                        sum_x += atom->x * atom->mass.get();
                        sum_y += atom->y * atom->mass.get();
                        sum_z += atom->z * atom->mass.get();
                        weight += atom->mass.get();
                        break;
                    case Mode::Geom:
                        sum_x += atom->x;
                        sum_y += atom->y;
                        sum_z += atom->z;
                        weight++;
                        break;
                    case Mode::Noop:
                        break;
                }
            }
        }
        if (weight != 0.0) {
            switch (mode) {
                case Mode::Mass:

                    std::cout << boost::format("Mass Center %8s%8s%8s\n") % "X(Ang)" % "Y(Ang)" % "Z(Ang)";
                    std::cout << boost::format("            %8.3f%8.3f%8.3f\n")
                                 % (sum_x / weight)
                                 % (sum_y / weight)
                                 % (sum_z / weight);

                    break;
                case Mode::Geom:
                    std::cout << boost::format("Geom Center %8s%8s%8s\n") % "X(Ang)" % "Y(Ang)" % "Z(Ang)";
                    std::cout << boost::format("            %8.3f%8.3f%8.3f\n")
                                 % (sum_x / weight)
                                 % (sum_y / weight)
                                 % (sum_z / weight);
                    break;
                case Mode::Noop:
                    break;
            }
        }

    }
}
