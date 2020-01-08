//
// Created by xiamr on 6/14/19.
//

#include "PrintTopolgy.hpp"

#include "dsl/center_selection_grammar.hpp"
#include "trajectory_reader/trajectoryreader.hpp"
#include "data_structure/frame.hpp"
#include "data_structure/atom.hpp"
#include "data_structure/forcefield.hpp"
#include "data_structure/molecule.hpp"

namespace qi = boost::spirit::qi;

void PrintTopolgy::action(const std::string &topology_filename) {
    TrajectoryReader reader;
    reader.add_topology(topology_filename);
    auto frame = reader.readTopology();
    if (forcefield.isValid()) {
        forcefield.assign_forcefield(frame);
    }

    enum class Mode {
        Mass,
        Geom,
        Noop,
    } mode;

    std::cout << std::left << "Total atom numbers : " << std::setw(10) << frame->atom_list.size()
              << "Total molecule numbers : " << frame->molecule_list.size() << '\n';
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
                      << " 3. ambermask\n"
                      << " 4. help\n"
                      << " 5. quit\n";
            continue;
        }

        std::cout << boost::format("%-6s %-7s %4s %-7s %4s %-6s  %8s %8s  %8s%8s%8s\n")
                     % "#Atom" % "Name" % "#Res" % "Name" % "#Mol" % "Type" % "Charge" % "Mass"
                     % "X(Ang)" % "Y(Ang)" % "Z(Ang)";
        double weight = 0;
        double sum_x = 0.0;
        double sum_y = 0.0;
        double sum_z = 0.0;
        const boost::format fmt{"%6d %-7s %4s %-7s %4s %-6s  %8s %8s  %8.3f%8.3f%8.3f\n"};
        for (auto &atom : frame->atom_list) {
            if (Atom::is_match(atom, ast)) {
                std::cout << boost::format(fmt)
                             % atom->seq % atom->atom_name
                             % (atom->residue_num ? std::to_string(atom->residue_num.get()) : "-")
                             % (atom->residue_name ? atom->residue_name.get() : "-")
                             % atom->molecule.lock()->sequence
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
