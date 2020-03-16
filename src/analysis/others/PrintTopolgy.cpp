//
// Created by xiamr on 6/14/19.
//

#include "PrintTopolgy.hpp"

#include "data_structure/atom.hpp"
#include "data_structure/forcefield.hpp"
#include "data_structure/frame.hpp"
#include "data_structure/molecule.hpp"
#include "dsl/center_selection_grammar.hpp"
#include "trajectory_reader/trajectoryreader.hpp"

namespace qi = boost::spirit::qi;

void PrintTopolgy::action(const std::string &topology_filename) {
    TrajectoryReader reader;

    auto t1 = std::chrono::high_resolution_clock::now();

    reader.set_topology(topology_filename);
    auto frame = reader.readTopology();
    if (forcefield.isValid()) {
        forcefield.assign_forcefield(frame);
    }

    enum class Mode {
        Mass,
        Geom,
        Noop,
    } mode;

    std::cout
        << std::left << "Total atom numbers : " << std::setw(10) << frame->atom_list.size()
        << "Total molecule numbers : " << std::setw(10) << frame->molecule_list.size() << "Load time : "
        << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t1).count()
        << " ms\n";

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

        std::cout << boost::format("%-6s %-7s %4s %-7s %4s %-6s  %8s %8s  %8s%8s%8s") % "#Atom" % "Name" % "#Res" %
                         "Name" % "#Mol" % "Type" % "Charge" % "Mass" % "X(Ang)" % "Y(Ang)" % "Z(Ang)";
        if (frame->atom_list.front()->lj_param.has_value())
            std::cout << boost::format("%14s%14s%14s%14s") % "LJ_SR c6" % "LJ_SR c12" % "sigma(nm)" % "eps(kJ/mol)";
        std::cout << '\n';

        double weight = 0;
        double sum_x = 0.0;
        double sum_y = 0.0;
        double sum_z = 0.0;
        const boost::format fmt{"%6d %-7s %4s %-7s %4s %-6s  %8s %8s  %8.3f%8.3f%8.3f"};

        for (auto &atom : frame->atom_list) {
            if (Atom::is_match(atom, ast)) {
                std::cout << boost::format(fmt) % atom->seq % atom->atom_name %
                                 (atom->residue_num ? std::to_string(atom->residue_num.get()) : "-") %
                                 (atom->residue_name ? atom->residue_name.get() : "-") %
                                 atom->molecule.lock()->sequence % atom->type_name %
                                 (atom->charge ? (boost::format("%8.4f") % atom->charge.get()).str() : "-") %
                                 (atom->mass ? (boost::format("%8.4f") % atom->mass.get()).str() : "-") % atom->x %
                                 atom->y % atom->z;
                if (atom->lj_param.has_value()) {
                    const auto c6 = (*atom->lj_param).c6;
                    const auto c12 = (*atom->lj_param).c12;
                    const auto sigma = c6 != 0.0 ? std::pow(c12 / c6, 1.0 / 6) : 0.0;
                    const auto epsilon = c12 != 0.0 ? c6 * c6 / (4 * c12) : 0.0;
                    std::cout << boost::format("%14e%14e%14e%14e") % c6 % c12 % sigma % epsilon;
                }
                std::cout << '\n';
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
                std::cout << boost::format("            %8.3f%8.3f%8.3f\n") % (sum_x / weight) % (sum_y / weight) %
                                 (sum_z / weight);

                break;
            case Mode::Geom:
                std::cout << boost::format("Geom Center %8s%8s%8s\n") % "X(Ang)" % "Y(Ang)" % "Z(Ang)";
                std::cout << boost::format("            %8.3f%8.3f%8.3f\n") % (sum_x / weight) % (sum_y / weight) %
                                 (sum_z / weight);
                break;
            case Mode::Noop:
                break;
            }
        }
    }
}
