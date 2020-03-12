//
// Created by xiamr on 8/29/19.
//

#include "GmxTopologyPrinter.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/range/adaptors.hpp>

#include "data_structure/atom.hpp"
#include "data_structure/forcefield.hpp"
#include "data_structure/frame.hpp"
#include "trajectory_reader/trajectoryreader.hpp"
#include "utils/common.hpp"

void GmxTopologyPrinter::print(const std::string &topolgy, const std::string &prm, const std::string &out) {
    TrajectoryReader reader;
    reader.add_trajectoy_file(topolgy);
    auto frame = reader.readOneFrame();
    forcefield.read(prm);
    forcefield.assign_forcefield(frame);

    std::ofstream os(out);

    printFrame(frame, os);
}

void GmxTopologyPrinter::printFrame(const std::shared_ptr<Frame> &frame, std::ofstream &os) {
    os << "[ defaults ]\n"
          ";   nbfunc  comb-rule   gen-pairs   fudgeLJ fudgeQQ\n"
          "    1           2           no          1       1\n"
          ";nbfunc=1 for LJ; 2 for LB rule and 3 for geometric mixing rule.\n"
          "\n";

    os << "[ atomtypes ]\n"
          ";   name    at_no   mass        charge  ptype   sigma   epsilon \n";

    std::unordered_map<std::string, std::tuple<int, double, double>> atomtypes;

    for (auto &atom : frame->atom_list) {
        auto name = std::string("A") + std::to_string(atom->typ);
        atomtypes[name] = {atom->getAtNo().value(), atom->mass.value(), atom->charge.value()};
    }

    for (auto &types : atomtypes) {
        os << boost::format(" %-6s  %-6s %10.5f %10.5f ") % types.first % std::get<0>(types.second) %
                  std::get<1>(types.second) % std::get<2>(types.second)
           << "  A       0.295   0.530\n";
    }

    os << "\n[ moleculetype ]\n"
          "; molname      nrexcl\n"
          "TinkerMOL         3 \n";

    os << "\n[ atoms ]\n"
          ";   nr  type  resi  res  atom  cgnr     charge      mass\n";

    std::shared_ptr<Molecule> last_mol;

    int mol_index = 0;

    for (auto &atom : frame->atom_list) {
        if (atom->molecule.lock() != last_mol) {
            last_mol = atom->molecule.lock();
            ++mol_index;
        }
        std::string resname;
        if (boost::starts_with(atom->residue_name.value(), "AMOEBA Water")) {
            resname = "WAT";
        } else {
            resname = split(atom->residue_name.value()).back();
        }
        os << boost::format("%6d %4s %5d %5s %5s %4d %10.6f %11.6f\n") % atom->seq %
                  (std::string("A") + std::to_string(atom->typ)) % mol_index % resname % atom->atom_name % atom->seq %
                  atom->charge.value() % atom->mass.value();
    }

    os << "\n[ bonds ]\n"
          ";   ai     aj funct   r             k\n";

    for (auto &atom : frame->atom_list) {
        for (auto i : atom->con_list) {
            if (i > atom->seq) {
                os << "   " << atom->seq << "  " << i << "      1    1.0100e-01    3.6317e+05 ; " << atom->atom_name
                   << " - " << frame->atom_map[i]->atom_name << '\n';
            }
        }
    }

    os << "\n[ system ]\n"
          " TinkerMOL\n"
          "\n"
          "[ molecules ]\n"
          "; Compound        nmols\n"
          " TinkerMOL         1\n";
}
