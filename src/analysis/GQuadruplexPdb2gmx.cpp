//
// Created by xiamr on 10/18/19.
//

#include "GQuadruplexPdb2gmx.hpp"
#include "common.hpp"
#include <boost/algorithm/string.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/numeric.hpp>
#include "RMSDCal.hpp"

namespace {
    class PDBAtom {
    public:
        PDBAtom(int atomNum, std::string atomName, double x, double y, double z)
                : atom_num(atomNum), atom_name(std::move(atomName)), x(x), y(y), z(z) {}

        int atom_num;
        std::string atom_name;
        double x, y, z;
    };

    class Residue {
    public:
        Residue(std::string resName, int residueNum) : res_name(std::move(resName)), residue_num(residueNum) {}

        std::string res_name;
        int residue_num;
        std::vector<PDBAtom> atoms;
    };

    std::ostream &operator<<(std::ostream &os, const Residue &residue) {
        int min = boost::min_element(residue.atoms, [](auto &lhs, auto &rhs) {
            return lhs.atom_num < rhs.atom_num;
        })->atom_num;
        for (auto i : boost::irange(std::size_t(0), residue.atoms.size())) {
            os << boost::format("ATOM  %5d %-4s %+3s  %4d    %8.3f%8.3f%8.3f  1.00  0.00\n")
                  % (i + min)
                  % residue.atoms[i].atom_name
                  % residue.res_name
                  % residue.residue_num
                  % residue.atoms[i].x
                  % residue.atoms[i].y
                  % residue.atoms[i].z;
        }
        return os;
    }

    std::map<std::string, std::string> replace{
            {"H5'",  "H5'1"},
            {"H2'",  "H2'1"},
            {"H5''", "H5'2"},
            {"H2''", "H2'2"},
            {"OP1",  "O1P"},
            {"OP2",  "O2P"},
            {"HO5'", "H5T"},
            {"HO3",  "H3T"}
    };

    enum class MODE {
        BONDEDTYPES,
        RESIDUE,
        ATOMS,
        BONDS,
        IMPROPERS,
    };

    std::map<std::string, MODE> mode_mapping{
            {"bondedtypes", MODE::BONDEDTYPES},
            {"atoms",       MODE::ATOMS},
            {"bonds",       MODE::BONDS},
            {"impropers",   MODE::IMPROPERS}
    };

    std::map<std::string, std::map<std::string, int>> parse_rtp(std::ifstream &&ifs) {
        std::map<std::string, std::map<std::string, int>> ret;
        std::map<std::string, std::map<std::string, int>>::iterator current_it;
        auto mode = MODE::BONDEDTYPES;
        std::string line;
        while (!ifs.eof()) {
            std::getline(ifs, line);
            boost::trim(line);
            if (line.empty()) continue;
            if (boost::starts_with(line, ";")) continue;
            if (boost::starts_with(line, "[") and boost::ends_with(line, "]")) {
                line.erase(line.length() - 1);
                line.erase(0);
                auto &key = line;
                auto it = mode_mapping.find(key);
                mode = MODE::RESIDUE;
                if (it != mode_mapping.end()) {
                    mode = it->second;
                } else {
                    auto[i, ok] = ret.emplace(key, std::map<std::string, int>{});
                    assert(ok);
                    current_it = i;
                }
            } else {
                if (mode == MODE::ATOMS) {
                    auto field = split(line);
                    assert(field.size() > 3);
                    auto atom_name = field[0];
                    auto atom_num = std::stoi(field[3]);
                    current_it->second.emplace(atom_name, atom_num);
                }
            }
        }
        return ret;
    }

    void match_and_change_residue_name(const std::vector<Residue> &residues, size_t i,
                                       Residue &current_res) {

        static std::unordered_set<std::string> accepted_residue{
                "GUA", "ADE", "THY", "CYT"
        };
        if (accepted_residue.count(current_res.res_name)) {
            if (i == 0) { //5'
                current_res.res_name = std::string("D") + current_res.res_name[0] + std::to_string(5);
            } else if (i == residues.size() - 1) { // 3'
                current_res.res_name = std::string("D") + current_res.res_name[0] + std::to_string(3);
            } else {
                current_res.res_name = std::string("D") + current_res.res_name[0];
            }
        } else {
            std::cerr << "ERROR !! Unexcepted Residue name " << current_res.res_name << '\n';
            exit(EXIT_FAILURE);
        }
    }


    std::vector<Residue> readPDB(std::ifstream &ifs) {
        std::string line;
        std::vector<Residue> residues;
        while (!ifs.eof()) {
            std::getline(ifs, line);
            if (boost::starts_with(line, "ATOM")) {
                auto atom_num = std::stoi(line.substr(6, 5));
                auto atom_name = line.substr(12, 4);
                boost::trim(atom_name);
                auto res_name = line.substr(17, 3);
                auto res_num = std::stoi(line.substr(22, 4));

                auto x = std::stod(line.substr(30, 8));
                auto y = std::stod(line.substr(38, 8));
                auto z = std::stod(line.substr(46, 8));

                if (residues.empty() || residues.back().residue_num != res_num) {
                    residues.emplace_back(Residue{res_name, res_num});
                }
                residues.back().atoms.emplace_back(atom_num, atom_name, x, y, z);
            }
        }
        return residues;
    }
}


void GQuadruplexPdb2gmx::convert() {
    std::string input_pdb = choose_file("Input PDB : ").extension("pdb").isExist(true);
    std::string output_pdb = choose_file("Output PDB : ").extension("pdb").isExist(false);
    std::string gmxrtp = choose_file("GROMACS RTP : ").extension("rtp").isExist(false);

    auto rtp_mapping = parse_rtp(std::ifstream(gmxrtp));

    std::ifstream ifs(input_pdb);

    std::vector<Residue> residues = readPDB(ifs);

    for (std::size_t i = 0; i < residues.size(); i++) {
        auto &current_res = residues[i];
        match_and_change_residue_name(residues, i, current_res);

        boost::sort(current_res.atoms, [&](auto &lhs, auto &rhs) {
            auto &res = rtp_mapping[current_res.res_name];
            auto lhs_it = replace.find(lhs.atom_name);
            auto rhs_it = replace.find(rhs.atom_name);
            auto lhs_num = res[lhs_it != replace.end() ? lhs_it->second : lhs.atom_name];
            auto rhs_num = res[rhs_it != replace.end() ? rhs_it->second : rhs.atom_name];
            return lhs_num < rhs_num;
        });
    }

    std::ofstream ofs(output_pdb);

    boost::copy(residues, std::ostream_iterator<Residue>(ofs));
}

namespace {
    void fill_array(double *x, double *y, double *z, std::vector<Residue> &residues) {
        for (auto &r : residues) {
            for (auto &atom : r.atoms) {
                *x++ = atom.x;
                *y++ = atom.y;
                *z++ = atom.z;
            }
        }
    }

    void move_position(double *x, double *y, double *z, std::vector<Residue> &residues, double mid[3]) {
        for (auto &r : residues) {
            for (auto &atom : r.atoms) {
                atom.x = *x++ + mid[0];
                atom.y = *y++ + mid[1];
                atom.z = *z++ + mid[2];
            }
        }
    }
}

void GQuadruplexPdb2gmx::superpose_and_move() {
    std::string input_pdb = choose_file("Input PDB : ").extension("pdb").isExist(true);
    std::string stub_pdb = choose_file("Input STUB PDB : ").extension("pdb").isExist(false);

    std::ifstream stub_ifs(stub_pdb);
    auto stub = readPDB(stub_ifs);

    std::ifstream pdb_ifs(input_pdb);
    auto residues = readPDB(pdb_ifs);

    double superpose_x[2], superpose_y[2], superpose_z[2];

    auto total_atom_size = boost::accumulate(residues, 0, [](auto result, auto &r) { return result + r.atoms.size(); });
    double x[total_atom_size], y[total_atom_size], z[total_atom_size];

    fill_array(superpose_x, superpose_y, superpose_z, stub);
    fill_array(x, y, z, residues);

    double mid[3];
    RMSDCal::center(total_atom_size, x, y, z, mid, 2);
    RMSDCal::center(2, superpose_x, superpose_y, superpose_z, mid, 2);

    RMSDCal::quatfit(2, superpose_x, superpose_y, superpose_z, total_atom_size, x, y, z, 2);

    move_position(x, y, z, residues, mid);

//    int start_seq = choose(1, "Enter Start Sequence Number for Atom : ");
//    int start_residue_number = choose(1, "Enter Start  Number for Residue : ");
//
//    for (auto &r : residues) {
//        r.residue_num = start_residue_number++;
//        for (auto &atom : r.atoms) {
//            atom.atom_num = start_seq++;
//        }
//    }

    std::string output_pdb = choose_file("Output PDB : ").extension("pdb").isExist(false);

    std::ofstream ofs(output_pdb);

    boost::copy(residues, std::ostream_iterator<Residue>(ofs));
}

