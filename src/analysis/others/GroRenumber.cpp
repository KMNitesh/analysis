#include "GroRenumber.hpp"

#include <Eigen/Eigen>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>

#include "utils/common.hpp"
namespace gmx {

#include "gromacs/math/vectypes.h"

} // namespace gmx

namespace {
struct GROStruct {
    std::string title;
    int atom_num;

    struct Atom {
        std::string name;
        Eigen::Vector3d coord;
    };

    struct Residue {
        Residue(int num, std::string name) : num(num), name(std::move(name)) {}
        int num;
        std::string name;
        std::vector<Atom> atoms;
    };

    std::vector<Residue> residues;
    gmx::matrix box;

    static GROStruct read_gro(std::istream &is);
};

GROStruct GROStruct::read_gro(std::istream &is) {
    GROStruct gro;
    std::getline(is, gro.title);
    boost::trim(gro.title);
    std::string line;
    std::getline(is, line);
    gro.atom_num = std::stoi(line);

    int lengths[] = {5, 5, 5, 5, 8, 8, 8};
    boost::offset_separator ofs(std::begin(lengths), std::end(lengths));
    for (int i = 0; i < gro.atom_num; ++i) {
        //   %5d%-5s%5s%5d%8.3f%8.3f%8.3f
        std::getline(is, line);
        boost::tokenizer tokenizer(line, ofs);
        auto it = std::begin(tokenizer);
        auto residue_num = std::stoi(*it++);
        std::string residue_name = *it++;
        boost::trim(residue_name);

        GROStruct::Atom atom;
        atom.name = *it++;
        boost::trim(atom.name);
        ++it;
        atom.coord[0] = std::stod(*it++);
        atom.coord[1] = std::stod(*it++);
        atom.coord[2] = std::stod(*it++);

        if (gro.residues.empty() or
            (gro.residues.back().name != residue_name or gro.residues.back().num != residue_num)) {
            gro.residues.emplace_back(residue_num, std::move(residue_name));
        }
        gro.residues.back().atoms.push_back(std::move(atom));
    }
    std::getline(is, line);
    auto fields = split(line);
    std::memset(gro.box, 0, 9 * sizeof(gmx::real));

    gro.box[0][0] = std::stod(fields[0]);
    gro.box[1][1] = std::stod(fields[1]);
    gro.box[2][2] = std::stod(fields[2]);

    // v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y);
    if (fields.size() > 3) {
        gro.box[0][1] = std::stod(fields[3]);
        gro.box[0][2] = std::stod(fields[4]);
        gro.box[1][0] = std::stod(fields[5]);
        gro.box[1][2] = std::stod(fields[6]);
        gro.box[2][0] = std::stod(fields[7]);
        gro.box[2][1] = std::stod(fields[8]);
    }
    return gro;
}

std::istream &operator>>(std::istream &is, GROStruct &gro) {
    gro = GROStruct::read_gro(is);
    return is;
}

std::ostream &operator<<(std::ostream &os, const GROStruct &gro) {
    os << gro.title << '\n';
    os << gro.atom_num << '\n';
    const boost::format fmt("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n");
    int atom_no = 0;
    // int residue_num = gro.residues.front().num;
    for (const auto &residue : gro.residues) {
        for (const auto &atom : residue.atoms) {
            os << boost::format(fmt) % (residue.num % 100000) % residue.name % atom.name % (++atom_no % 100000) %
                      atom.coord[0] % atom.coord[1] % atom.coord[2];
        }
        // ++residue_num;
    }
    os << boost::format(" %9.5f %9.5f %9.5f\n") % gro.box[0][0] % gro.box[1][1] % gro.box[2][2];
    return os;
}

} // namespace

void GroRenumber::process() {
    std::string gro_filename = choose_file("gro file : ").extension("gro").isExist(true);
    std::string output_gro_filename = choose_file("output gro file : ").extension("gro").isExist(false);

    GROStruct gro;
    std::ifstream ifs(gro_filename);

    ifs >> gro;

    std::ofstream ofs(output_gro_filename);

    ofs << gro;
}
