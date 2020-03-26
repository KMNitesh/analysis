
#include "TprReader.hpp"

#include "data_structure/atom.hpp"
#include "data_structure/frame.hpp"
#include "data_structure/molecule.hpp"
#include "topology_utils.hpp"
#include "utils/common.hpp"

namespace gmx {

#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trnio.h"
#include "gromacs/fileio/xtcio.h"
#include "gromacs/legacyheaders/types/enums.h"
#include "gromacs/legacyheaders/types/ifunc.h"
#include "gromacs/legacyheaders/types/inputrec.h"
#include "gromacs/legacyheaders/types/simple.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/smalloc.h"

} // namespace gmx

std::shared_ptr<Frame> TprReader::read(const std::string &filename) {
    auto frame = std::make_shared<Frame>();

    gmx::t_state state;
    gmx::t_tpxheader tpx;
    gmx::t_inputrec ir;
    gmx::gmx_mtop_t mtop;
    gmx::t_topology top;

    gmx::read_tpxheader(filename.c_str(), &tpx, TRUE, nullptr, nullptr);

    gmx::read_tpx_state(filename.c_str(), tpx.bIr ? &ir : nullptr, &state, nullptr, tpx.bTop ? &mtop : nullptr);

    top = gmx::gmx_mtop_t_to_t_topology(&mtop);

    auto &atoms = top.atoms;

    if (!tpx.bX) {
        std::cerr << "Tpr Topolgy do not have coordinates\n";
        exit(EXIT_FAILURE);
    }

    frame->title = *top.name;
    frame->box = PBCBox(state.box);

    const auto &ffparams = mtop.ffparams;

    for (int i = 0; i < atoms.nr; i++) {
        auto atom = std::make_shared<Atom>();
        atom->seq = i + 1;
        atom->atom_name = (*(atoms.atomname[i]));
        atom->type_name = (*(atoms.atomtype[i]));
        atom->charge = atoms.atom[i].q;
        atom->mass = atoms.atom[i].m;
        atom->typ = atoms.atom[i].type;

        const auto &lj = ffparams.iparams[atom->typ * (ffparams.atnr + 1)].lj;
        atom->lj_param = lj_t{lj.c6, lj.c12};

        atom->x = 10 * state.x[i][XX];
        atom->y = 10 * state.x[i][YY];
        atom->z = 10 * state.x[i][ZZ];

        atom->residue_name = *atoms.resinfo[atoms.atom[i].resind].name;
        atom->residue_num = atoms.resinfo[atoms.atom[i].resind].nr;
        frame->atom_list.push_back(atom);
        frame->atom_map[atom->seq] = atom;
    }

    for (auto t :
         {gmx::F_BONDS, gmx::F_G96BONDS, gmx::F_MORSE, gmx::F_CUBICBONDS, gmx::F_CONNBONDS, gmx::F_HARMONIC,
          gmx::F_FENEBONDS, gmx::F_TABBONDS, gmx::F_TABBONDSNC, gmx::F_RESTRBONDS, gmx::F_CONSTR, gmx::F_CONSTRNC}) {
        auto ilist = &top.idef.il[t];
        gmx::t_iatom *iatoms = ilist->iatoms;
        for (int i = 0; i < ilist->nr; i += 3) {
            auto type = *(iatoms++);
            auto ftype = top.idef.functype[type];
            auto nratoms = gmx::interaction_function[ftype].nratoms;

            if (nratoms != 2) {
                std::cerr << "Code inpsect > " << __FILE__ << ':' << __LINE__ << ' ' << __PRETTY_FUNCTION__
                          << std::endl;
                std::exit(EXIT_FAILURE);
            }

            int atom_num1 = *(iatoms++) + 1;
            int atom_num2 = *(iatoms++) + 1;

            auto atom1 = frame->atom_map[atom_num1];
            auto atom2 = frame->atom_map[atom_num2];

            atom1->con_list.push_back(atom_num2);
            atom2->con_list.push_back(atom_num1);
        }
    }
    auto ilist = &top.idef.il[gmx::F_SETTLE];
    auto iatoms = ilist->iatoms;
    for (int i = 0; i < ilist->nr;) {
        auto type = *(iatoms++);
        auto ftype = top.idef.functype[type];

        auto nratoms = gmx::interaction_function[ftype].nratoms;

        int atom_num1 = *(iatoms++) + 1;
        auto atom1 = frame->atom_map[atom_num1];

        for (int j = 1; j < nratoms; j++) {
            int atom_num2 = *(iatoms++) + 1;
            auto atom2 = frame->atom_map[atom_num2];

            atom1->con_list.push_back(atom_num2);
            atom2->con_list.push_back(atom_num1);
        }

        i += 1 + nratoms;
    }

    for (int k = 0; k < top.idef.il[gmx::F_BONDS].nr;) {
        auto type = top.idef.il[gmx::F_BONDS].iatoms[k++];
        auto ai = top.idef.il[gmx::F_BONDS].iatoms[k++];
        auto aj = top.idef.il[gmx::F_BONDS].iatoms[k++];
        const auto &harmonic = top.idef.iparams[type].harmonic;
        frame->f_bond_params.emplace(std::array{ai + 1, aj + 1},
                                     Frame::harmonic{harmonic.krA * 0.01, harmonic.rA * 10});
    }
    for (int k = 0; k < top.idef.il[gmx::F_ANGLES].nr;) {
        auto type = top.idef.il[gmx::F_ANGLES].iatoms[k++];
        auto ai = top.idef.il[gmx::F_ANGLES].iatoms[k++];
        auto aj = top.idef.il[gmx::F_ANGLES].iatoms[k++];
        auto ak = top.idef.il[gmx::F_ANGLES].iatoms[k++];
        const auto &harmonic = top.idef.iparams[type].harmonic;
        frame->f_angle_params.emplace(std::array{ai + 1, aj + 1, ak + 1}, Frame::harmonic{harmonic.krA, harmonic.rA});
    }
    for (int k = 0; k < top.idef.il[gmx::F_PDIHS].nr;) {
        auto type = top.idef.il[gmx::F_PDIHS].iatoms[k++];
        auto ai = top.idef.il[gmx::F_PDIHS].iatoms[k++];
        auto aj = top.idef.il[gmx::F_PDIHS].iatoms[k++];
        auto ak = top.idef.il[gmx::F_PDIHS].iatoms[k++];
        auto an = top.idef.il[gmx::F_PDIHS].iatoms[k++];

        const auto &pdihs = top.idef.iparams[type].pdihs;
        frame->f_dihedral_params.emplace(std::array{ai + 1, aj + 1, ak + 1, an + 1},
                                         Frame::pdihs{pdihs.phiA, pdihs.cpA, pdihs.mult});
    }

    topology_utils::assgin_atom_to_molecule(frame);
    frame->build_graph();

    if (atoms.nres != boost::numeric_cast<int>(frame->molecule_list.size())) {
        std::cout << "Residue numbers and Molecule numbers not match\n";
        //        exit(EXIT_FAILURE);
    }

    return frame;
}
