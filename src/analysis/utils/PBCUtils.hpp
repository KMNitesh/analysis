//
// Created by xiamr on 6/29/19.
//

#ifndef TINKER_PBCUTILS_HPP
#define TINKER_PBCUTILS_HPP

#include <memory>
#include "ana_module/Trajconv.hpp"

class Frame;

class PBCUtils {

public:
    void do_move_center_basedto_atom(AmberMask &mask, std::shared_ptr<Frame> &frame) const;

    void do_move_center_basedto_molecule(AmberMask &mask, std::shared_ptr<Frame> &frame) const;

    void do_molecule_aggregate(std::shared_ptr<Frame> &frame) const;

    void doPBC(Trajconv::PBCType pbc_mode, AmberMask &mask, std::shared_ptr<Frame> &frame) const;

    static std::shared_ptr<Atom> find_atom(AmberMask &mask, std::shared_ptr<Frame> &frame);

    static std::shared_ptr<Molecule> find_molecule(AmberMask &mask, std::shared_ptr<Frame> &frame);
};


#endif //TINKER_PBCUTILS_HPP
