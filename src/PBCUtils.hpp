//
// Created by xiamr on 6/29/19.
//

#ifndef TINKER_PBCUTILS_HPP
#define TINKER_PBCUTILS_HPP

#include <memory>
#include "Trajconv.hpp"

class Frame;

class PBCUtils {

public:
    void do_move_center_basedto_atom(int num, std::shared_ptr<Frame> &frame) const;

    void do_move_center_basedto_molecule(int num, std::shared_ptr<Frame> &frame) const;

    void do_molecule_aggregate(std::shared_ptr<Frame> &frame) const;

    void doPBC(Trajconv::PBCType pbc_mode, int num, std::shared_ptr<Frame> &frame) const;
};


#endif //TINKER_PBCUTILS_HPP
