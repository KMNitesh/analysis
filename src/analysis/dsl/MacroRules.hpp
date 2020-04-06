
#ifndef TINKER_MACRORULES_HPP
#define TINKER_MACRORULES_HPP

#include "dsl/AmberMask.hpp"

namespace AMBERMASK {

extern AmberMaskAST::select_ranges protein, dna, rna, water;

extern AmberMask System, Protein, Protein_H, Backbone, MainChain, MainChain_plus_Cb, MainChain_plus_H, C_alpha,
    SideChain, SideChain_H, DNA, DNA_H, RNA, RNA_H, Water;

} // namespace AMBERMASK

#endif // TINKER_MACRORULES_HPP