
#ifndef TINKER_MACRORULES_HPP
#define TINKER_MACRORULES_HPP

#include "dsl/AmberMask.hpp"

namespace AMBERMASK {

using std::make_shared;

inline static auto protein = AmberMaskAST::select_ranges{
    "ABU",   "ACE",   "AIB",   "ALA",   "ARG",   "ARGN",  "ASN",   "ASN1",  "ASP",  "ASP1",  "ASPH",  "ASPP",  "ASH",
    "CT3",   "CYS",   "CYS1",  "CYS2",  "CYSH",  "DALA",  "GLN",   "GLU",   "GLUH", "GLUP",  "GLH",   "GLY",   "HIS",
    "HIS1",  "HISA",  "HISB",  "HISH",  "HISD",  "HISE",  "HISP",  "HSD",   "HSE",  "HSP",   "HYP",   "ILE",   "LEU",
    "LSN",   "LYS",   "LYSH",  "MELEU", "MET",   "MEVAL", "NAC",   "NME",   "NHE",  "NH2",   "PHE",   "PHEH",  "PHEU",
    "PHL",   "PRO",   "SER",   "THR",   "TRP",   "TRPH",  "TRPU",  "TYR",   "TYRH", "TYRU",  "VAL",   "PGLU",  "HID",
    "HIE",   "HIP",   "LYP",   "LYN",   "CYN",   "CYM",   "CYX",   "DAB",   "ORN",  "NALA",  "NGLY",  "NSER",  "NTHR",
    "NLEU",  "NILE",  "NVAL",  "NASN",  "NGLN",  "NARG",  "NHID",  "NHIE",  "NHIP", "NHISD", "NHISE", "NHISH", "NTRP",
    "NPHE",  "NTYR",  "NGLU",  "NASP",  "NLYS",  "NORN",  "NDAB",  "NLYSN", "NPRO", "NHYP",  "NCYS",  "NCYS2", "NMET",
    "NASPH", "NGLUH", "CALA",  "CGLY",  "CSER",  "CTHR",  "CLEU",  "CILE",  "CVAL", "CASN",  "CGLN",  "CARG",  "CHID",
    "CHIE",  "CHIP",  "CHISD", "CHISE", "CHISH", "CTRP",  "CPHE",  "CTYR",  "CGLU", "CASP",  "CLYS",  "CORN",  "CDAB",
    "CLYSN", "CPRO",  "CHYP",  "CCYS",  "CCYS2", "CMET",  "CASPH", "CGLUH"};

inline static auto dna = AmberMaskAST::select_ranges{"DA5", "DA", "DA3", "DAN", "DT5", "DT", "DT3", "DTN",
                                                     "DG5", "DG", "DG3", "DGN", "DC5", "DC", "DC3", "DCN"};

inline static auto rna =
    AmberMaskAST::select_ranges{"A",   "U",  "C",   "G",   "RA5", "RA", "RA3", "RAN", "RU5", "RU",  "RU3", "RUN",
                                "RG5", "RG", "RG3", "RGN", "RC5", "RC", "RC3", "RCN", "RT5", "RT3", "RTN"};

inline static auto water = AmberMaskAST::select_ranges{"SOL", "WAT", "HOH", "OHH", "TIP", "T3P", "T4P", "T5P", "T3H"};

inline static AmberMask System = make_shared<AmberMaskAST::atom_name_nums>(AmberMaskAST::select_ranges{"*"});

inline static AmberMask Protein = make_shared<AmberMaskAST::residue_name_nums>(protein);

inline static AmberMask Not_Hydrogen = make_shared<AmberMaskAST::Operator>(
    AmberMaskAST::Op::NOT, make_shared<AmberMaskAST::atom_name_nums>(AmberMaskAST::select_ranges{"H*"}));

inline static AmberMask Protein_H = make_shared<AmberMaskAST::Operator>(AmberMaskAST::Op::AND, Protein, Not_Hydrogen);

inline static AmberMask Backbone = make_shared<AmberMaskAST::Operator>(
    AmberMaskAST::Op::AND, Protein,
    make_shared<AmberMaskAST::atom_name_nums>(AmberMaskAST::select_ranges{"CA", "C", "N"}));

inline static AmberMask MainChain = make_shared<AmberMaskAST::Operator>(
    AmberMaskAST::Op::AND, Protein,
    make_shared<AmberMaskAST::atom_name_nums>(AmberMaskAST::select_ranges{"CA", "C", "N", "O"}));

inline static AmberMask MainChain_plus_Cb = make_shared<AmberMaskAST::Operator>(
    AmberMaskAST::Op::AND, Protein,
    make_shared<AmberMaskAST::atom_name_nums>(AmberMaskAST::select_ranges{"CA", "C", "N", "O", "CB"}));

inline static AmberMask MainChain_plus_H = make_shared<AmberMaskAST::Operator>(
    AmberMaskAST::Op::AND, Protein,
    make_shared<AmberMaskAST::atom_name_nums>(AmberMaskAST::select_ranges{"CA", "C", "N", "O", "H"}));

inline static AmberMask C_alpha = make_shared<AmberMaskAST::Operator>(
    AmberMaskAST::Op::AND, Protein, make_shared<AmberMaskAST::atom_name_nums>(AmberMaskAST::select_ranges{"CA"}));

inline static AmberMask SideChain = make_shared<AmberMaskAST::Operator>(
    AmberMaskAST::Op::AND, Protein,
    make_shared<AmberMaskAST::Operator>(
        AmberMaskAST::Op::NOT,
        make_shared<AmberMaskAST::atom_name_nums>(AmberMaskAST::select_ranges{"CA", "C", "O", "N", "H"})));

inline static AmberMask SideChain_H =
    make_shared<AmberMaskAST::Operator>(AmberMaskAST::Op::AND, SideChain, Not_Hydrogen);

inline static AmberMask DNA = make_shared<AmberMaskAST::residue_name_nums>(dna);

inline static AmberMask DNA_H = make_shared<AmberMaskAST::Operator>(AmberMaskAST::Op::AND, DNA, Not_Hydrogen);

inline static AmberMask RNA = make_shared<AmberMaskAST::residue_name_nums>(rna);

inline static AmberMask RNA_H = make_shared<AmberMaskAST::Operator>(AmberMaskAST::Op::AND, RNA, Not_Hydrogen);

inline static AmberMask Water = make_shared<AmberMaskAST::residue_name_nums>(water);

} // namespace AMBERMASK

#endif // TINKER_MACRORULES_HPP