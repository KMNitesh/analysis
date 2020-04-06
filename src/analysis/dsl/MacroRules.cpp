
#include "dsl/MacroRules.hpp"

namespace AMBERMASK {

using std::make_shared;
using namespace AmberMaskAST;

select_ranges protein{
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

select_ranges dna{"DA5", "DA", "DA3", "DAN", "DT5", "DT", "DT3", "DTN",
                  "DG5", "DG", "DG3", "DGN", "DC5", "DC", "DC3", "DCN"};

select_ranges rna{"A",   "U",  "C",   "G",   "RA5", "RA", "RA3", "RAN", "RU5", "RU",  "RU3", "RUN",
                  "RG5", "RG", "RG3", "RGN", "RC5", "RC", "RC3", "RCN", "RT5", "RT3", "RTN"};
select_ranges water{"SOL", "WAT", "HOH", "OHH", "TIP", "T3P", "T4P", "T5P", "T3H"};

AmberMask System = make_shared<atom_name_nums>(select_ranges{"*"});

AmberMask Protein = make_shared<residue_name_nums>(protein);

AmberMask Not_Hydrogen = make_shared<Operator>(Op::NOT, make_shared<atom_name_nums>(select_ranges{"H*"}));

AmberMask Protein_H = make_shared<Operator>(Op::AND, Protein, Not_Hydrogen);

AmberMask Backbone =
    make_shared<Operator>(Op::AND, Protein, make_shared<atom_name_nums>(select_ranges{"CA", "C", "N"}));

AmberMask MainChain = make_shared<Operator>(
    Op::AND, Protein, make_shared<atom_name_nums>(select_ranges{"CA", "C", "N", "O", "OC1", "OC2"}));

AmberMask MainChain_plus_Cb = make_shared<Operator>(
    Op::AND, Protein, make_shared<atom_name_nums>(select_ranges{"CA", "C", "N", "O", "OC1", "OC2", "CB"}));

AmberMask MainChain_plus_H = make_shared<Operator>(
    Op::AND, Protein,
    make_shared<atom_name_nums>(select_ranges{"CA", "C", "N", "O", "OC1", "OC2", "H", "H1", "H2", "H3"}));

AmberMask C_alpha = make_shared<Operator>(Op::AND, Protein, make_shared<atom_name_nums>(select_ranges{"CA"}));

AmberMask SideChain = make_shared<Operator>(
    Op::AND, Protein,
    make_shared<Operator>(
        Op::NOT, make_shared<atom_name_nums>(select_ranges{"CA", "C", "O", "N", "H", "OC1", "OC2", "H1", "H2", "H3"})));

AmberMask SideChain_H = make_shared<Operator>(Op::AND, SideChain, Not_Hydrogen);

AmberMask DNA = make_shared<residue_name_nums>(dna);

AmberMask DNA_H = make_shared<Operator>(Op::AND, DNA, Not_Hydrogen);

AmberMask RNA = make_shared<residue_name_nums>(rna);

AmberMask RNA_H = make_shared<Operator>(Op::AND, RNA, Not_Hydrogen);

AmberMask Water = make_shared<residue_name_nums>(water);

} // namespace AMBERMASK
