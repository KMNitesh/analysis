#ifndef TINKER_PRMTOPPARSER_HPP
#define TINKER_PRMTOPPARSER_HPP

#include <boost/fusion/sequence.hpp>
#include <boost/optional.hpp>

#include "utils/std.hpp"

struct PrmtopStruct final {
    PrmtopStruct() = default;

    PrmtopStruct(const PrmtopStruct &) = delete;

    PrmtopStruct &operator=(const PrmtopStruct &) = delete;

    PrmtopStruct(PrmtopStruct &&) = default;

    PrmtopStruct &operator=(PrmtopStruct &&) = default;

    std::string version;
    std::string title;

    std::vector<int> pointers;
    enum Pointer {
        NATOM,     //  total number of atoms
        NTYPES,    // total number of distinct atom types
        NBONH,     // number of bonds containing hydrogen
        MBONA,     // number of bonds not containing hydrogen
        NTHETH,    // number of angles containing hydrogen
        MTHETA,    // number of angles not containing hydrogen
        NPHIH,     // number of dihedrals containing hydrogen
        MPHIA,     // number of dihedrals not containing hydrogen
        NHPARM,    // currently not used
        NPARM,     // used to determine if addles created prmtop
        NNB,       // number of excluded atoms
        NRES,      // number of residues
        NBONA,     // MBONA + number of constraint bonds
        NTHETA,    // MTHETA + number of constraint angles
        NPHIA,     // MPHIA + number of constraint dihedrals
        NUMBND,    // number of unique bond types
        NUMANG,    // number of unique angle types
        NPTRA,     // number of unique dihedral types
        NATYP,     // number of atom types in parameter file, see SOLTY below
        NPHB,      // number of distinct 10-12 hydrogen bond pair types
        IFPERT,    // set to 1 if perturbation info is to be read in
        NBPER,     // number of bonds to be perturbed
        NGPER,     // number of angles to be perturbed
        NDPER,     // number of dihedrals to be perturbed
        MBPER,     // number of bonds with atoms completely in perturbed group
        MGPER,     // number of angles with atoms completely in perturbed group
        MDPER,     // number of dihedrals with atoms completely in perturbed groups
        IFBOX,     // set to 1 if standard periodic box, 2 when truncated octahedral
        NMXRS,     // number of atoms in the largest residue
        IFCAP,     // set to 1 if the CAP option from edit was specified
        NUMEXTRA,  // number of extra points found in topology
        NCOPY      // number of PIMD slices / number of beads
    };

    std::vector<std::string> atom_name;
    std::vector<double> charge;
    std::vector<int> atomic_number;
    std::vector<double> mass;
    std::vector<int> atom_type_index;
    std::vector<int> number_excluded_atoms;
    std::vector<int> nonbonded_parm_index;
    std::vector<std::string> residue_label;
    std::vector<int> residue_pointer;
    std::vector<double> bond_force_constant;
    std::vector<double> bond_equil_value;
    std::vector<double> angle_force_constant;
    std::vector<double> angle_equil_value;
    std::vector<double> dihedral_force_constant;
    std::vector<double> dihedral_periodicity;
    std::vector<double> dihedral_phase;
    std::vector<double> scee_scale_factor;
    std::vector<double> scnb_scale_factor;
    std::vector<double> solty;
    std::vector<double> lennard_jones_acoef;
    std::vector<double> lennard_jones_bcoef;
    std::vector<int> bonds_inc_hydrogen;
    std::vector<int> bonds_without_hydrogen;
    std::vector<int> angles_inc_hydrogen;
    std::vector<int> angles_without_hydrogen;
    std::vector<int> dihedrals_inc_hydrogen;
    std::vector<int> dihedrals_without_hydrogen;
    std::vector<int> excluded_atoms_list;
    std::vector<double> hbond_acoef;
    std::vector<double> hbond_bcoef;
    std::vector<std::string> amber_atom_type;
    std::vector<std::string> tree_chain_classification;
    std::vector<int> join_array;
    std::vector<int> irotat;
    boost::optional<boost::fusion::vector<std::vector<int>, std::vector<int>, std::vector<double>>>
        solvent_pointers_atoms_per_molecule_box_dimensions;
    std::string radius_set;
    std::vector<double> radii;
    std::vector<double> screen;

    int ipol;
};

class PrmtopParser {
public:
    static boost::optional<PrmtopStruct> parse(std::istream &is);
};

#endif  // TINKER_PRMTOPPARSER_HPP
