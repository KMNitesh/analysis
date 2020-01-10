#ifndef TINKER_PRMTOPGRAMMAR_HPP
#define TINKER_PRMTOPGRAMMAR_HPP

#include <boost/algorithm/string.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/fusion/include/io.hpp>
#include <boost/fusion/adapted/adt/adapt_adt.hpp>
#include <boost/fusion/include/adapt_adt.hpp>
#include <boost/phoenix.hpp>
#include "trajectory_reader/PrmtopParser.hpp"


BOOST_PHOENIX_ADAPT_FUNCTION(void, trim, boost::trim, 1)

BOOST_PHOENIX_ADAPT_FUNCTION(int, stoi, std::stoi, 1)

struct field_f {
    template<typename... A>
    struct result {
        typedef int type;
    };

    int operator()(const PrmtopStruct &v, int i) const {
        return v.pointers[i];
    }
};


template<typename... _Args>
inline auto field_(_Args &&... __args) {
    return boost::phoenix::function<field_f>()(std::forward<_Args>(__args)...);
}

using Optional = boost::optional<boost::fusion::vector<std::vector<int>, std::vector<int>, std::vector<double>>>;

BOOST_FUSION_ADAPT_STRUCT(PrmtopStruct, (std::string, version)
        (std::string, title)
        (std::vector<int>, pointers)
        (std::vector<std::string>, atom_name)
        (std::vector<double>, charge)
        (std::vector<int>, atomic_number)
        (std::vector<double>, mass)
        (std::vector<int>, atom_type_index)
        (std::vector<int>, number_excluded_atoms)
        (std::vector<int>, nonbonded_parm_index)
        (std::vector<std::string>, residue_label)
        (std::vector<int>, residue_pointer)
        (std::vector<double>, bond_force_constant)
        (std::vector<double>, bond_equil_value)
        (std::vector<double>, angle_force_constant)
        (std::vector<double>, angle_equil_value)
        (std::vector<double>, dihedral_force_constant)
        (std::vector<double>, dihedral_periodicity)
        (std::vector<double>, dihedral_phase)
        (std::vector<double>, scee_scale_factor)
        (std::vector<double>, scnb_scale_factor)
        (std::vector<double>, solty)
        (std::vector<double>, lennard_jones_acoef)
        (std::vector<double>, lennard_jones_bcoef)
        (std::vector<int>, bonds_inc_hydrogen)
        (std::vector<int>, bonds_without_hydrogen)
        (std::vector<int>, angles_inc_hydrogen)
        (std::vector<int>, angles_without_hydrogen)
        (std::vector<int>, dihedrals_inc_hydrogen)
        (std::vector<int>, dihedrals_without_hydrogen)
        (std::vector<int>, excluded_atoms_list)
        (std::vector<double>, hbond_acoef)
        (std::vector<double>, hbond_bcoef)
        (std::vector<std::string>, amber_atom_type)
        (std::vector<std::string>, tree_chain_classification)
        (std::vector<int>, join_array)
        (std::vector<int>, irotat)
        (Optional, solvent_pointers_atoms_per_molecule_box_dimensions)
        (std::string, radius_set)
        (std::vector<double>, radii)
        (std::vector<double>, screen)
        (int, ipol)
)


template<typename Iterator, typename Skipper>
struct PrmtopGrammar : boost::spirit::qi::grammar<Iterator, PrmtopStruct(), Skipper> {
    PrmtopGrammar();

    boost::spirit::qi::rule<Iterator, std::string(), Skipper> version;
    boost::spirit::qi::rule<Iterator, std::string(), Skipper> title;
    boost::spirit::qi::rule<Iterator, std::vector<int>(), Skipper> pointers;
    boost::spirit::qi::rule<Iterator, int(int), boost::spirit::ascii::space_type> fix_int;
    boost::spirit::qi::rule<Iterator, std::string(int), boost::spirit::ascii::space_type> string_;
    boost::spirit::qi::rule<Iterator, std::vector<std::string>(int), Skipper> atom_name;
    boost::spirit::qi::rule<Iterator, std::vector<double>(int), Skipper> charge;
    boost::spirit::qi::rule<Iterator, std::vector<int>(int), Skipper> atomic_number;
    boost::spirit::qi::rule<Iterator, std::vector<double>(int), Skipper> mass;
    boost::spirit::qi::rule<Iterator, std::vector<int>(int), Skipper> atom_type_index;
    boost::spirit::qi::rule<Iterator, std::vector<int>(int), Skipper> number_excluded_atoms;
    boost::spirit::qi::rule<Iterator, std::vector<int>(int), Skipper> nonbonded_parm_index;
    boost::spirit::qi::rule<Iterator, std::vector<std::string>(int), Skipper> residue_label;
    boost::spirit::qi::rule<Iterator, std::vector<int>(int), Skipper> residue_pointer;
    boost::spirit::qi::rule<Iterator, std::vector<double>(int), Skipper> bond_force_constant;
    boost::spirit::qi::rule<Iterator, std::vector<double>(int), Skipper> bond_equil_value;
    boost::spirit::qi::rule<Iterator, std::vector<double>(int), Skipper> angle_force_constant;
    boost::spirit::qi::rule<Iterator, std::vector<double>(int), Skipper> angle_equil_value;
    boost::spirit::qi::rule<Iterator, std::vector<double>(int), Skipper> dihedral_force_constant;
    boost::spirit::qi::rule<Iterator, std::vector<double>(int), Skipper> dihedral_periodicity;
    boost::spirit::qi::rule<Iterator, std::vector<double>(int), Skipper> dihedral_phase;
    boost::spirit::qi::rule<Iterator, std::vector<double>(int), Skipper> scee_scale_factor;
    boost::spirit::qi::rule<Iterator, std::vector<double>(int), Skipper> scnb_scale_factor;
    boost::spirit::qi::rule<Iterator, std::vector<double>(int), Skipper> solty;
    boost::spirit::qi::rule<Iterator, std::vector<double>(int), Skipper> lennard_jones_acoef;
    boost::spirit::qi::rule<Iterator, std::vector<double>(int), Skipper> lennard_jones_bcoef;
    boost::spirit::qi::rule<Iterator, std::vector<int>(int), Skipper> bonds_inc_hydrogen;
    boost::spirit::qi::rule<Iterator, std::vector<int>(int), Skipper> bonds_without_hydrogen;
    boost::spirit::qi::rule<Iterator, std::vector<int>(int), Skipper> angles_inc_hydrogen;
    boost::spirit::qi::rule<Iterator, std::vector<int>(int), Skipper> angles_without_hydrogen;
    boost::spirit::qi::rule<Iterator, std::vector<int>(int), Skipper> dihedrals_inc_hydrogen;
    boost::spirit::qi::rule<Iterator, std::vector<int>(int), Skipper> dihedrals_without_hydrogen;
    boost::spirit::qi::rule<Iterator, std::vector<int>(int), Skipper> excluded_atoms_list;
    boost::spirit::qi::rule<Iterator, std::vector<double>(int), Skipper> hbond_acoef;
    boost::spirit::qi::rule<Iterator, std::vector<double>(int), Skipper> hbond_bcoef;
    boost::spirit::qi::rule<Iterator, void(int), Skipper> hbcut; // hbond_bcoef
    boost::spirit::qi::rule<Iterator, std::vector<std::string>(int), Skipper> amber_atom_type;
    boost::spirit::qi::rule<Iterator, std::vector<std::string>(int), Skipper> tree_chain_classification;
    boost::spirit::qi::rule<Iterator, std::vector<int>(int), Skipper> join_array;
    boost::spirit::qi::rule<Iterator, std::vector<int>(int), Skipper> irotat;

    boost::spirit::qi::rule<Iterator, std::vector<int>(int), Skipper> solvent_pointers;
    boost::spirit::qi::rule<Iterator, std::vector<int>(int), Skipper> atoms_per_molecule;
    boost::spirit::qi::rule<Iterator, std::vector<double>(int), Skipper> box_dimensions;

    boost::spirit::qi::rule<Iterator,
            boost::fusion::vector<std::vector<int>, std::vector<int>, std::vector<double>>(),
            Skipper> solvent_pointers_atoms_per_molecule_box_dimensions;

    boost::spirit::qi::rule<Iterator, std::string(), Skipper> radius_set;
    boost::spirit::qi::rule<Iterator, std::vector<double>(int), Skipper> radii;
    boost::spirit::qi::rule<Iterator, std::vector<double>(int), Skipper> screen;
    boost::spirit::qi::rule<Iterator, int(), Skipper> ipol;

    boost::spirit::qi::rule<Iterator, PrmtopStruct(), Skipper> root;
};


template<typename Iterator, typename Skipper>
inline PrmtopGrammar<Iterator, Skipper>::PrmtopGrammar() : PrmtopGrammar::base_type(root) {
    using namespace boost::spirit;
    using qi::lexeme;
    using qi::char_;
    using qi::eol;
    using qi::repeat;
    using qi::uint_;
    using qi::skip;
    using qi::_1;
    using qi::_val;
    using qi::_r1;
    using qi::as_string;

    using boost::phoenix::ref;
    using boost::phoenix::at_c;
    using boost::phoenix::at;
    using boost::phoenix::bind;

    using Prmtop = PrmtopStruct;

    version %= "%VERSION" >> lexeme[+(char_ - eol)][trim(_1)] >> eol;

    title %= "%FLAG TITLE"
            >> eol
            >> "%FORMAT(20a4)" >> eol
            >> lexeme[+(char_ - eol)][trim(_1)] >> eol;

    fix_int = as_string[no_skip[-eol >> repeat(_r1)[char_]]][_val = stoi(_1)];

    pointers = "%FLAG POINTERS"
            >> eol
            >> "%FORMAT(10I8)" >> eol
            >> skip(ascii::space)[repeat(31)[fix_int(8)]] >> eol;

    string_ %= no_skip[-eol >> repeat(_r1)[char_]][trim(_1)];

#define SECTION(X, Y, Z) X >> eol >> Y >> eol >> skip(ascii::space)[repeat(_r1)[Z]] >> eol

    atom_name = SECTION("%FLAG ATOM_NAME", "%FORMAT(20a4)", string_(4));

    charge = SECTION("%FLAG CHARGE", "%FORMAT(5E16.8)", double_);

    atomic_number = SECTION("%FLAG ATOMIC_NUMBER", "%FORMAT(10I8)", fix_int(8));

    mass = SECTION("%FLAG MASS", "%FORMAT(5E16.8)", double_);

    atom_type_index = SECTION("%FLAG ATOM_TYPE_INDEX", "%FORMAT(10I8)", fix_int(8));

    number_excluded_atoms = SECTION("%FLAG NUMBER_EXCLUDED_ATOMS", "%FORMAT(10I8)", fix_int(8));

    nonbonded_parm_index = SECTION("%FLAG NONBONDED_PARM_INDEX", "%FORMAT(10I8)", fix_int(8));

    residue_label = SECTION("%FLAG RESIDUE_LABEL", "%FORMAT(20a4)", string_(4));

    residue_pointer = SECTION("%FLAG RESIDUE_POINTER", "%FORMAT(10I8)", fix_int(8));

    bond_force_constant = SECTION("%FLAG BOND_FORCE_CONSTANT", "%FORMAT(5E16.8)", double_);

    bond_equil_value = SECTION("%FLAG BOND_EQUIL_VALUE", "%FORMAT(5E16.8)", double_);

    angle_force_constant = SECTION("%FLAG ANGLE_FORCE_CONSTANT", "%FORMAT(5E16.8)", double_);

    angle_equil_value = SECTION("%FLAG ANGLE_EQUIL_VALUE", "%FORMAT(5E16.8)", double_);

    dihedral_force_constant = SECTION("%FLAG DIHEDRAL_FORCE_CONSTANT", "%FORMAT(5E16.8)", double_);

    dihedral_periodicity = SECTION("%FLAG DIHEDRAL_PERIODICITY", "%FORMAT(5E16.8)", double_);

    dihedral_phase = SECTION("%FLAG DIHEDRAL_PHASE", "%FORMAT(5E16.8)", double_);

    scee_scale_factor = SECTION("%FLAG SCEE_SCALE_FACTOR", "%FORMAT(5E16.8)", double_);

    scnb_scale_factor = SECTION("%FLAG SCNB_SCALE_FACTOR", "%FORMAT(5E16.8)", double_);

    solty = SECTION("%FLAG SOLTY", "%FORMAT(5E16.8)", double_);

    lennard_jones_acoef = SECTION("%FLAG LENNARD_JONES_ACOEF", "%FORMAT(5E16.8)", double_);

    lennard_jones_bcoef = SECTION("%FLAG LENNARD_JONES_BCOEF", "%FORMAT(5E16.8)", double_);

    bonds_inc_hydrogen = SECTION("%FLAG BONDS_INC_HYDROGEN", "%FORMAT(10I8)", fix_int(8));

    bonds_without_hydrogen = SECTION("%FLAG BONDS_WITHOUT_HYDROGEN", "%FORMAT(10I8)", fix_int(8));

    angles_inc_hydrogen = SECTION("%FLAG ANGLES_INC_HYDROGEN", "%FORMAT(10I8)", fix_int(8));

    angles_without_hydrogen = SECTION("%FLAG ANGLES_WITHOUT_HYDROGEN", "%FORMAT(10I8)", fix_int(8));

    dihedrals_inc_hydrogen = SECTION("%FLAG DIHEDRALS_INC_HYDROGEN", "%FORMAT(10I8)", fix_int(8));

    dihedrals_without_hydrogen = SECTION("%FLAG DIHEDRALS_WITHOUT_HYDROGEN", "%FORMAT(10I8)", fix_int(8));

    excluded_atoms_list = SECTION("%FLAG EXCLUDED_ATOMS_LIST", "%FORMAT(10I8)", fix_int(8));

    hbond_acoef = SECTION("%FLAG HBOND_ACOEF", "%FORMAT(5E16.8)", double_);

    hbond_bcoef = SECTION("%FLAG HBOND_BCOEF", "%FORMAT(5E16.8)", double_);

    hbcut = omit[SECTION("%FLAG HBCUT", "%FORMAT(5E16.8)", double_)];

    amber_atom_type = SECTION("%FLAG AMBER_ATOM_TYPE", "%FORMAT(20a4)", string_(4));

    tree_chain_classification = SECTION("%FLAG TREE_CHAIN_CLASSIFICATION", "%FORMAT(20a4)", string_(4));

    join_array = SECTION("%FLAG JOIN_ARRAY", "%FORMAT(10I8)", fix_int(8));

    irotat = SECTION("%FLAG IROTAT", "%FORMAT(10I8)", fix_int(8));

    solvent_pointers = SECTION("%FLAG SOLVENT_POINTERS", "%FORMAT(3I8)", fix_int(8));
    atoms_per_molecule = SECTION("%FLAG ATOMS_PER_MOLECULE", "%FORMAT(10I8)", fix_int(8));
    box_dimensions = SECTION("%FLAG BOX_DIMENSIONS", "%FORMAT(5E16.8)", double_);

    solvent_pointers_atoms_per_molecule_box_dimensions %=
            solvent_pointers(3) >> atoms_per_molecule(at(at_c<0>(ref(_val)), 1)) >> box_dimensions(4);

    radius_set = "%FLAG RADIUS_SET" >> eol >> "%FORMAT(1a80)" >> eol >> no_skip[repeat(80)[char_]] >> eol;

    radii = SECTION("%FLAG RADII", "%FORMAT(5E16.8)", double_);

    screen = SECTION("%FLAG SCREEN", "%FORMAT(5E16.8)", double_);

    ipol = "%FLAG IPOL" >> eol >> "%FORMAT(1I8)" >> eol >> int_ >> eol;

    root = version
            >> title >> pointers
            >> atom_name(field_(ref(_val), Prmtop::NATOM))
            >> charge(field_(ref(_val), Prmtop::NATOM))
            >> atomic_number(field_(ref(_val), Prmtop::NATOM))
            >> mass(field_(ref(_val), Prmtop::NATOM))
            >> atom_type_index(field_(ref(_val), Prmtop::NATOM))
            >> number_excluded_atoms(field_(ref(_val), Prmtop::NATOM))
            >> nonbonded_parm_index(field_(ref(_val), Prmtop::NTYPES) * field_(ref(_val), Prmtop::NTYPES))
            >> residue_label(field_(ref(_val), Prmtop::NRES))
            >> residue_pointer(field_(ref(_val), Prmtop::NRES))
            >> bond_force_constant(field_(ref(_val), Prmtop::NUMBND))
            >> bond_equil_value(field_(ref(_val), Prmtop::NUMBND))
            >> angle_force_constant(field_(ref(_val), Prmtop::NUMANG))
            >> angle_equil_value(field_(ref(_val), Prmtop::NUMANG))
            >> dihedral_force_constant(field_(ref(_val), Prmtop::NPTRA))
            >> dihedral_periodicity(field_(ref(_val), Prmtop::NPTRA))
            >> dihedral_phase(field_(ref(_val), Prmtop::NPTRA))
            >> scee_scale_factor(field_(ref(_val), Prmtop::NPTRA))
            >> scnb_scale_factor(field_(ref(_val), Prmtop::NPTRA))
            >> solty(field_(ref(_val), Prmtop::NATYP))
            >> lennard_jones_acoef(field_(ref(_val), Prmtop::NTYPES) * (field_(ref(_val), Prmtop::NTYPES) + 1) / 2)
            >> lennard_jones_bcoef(field_(ref(_val), Prmtop::NTYPES) * (field_(ref(_val), Prmtop::NTYPES) + 1) / 2)
            >> bonds_inc_hydrogen(3 * field_(ref(_val), Prmtop::NBONH))
            >> bonds_without_hydrogen(3 * field_(ref(_val), Prmtop::NBONA))
            >> angles_inc_hydrogen(4 * field_(ref(_val), Prmtop::NTHETH))
            >> angles_without_hydrogen(4 * field_(ref(_val), Prmtop::NTHETA))
            >> dihedrals_inc_hydrogen(5 * field_(ref(_val), Prmtop::NPHIH))
            >> dihedrals_without_hydrogen(5 * field_(ref(_val), Prmtop::NPHIA))
            >> excluded_atoms_list(field_(ref(_val), Prmtop::NNB))
            >> hbond_acoef(field_(ref(_val), Prmtop::NPHB))
            >> hbond_bcoef(field_(ref(_val), Prmtop::NPHB))
            >> hbcut(field_(ref(_val), Prmtop::NPHB))
            >> amber_atom_type(field_(ref(_val), Prmtop::NATOM))
            >> tree_chain_classification(field_(ref(_val), Prmtop::NATOM))
            >> join_array(field_(ref(_val), Prmtop::NATOM))
            >> irotat(field_(ref(_val), Prmtop::NATOM))
            >> (eps(field_(ref(_val), Prmtop::IFBOX) > 0) >> solvent_pointers_atoms_per_molecule_box_dimensions | eps)
            >> radius_set
            >> radii(field_(ref(_val), Prmtop::NATOM))
            >> screen(field_(ref(_val), Prmtop::NATOM))
            >> ipol;
}


#endif //TINKER_PRMTOPGRAMMAR_HPP
