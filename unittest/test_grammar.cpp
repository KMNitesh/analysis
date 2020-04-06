
#include "data_structure/atom.hpp"
#include "data_structure/molecule.hpp"
#include "dsl/grammar.hpp"
#include <gmock/gmock.h>

using namespace testing;

class TestGrammar : public Test {
public:
    TestGrammar() : atom(std::make_shared<Atom>()) {}

    bool result(bool need_suceess = true) {
        auto it = input_string.begin();
        bool ret;
        bool status;
        try {
            status = qi::phrase_parse(it, input_string.end(), grammar, qi::ascii::space, mask);
            ret = status && it == input_string.end() && is_match(atom, mask);
        } catch (...) {
            ret = false;
            status = false;
        }

        ret = need_suceess == ret;
        if (!ret) {
            boost::apply_visitor(AmberMaskAST::print(), mask);
            if (!(status and (it == input_string.end()))) {
                std::cout << "error-pos : " << std::endl;
                std::cout << input_string << std::endl;
                for (auto iter = input_string.begin(); iter != it; ++iter)
                    std::cout << " ";
                std::cout << "^" << std::endl;
            }
        }
        return ret;
    }

    Grammar<std::string::iterator, qi::ascii::space_type> grammar;
    AmberMask mask;
    std::shared_ptr<Atom> atom;
    std::string input_string;
};

TEST_F(TestGrammar, test_residue_name1) {
    input_string = ":ASP";
    atom->residue_name = "ASP";
    atom->residue_num = 1;

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_residue_name2) {
    input_string = ":A*P";
    atom->residue_name = "ASP";
    atom->residue_num = 1;

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_residue_name3) {
    input_string = ":A?P";
    atom->residue_name = "ASP";
    atom->residue_num = 1;

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_residue_name4) {
    input_string = ":A=P";
    atom->residue_name = "ASP";
    atom->residue_num = 1;

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_residue_name5) {
    input_string = ":*SP";
    atom->residue_name = "ASP";
    atom->residue_num = 1;

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_residue_name6) {
    input_string = ":?SP";
    atom->residue_name = "ASP";
    atom->residue_num = 1;

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_residue_name7) {
    input_string = ":=SP";
    atom->residue_name = "ASP";
    atom->residue_num = 1;

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_residue_number1) {
    input_string = ":1";
    atom->residue_name = "ASP";
    atom->residue_num = 1;

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_residue_number2) {
    input_string = ":1*";
    atom->residue_name = "ASP";
    atom->residue_num = 10;

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_residue_number3) {
    input_string = ":1?";
    atom->residue_name = "ASP";
    atom->residue_num = 10;

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_residue_number4) {
    input_string = ":1=";
    atom->residue_name = "ASP";
    atom->residue_num = 10;

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_residue_number5) {
    input_string = ":*0";
    atom->residue_name = "ASP";
    atom->residue_num = 10;

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_residue_number6) {
    input_string = ":?0";
    atom->residue_name = "ASP";
    atom->residue_num = 10;

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_residue_number7) {
    input_string = ":=0";
    atom->residue_name = "ASP";
    atom->residue_num = 10;

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_residue_number_range) {
    input_string = ":1-10";
    atom->residue_name = "ASP";
    atom->residue_num = 5;

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_residue_number_range2) {
    input_string = ":1-10#2";
    atom->residue_name = "ASP";
    atom->residue_num = 5;

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_residue_number_range3) {
    input_string = ":1-10#3";
    atom->residue_name = "ASP";
    atom->residue_num = 5;

    ASSERT_TRUE(result(false));
}

TEST_F(TestGrammar, test_molecule_number1) {
    input_string = "$1";
    auto molecule = std::make_shared<Molecule>();
    atom->molecule = molecule;
    molecule->sequence = 1;

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_molecule_number2) {
    input_string = "$1*";
    auto molecule = std::make_shared<Molecule>();
    atom->molecule = molecule;
    molecule->sequence = 1;

    ASSERT_TRUE(result(false));
}

TEST_F(TestGrammar, test_molecule_number3) {
    input_string = "$1?";
    auto molecule = std::make_shared<Molecule>();
    atom->molecule = molecule;
    molecule->sequence = 1;

    ASSERT_TRUE(result(false));
}

TEST_F(TestGrammar, test_molecule_number4) {
    input_string = "$1=";
    auto molecule = std::make_shared<Molecule>();
    atom->molecule = molecule;
    molecule->sequence = 1;

    ASSERT_TRUE(result(false));
}

TEST_F(TestGrammar, test_molecule_number5) {
    input_string = "$*0";
    auto molecule = std::make_shared<Molecule>();
    atom->molecule = molecule;
    molecule->sequence = 1;

    ASSERT_TRUE(result(false));
}

TEST_F(TestGrammar, test_molecule_number6) {
    input_string = "$?0";
    auto molecule = std::make_shared<Molecule>();
    atom->molecule = molecule;
    molecule->sequence = 1;

    ASSERT_TRUE(result(false));
}

TEST_F(TestGrammar, test_molecule_number7) {
    input_string = "$=0";
    auto molecule = std::make_shared<Molecule>();
    atom->molecule = molecule;
    molecule->sequence = 1;

    ASSERT_TRUE(result(false));
}

TEST_F(TestGrammar, test_molecule_number_range) {
    input_string = "$1-10";
    auto molecule = std::make_shared<Molecule>();
    atom->molecule = molecule;
    molecule->sequence = 1;

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_molecule_number_range2) {
    input_string = "$1-10#2";
    auto molecule = std::make_shared<Molecule>();
    atom->molecule = molecule;
    molecule->sequence = 3;

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_molecule_number_range3) {
    input_string = "$1-10#3";
    auto molecule = std::make_shared<Molecule>();
    atom->molecule = molecule;
    molecule->sequence = 5;

    ASSERT_TRUE(result(false));
}

TEST_F(TestGrammar, test_atom_name) {
    input_string = "@CA";
    atom->atom_name = "CA";

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_atom_name2) {
    input_string = "@C";
    atom->atom_name = "CA";

    ASSERT_TRUE(result(false));
}

TEST_F(TestGrammar, test_atom_number) {
    input_string = "@10";
    atom->seq = 10;

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_atom_number2) {
    input_string = "@1,2,5,9";
    atom->seq = 5;

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_atom_number_seq_add_range) {
    input_string = "@1,2-5#2,9";
    atom->seq = 5;

    ASSERT_TRUE(result(false));
}

TEST_F(TestGrammar, test_atom_number_range) {
    input_string = "@1-10";
    atom->seq = 5;

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_atom_number_range_with_step) {
    input_string = "@1-10#2";
    atom->seq = 5;

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_atom_number_range_with_step2) {
    input_string = "@1-10#3";
    atom->seq = 5;

    ASSERT_TRUE(result(false));
}

TEST_F(TestGrammar, test_atom_type_name) {
    input_string = "@%CA";
    atom->type_name = "CA";

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_atom_type_number) {
    input_string = "@%10";
    atom->typ = 10;

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_atom_type_number_range) {
    input_string = "@%1-10";
    atom->typ = 5;

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_atom_type_number_range_with_step) {
    input_string = "@%1-10#4";
    atom->typ = 5;

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_atom_symbol_name) {
    input_string = "@/C";
    atom->atom_symbol = "C";

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_not) {
    input_string = "!@10";
    atom->seq = 6;

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_and) {
    input_string = ":1 & @10";
    atom->residue_name = "ASP";
    atom->residue_num = 1;
    atom->seq = 10;

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_or) {
    input_string = ":1 | @10";
    atom->residue_name = "ASP";
    atom->residue_num = 2;
    atom->seq = 10;

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_and_plus_or) {
    input_string = ":1 & @C | @%CA";
    atom->residue_name = "ASP";
    atom->residue_num = 2;
    atom->atom_name = "CB";
    atom->type_name = "CA";

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_parentheses1) {
    input_string = ":1 & ( @C | @%CA )";
    atom->residue_name = "ASP";
    atom->residue_num = 1;
    atom->atom_name = "CB";
    atom->type_name = "CA";

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, test_parentheses2) {
    input_string = ":1-10#1 & ( @C | @%CA )";
    atom->residue_name = "ASP";
    atom->residue_num = 9;
    atom->atom_name = "CB";
    atom->type_name = "CA";

    ASSERT_TRUE(result());
}

TEST_F(TestGrammar, zero_step) {
    input_string = ":1-10#0 & ( @C | @%CA )";
    atom->residue_name = "ASP";
    atom->residue_num = 9;
    atom->atom_name = "CB";
    atom->type_name = "CA";

    ASSERT_TRUE(result(false));
}

class TestGrammarMacro : public Test {
protected:
    void parse() {
        auto it = input_string.begin();
        ASSERT_TRUE(qi::phrase_parse(it, input_string.end(), grammar, qi::ascii::space, mask) and
                    it == input_string.end());
    }

    Grammar<std::string::iterator, qi::ascii::space_type> grammar;
    AmberMask mask;
    std::string input_string;
};

TEST_F(TestGrammarMacro, System) {
    input_string = "System";
    parse();
    AmberMask node = make_shared<AmberMaskAST::atom_name_nums>(AmberMaskAST::select_ranges{"*"});
    ASSERT_THAT(mask, Eq(node));
}

TEST_F(TestGrammarMacro, All) {
    input_string = "All";
    parse();
    AmberMask node = make_shared<AmberMaskAST::atom_name_nums>(AmberMaskAST::select_ranges{"*"});
    ASSERT_THAT(mask, Eq(node));
}

TEST_F(TestGrammarMacro, Protein) {
    input_string = "Protein";
    parse();
    AmberMask node = make_shared<AmberMaskAST::residue_name_nums>(AMBERMASK::protein);
    ASSERT_THAT(mask, Eq(node));
}

TEST_F(TestGrammarMacro, Protein_H) {
    input_string = "Protein-H";
    parse();
    AmberMask node = make_shared<AmberMaskAST::Operator>(
        AmberMaskAST::Op::AND, make_shared<AmberMaskAST::residue_name_nums>(AMBERMASK::protein),
        make_shared<AmberMaskAST::Operator>(
            AmberMaskAST::Op::NOT, make_shared<AmberMaskAST::atom_name_nums>(AmberMaskAST::select_ranges{"H*"})));
    ASSERT_THAT(mask, Eq(node));
}

TEST_F(TestGrammarMacro, Backbone) {
    input_string = "Backbone";
    parse();
    AmberMask node = make_shared<AmberMaskAST::Operator>(
        AmberMaskAST::Op::AND, make_shared<AmberMaskAST::residue_name_nums>(AMBERMASK::protein),
        make_shared<AmberMaskAST::atom_name_nums>(AmberMaskAST::select_ranges{"CA", "C", "N"}));
    ASSERT_THAT(mask, Eq(node));
}

TEST_F(TestGrammarMacro, MainChain) {
    input_string = "MainChain";
    parse();
    AmberMask node = make_shared<AmberMaskAST::Operator>(
        AmberMaskAST::Op::AND, make_shared<AmberMaskAST::residue_name_nums>(AMBERMASK::protein),
        make_shared<AmberMaskAST::atom_name_nums>(AmberMaskAST::select_ranges{"CA", "C", "N", "O", "OC1", "OC2"}));
    ASSERT_THAT(mask, Eq(node));
}

TEST_F(TestGrammarMacro, MainChainPlusCb) {
    input_string = "MainChain+Cb";
    parse();
    AmberMask node = make_shared<AmberMaskAST::Operator>(
        AmberMaskAST::Op::AND, make_shared<AmberMaskAST::residue_name_nums>(AMBERMASK::protein),
        make_shared<AmberMaskAST::atom_name_nums>(
            AmberMaskAST::select_ranges{"CA", "C", "N", "O", "OC1", "OC2", "CB"}));
    ASSERT_THAT(mask, Eq(node));
}

TEST_F(TestGrammarMacro, MainChainPlusH) {
    input_string = "MainChain+H";
    parse();
    AmberMask node = make_shared<AmberMaskAST::Operator>(
        AmberMaskAST::Op::AND, make_shared<AmberMaskAST::residue_name_nums>(AMBERMASK::protein),
        make_shared<AmberMaskAST::atom_name_nums>(
            AmberMaskAST::select_ranges{"CA", "C", "N", "O", "OC1", "OC2", "H", "H1", "H2", "H3"}));
    ASSERT_THAT(mask, Eq(node));
}

TEST_F(TestGrammarMacro, C_alpha) {
    input_string = "C-alpha";
    parse();
    AmberMask node = make_shared<AmberMaskAST::Operator>(
        AmberMaskAST::Op::AND, make_shared<AmberMaskAST::residue_name_nums>(AMBERMASK::protein),
        make_shared<AmberMaskAST::atom_name_nums>(AmberMaskAST::select_ranges{"CA"}));
    ASSERT_THAT(mask, Eq(node));
}

TEST_F(TestGrammarMacro, SideChain) {
    input_string = "SideChain";
    parse();
    AmberMask node = make_shared<AmberMaskAST::Operator>(
        AmberMaskAST::Op::AND, make_shared<AmberMaskAST::residue_name_nums>(AMBERMASK::protein),
        make_shared<AmberMaskAST::Operator>(AmberMaskAST::Op::NOT,
                                            make_shared<AmberMaskAST::atom_name_nums>(AmberMaskAST::select_ranges{
                                                "CA", "C", "O", "N", "H", "OC1", "OC2", "H1", "H2", "H3"})));
    ASSERT_THAT(mask, Eq(node));
}

TEST_F(TestGrammarMacro, SideChain_H) {
    input_string = "SideChain-H";
    parse();
    AmberMask node = make_shared<AmberMaskAST::Operator>(
        AmberMaskAST::Op::AND,
        make_shared<AmberMaskAST::Operator>(
            AmberMaskAST::Op::AND, make_shared<AmberMaskAST::residue_name_nums>(AMBERMASK::protein),
            make_shared<AmberMaskAST::Operator>(AmberMaskAST::Op::NOT,
                                                make_shared<AmberMaskAST::atom_name_nums>(AmberMaskAST::select_ranges{
                                                    "CA", "C", "O", "N", "H", "OC1", "OC2", "H1", "H2", "H3"}))),
        make_shared<AmberMaskAST::Operator>(
            AmberMaskAST::Op::NOT, make_shared<AmberMaskAST::atom_name_nums>(AmberMaskAST::select_ranges{"H*"}))

    );
    ASSERT_THAT(mask, Eq(node));
}

TEST_F(TestGrammarMacro, DNA) {
    input_string = "DNA";
    parse();

    AmberMask node = make_shared<AmberMaskAST::residue_name_nums>(AMBERMASK::dna);

    ASSERT_THAT(mask, Eq(node));
}

TEST_F(TestGrammarMacro, DNA_H) {
    input_string = "DNA-H";
    parse();

    AmberMask node = make_shared<AmberMaskAST::Operator>(
        AmberMaskAST::Op::AND, make_shared<AmberMaskAST::residue_name_nums>(AMBERMASK::dna),
        make_shared<AmberMaskAST::Operator>(
            AmberMaskAST::Op::NOT, make_shared<AmberMaskAST::atom_name_nums>(AmberMaskAST::select_ranges{"H*"})));

    ASSERT_THAT(mask, Eq(node));
}

TEST_F(TestGrammarMacro, RNA) {
    input_string = "RNA";
    parse();

    AmberMask node = make_shared<AmberMaskAST::residue_name_nums>(AMBERMASK::rna);

    ASSERT_THAT(mask, Eq(node));
}

TEST_F(TestGrammarMacro, RNA_H) {
    input_string = "RNA-H";
    parse();

    AmberMask node = make_shared<AmberMaskAST::Operator>(
        AmberMaskAST::Op::AND, make_shared<AmberMaskAST::residue_name_nums>(AMBERMASK::rna),
        make_shared<AmberMaskAST::Operator>(
            AmberMaskAST::Op::NOT, make_shared<AmberMaskAST::atom_name_nums>(AmberMaskAST::select_ranges{"H*"})));

    ASSERT_THAT(mask, Eq(node));
}

TEST_F(TestGrammarMacro, Water) {
    input_string = "Water";
    parse();

    AmberMask node = make_shared<AmberMaskAST::residue_name_nums>(AMBERMASK::water);

    ASSERT_THAT(mask, Eq(node));
}
