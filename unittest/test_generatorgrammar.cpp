//
// Created by xiamr on 7/8/19.
//

#include <gmock/gmock.h>

#include "dsl/GeneratorGrammar.hpp"

using namespace std;
using namespace testing;

class GeneratorGrammarTest : public Test {

protected:
    std::string generated;
    AmberMaskAST::select_ranges selections;
};

TEST_F(GeneratorGrammarTest, AtomElementNames) {

    AmberMask node = make_shared<AmberMaskAST::atom_element_names>(std::vector<AmberMaskAST::Name>{"N", "O"});

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("@/N,O"));
}

TEST_F(GeneratorGrammarTest, ResidueNames) {

    AmberMask node = make_shared<AmberMaskAST::residue_name_nums>(
        AmberMaskAST::select_ranges{string{"TOL"}, AmberMaskAST::numItemType(10, boost::optional<pair<uint, int>>{})});

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq(":TOL,10"));
}

TEST_F(GeneratorGrammarTest, ResidueNum) {

    AmberMask node = make_shared<AmberMaskAST::residue_name_nums>(
        AmberMaskAST::select_ranges{{AmberMaskAST::numItemType(10, boost::optional<pair<uint, int>>{})}});

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq(":10"));
}

TEST_F(GeneratorGrammarTest, MoleculeNum) {

    AmberMask node = make_shared<AmberMaskAST::molecule_nums>(std::vector<AmberMaskAST::molecule_nums::Attr>{
        {AmberMaskAST::molecule_nums::Attr(10, boost::optional<boost::fusion::vector<uint, boost::optional<int>>>{})}});

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("$10"));
}

TEST_F(GeneratorGrammarTest, MoleculeNumRangeWithoutStep) {

    AmberMask node =
        make_shared<AmberMaskAST::molecule_nums>(std::vector<AmberMaskAST::molecule_nums::Attr>{{AmberMaskAST::molecule_nums::Attr(
            10, boost::optional<boost::fusion::vector<uint, boost::optional<int>>>{
                    boost::fusion::vector<uint, boost::optional<int>>{20, boost::optional<int>{}}})}});

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("$10-20"));
}

TEST_F(GeneratorGrammarTest, MoleculeNumRangeWithStep) {

    AmberMask node = make_shared<AmberMaskAST::molecule_nums>(std::vector<AmberMaskAST::molecule_nums::Attr>{
        {AmberMaskAST::molecule_nums::Attr(10, boost::optional<boost::fusion::vector<uint, boost::optional<int>>>{
                                           boost::fusion::vector<uint, boost::optional<int>>{20, 2}})}});

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("$10-20#2"));
}

TEST_F(GeneratorGrammarTest, MoleculeNumRangeWithStepOne) {

    AmberMask node = make_shared<AmberMaskAST::molecule_nums>(std::vector<AmberMaskAST::molecule_nums::Attr>{
        {AmberMaskAST::molecule_nums::Attr(10, boost::optional<boost::fusion::vector<uint, boost::optional<int>>>{
                                           boost::fusion::vector<uint, boost::optional<int>>{20, 1}})}});

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("$10-20"));
}

TEST_F(GeneratorGrammarTest, AtomNameNums) {

    AmberMask node = make_shared<AmberMaskAST::atom_name_nums>(
        AmberMaskAST::select_ranges{AmberMaskAST::numItemType(10, make_pair<uint>(20, 2)), string{"OW"}, string{"HW"}});

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("@10-20#2,OW,HW"));
}

TEST_F(GeneratorGrammarTest, AtomType) {

    AmberMask node = make_shared<AmberMaskAST::atom_types>(
        AmberMaskAST::select_ranges{AmberMaskAST::numItemType(10, make_pair<uint>(20, 2)), string{"OW"}, string{"HW"}});

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("@%10-20#2,OW,HW"));
}

TEST_F(GeneratorGrammarTest, AtomTypeWithStepOne) {

    AmberMask node = make_shared<AmberMaskAST::atom_types>(
        AmberMaskAST::select_ranges{AmberMaskAST::numItemType(10, make_pair<uint>(20, 1)), string{"OW"}, string{"HW"}});

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("@%10-20,OW,HW"));
}

TEST_F(GeneratorGrammarTest, NotOperator) {

    AmberMask node = make_shared<AmberMaskAST::Operator>(
        AmberMaskAST::Op::NOT, make_shared<AmberMaskAST::atom_name_nums>(AmberMaskAST::select_ranges{AmberMaskAST::Name{"OW"}, AmberMaskAST::Name{"HW"}}));

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("!@OW,HW"));
}

TEST_F(GeneratorGrammarTest, ANDOperator) {

    selections = {AmberMaskAST::Name{"OW"}, AmberMaskAST::Name{"HW"}};

    AmberMask node = make_shared<AmberMaskAST::Operator>(AmberMaskAST::Op::AND, make_shared<AmberMaskAST::residue_name_nums>(selections),
                                                  make_shared<AmberMaskAST::atom_name_nums>(selections));

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq(":OW,HW&@OW,HW"));
}

TEST_F(GeneratorGrammarTest, NOT_AND_Operator) {

    selections = {AmberMaskAST::Name{"OW"}, AmberMaskAST::Name{"HW"}};

    AmberMask node = make_shared<AmberMaskAST::Operator>(
        AmberMaskAST::Op::NOT, make_shared<AmberMaskAST::Operator>(AmberMaskAST::Op::AND, make_shared<AmberMaskAST::residue_name_nums>(selections),
                                                   make_shared<AmberMaskAST::atom_name_nums>(selections)));

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("!(:OW,HW&@OW,HW)"));
}

TEST_F(GeneratorGrammarTest, OROperator) {

    selections = {string{"OW"}, string{"HW"}};

    AmberMask node = make_shared<AmberMaskAST::Operator>(AmberMaskAST::Op::OR, make_shared<AmberMaskAST::residue_name_nums>(selections),
                                                  make_shared<AmberMaskAST::atom_name_nums>(selections));

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq(":OW,HW|@OW,HW"));
}

TEST_F(GeneratorGrammarTest, NOT_OR_Operator) {

    selections = {AmberMaskAST::Name{"OW"}, AmberMaskAST::Name{"HW"}};

    AmberMask node = make_shared<AmberMaskAST::Operator>(
        AmberMaskAST::Op::NOT, make_shared<AmberMaskAST::Operator>(AmberMaskAST::Op::OR, make_shared<AmberMaskAST::residue_name_nums>(selections),
                                                   make_shared<AmberMaskAST::atom_name_nums>(selections)));

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("!(:OW,HW|@OW,HW)"));
}

TEST_F(GeneratorGrammarTest, AND_OR_Operator) {

    selections = {AmberMaskAST::Name{"OW"}, AmberMaskAST::Name{"HW"}};

    AmberMask node = make_shared<AmberMaskAST::Operator>(
        AmberMaskAST::Op::AND,
        make_shared<AmberMaskAST::Operator>(AmberMaskAST::Op::OR, make_shared<AmberMaskAST::atom_types>(selections),
                                    make_shared<AmberMaskAST::atom_name_nums>(selections)),
        make_shared<AmberMaskAST::Operator>(AmberMaskAST::Op::OR, make_shared<AmberMaskAST::residue_name_nums>(selections),
                                    make_shared<AmberMaskAST::atom_name_nums>(selections)));

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("(@%OW,HW|@OW,HW)&(:OW,HW|@OW,HW)"));
}

TEST_F(GeneratorGrammarTest, OR_AND_Operator) {

    selections = {AmberMaskAST::Name{"OW"}, AmberMaskAST::Name{"HW"}};

    AmberMask node = make_shared<AmberMaskAST::Operator>(
        AmberMaskAST::Op::OR,
        make_shared<AmberMaskAST::Operator>(AmberMaskAST::Op::AND, make_shared<AmberMaskAST::atom_types>(selections),
                                    make_shared<AmberMaskAST::atom_name_nums>(selections)),
        make_shared<AmberMaskAST::Operator>(AmberMaskAST::Op::AND, make_shared<AmberMaskAST::residue_name_nums>(selections),
                                    make_shared<AmberMaskAST::atom_name_nums>(selections)));

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("@%OW,HW&@OW,HW|:OW,HW&@OW,HW"));
}

TEST_F(GeneratorGrammarTest, NOT_OR_AND_Operator) {

    selections = {AmberMaskAST::Name{"OW"}, AmberMaskAST::Name{"HW"}};

    AmberMask node = make_shared<AmberMaskAST::Operator>(
        AmberMaskAST::Op::NOT, make_shared<AmberMaskAST::Operator>(
                           AmberMaskAST::Op::OR,
                           make_shared<AmberMaskAST::Operator>(AmberMaskAST::Op::AND, make_shared<AmberMaskAST::atom_types>(selections),
                                                       make_shared<AmberMaskAST::atom_name_nums>(selections)),
                           make_shared<AmberMaskAST::Operator>(AmberMaskAST::Op::AND, make_shared<AmberMaskAST::residue_name_nums>(selections),
                                                       make_shared<AmberMaskAST::atom_name_nums>(selections))));

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("!(@%OW,HW&@OW,HW|:OW,HW&@OW,HW)"));
}

TEST_F(GeneratorGrammarTest, OR_OR_Operator) {

    selections = {AmberMaskAST::Name{"OW"}, AmberMaskAST::Name{"HW"}};

    AmberMask node = make_shared<AmberMaskAST::Operator>(
        AmberMaskAST::Op::OR,
        make_shared<AmberMaskAST::Operator>(AmberMaskAST::Op::OR, make_shared<AmberMaskAST::atom_types>(selections),
                                    make_shared<AmberMaskAST::atom_name_nums>(selections)),
        make_shared<AmberMaskAST::Operator>(AmberMaskAST::Op::OR, make_shared<AmberMaskAST::residue_name_nums>(selections),
                                    make_shared<AmberMaskAST::atom_name_nums>(selections)));

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("@%OW,HW|@OW,HW|:OW,HW|@OW,HW"));
}

TEST_F(GeneratorGrammarTest, AND_AND_Operator) {

    selections = {AmberMaskAST::Name{"OW"}, AmberMaskAST::Name{"HW"}};

    AmberMask node = make_shared<AmberMaskAST::Operator>(
        AmberMaskAST::Op::AND,
        make_shared<AmberMaskAST::Operator>(AmberMaskAST::Op::AND, make_shared<AmberMaskAST::atom_types>(selections),
                                    make_shared<AmberMaskAST::atom_name_nums>(selections)),
        make_shared<AmberMaskAST::Operator>(AmberMaskAST::Op::AND, make_shared<AmberMaskAST::residue_name_nums>(selections),
                                    make_shared<AmberMaskAST::atom_name_nums>(selections)));

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("@%OW,HW&@OW,HW&:OW,HW&@OW,HW"));
}

TEST_F(GeneratorGrammarTest, NOT_NOT_Operator) {

    selections = {AmberMaskAST::Name{"OW"}, AmberMaskAST::Name{"HW"}};

    AmberMask node = make_shared<AmberMaskAST::Operator>(
        AmberMaskAST::Op::NOT, make_shared<AmberMaskAST::Operator>(AmberMaskAST::Op::NOT, make_shared<AmberMaskAST::atom_types>(selections)));

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("!!@%OW,HW"));
}

TEST_F(GeneratorGrammarTest, isSystem) {
    ASSERT_TRUE(format_node(AMBERMASK::System, generated));
    ASSERT_THAT(generated, StrEq("System"));
}

TEST_F(GeneratorGrammarTest, isProtein) {
    ASSERT_TRUE(format_node(AMBERMASK::Protein, generated));
    ASSERT_THAT(generated, StrEq("Protein"));
}

TEST_F(GeneratorGrammarTest, isProtein_H) {
    ASSERT_TRUE(format_node(AMBERMASK::Protein_H, generated));
    ASSERT_THAT(generated, StrEq("Protein-H"));
}

TEST_F(GeneratorGrammarTest, isBackbone) {
    ASSERT_TRUE(format_node(AMBERMASK::Backbone, generated));
    ASSERT_THAT(generated, StrEq("Backbone"));
}

TEST_F(GeneratorGrammarTest, isMainChain) {
    ASSERT_TRUE(format_node(AMBERMASK::MainChain, generated));
    ASSERT_THAT(generated, StrEq("MainChain"));
}

TEST_F(GeneratorGrammarTest, isMainChain_plus_Cb) {
    ASSERT_TRUE(format_node(AMBERMASK::MainChain_plus_Cb, generated));
    ASSERT_THAT(generated, StrEq("MainChain+Cb"));
}

TEST_F(GeneratorGrammarTest, isMainChain_plus_H) {
    ASSERT_TRUE(format_node(AMBERMASK::MainChain_plus_H, generated));
    ASSERT_THAT(generated, StrEq("MainChain+H"));
}

TEST_F(GeneratorGrammarTest, isC_alpha) {
    ASSERT_TRUE(format_node(AMBERMASK::C_alpha, generated));
    ASSERT_THAT(generated, StrEq("C-alpha"));
}

TEST_F(GeneratorGrammarTest, isSideChain) {
    ASSERT_TRUE(format_node(AMBERMASK::SideChain, generated));
    ASSERT_THAT(generated, StrEq("SideChain"));
}

TEST_F(GeneratorGrammarTest, isSideChain_H) {
    ASSERT_TRUE(format_node(AMBERMASK::SideChain_H, generated));
    ASSERT_THAT(generated, StrEq("SideChain-H"));
}

TEST_F(GeneratorGrammarTest, isDNA) {
    ASSERT_TRUE(format_node(AMBERMASK::DNA, generated));
    ASSERT_THAT(generated, StrEq("DNA"));
}

TEST_F(GeneratorGrammarTest, isDNA_H) {
    ASSERT_TRUE(format_node(AMBERMASK::DNA_H, generated));
    ASSERT_THAT(generated, StrEq("DNA-H"));
}

TEST_F(GeneratorGrammarTest, isRNA) {
    ASSERT_TRUE(format_node(AMBERMASK::RNA, generated));
    ASSERT_THAT(generated, StrEq("RNA"));
}

TEST_F(GeneratorGrammarTest, isRNA_H) {
    ASSERT_TRUE(format_node(AMBERMASK::RNA_H, generated));
    ASSERT_THAT(generated, StrEq("RNA-H"));
}

TEST_F(GeneratorGrammarTest, isWater) {
    ASSERT_TRUE(format_node(AMBERMASK::Water, generated));
    ASSERT_THAT(generated, StrEq("Water"));
}
