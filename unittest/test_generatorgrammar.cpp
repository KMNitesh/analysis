//
// Created by xiamr on 7/8/19.
//

#include <gmock/gmock.h>

#include "dsl/GeneratorGrammar.hpp"

using namespace std;
using namespace testing;

class GereratorGrammarTest : public Test {

protected:
    std::string generated;
    Atom::select_ranges selections;
};

TEST_F(GereratorGrammarTest, AtomElementNames) {

    Atom::Node node = make_shared<Atom::atom_element_names>(std::vector<string>{"N", "O"});

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("@/N,O"));
}

TEST_F(GereratorGrammarTest, ResidueNames) {


    Atom::Node node = make_shared<Atom::residue_name_nums>(
            Atom::select_ranges{string{"TOL"}, Atom::numItemType(10, boost::optional<pair<uint, int>>{})});

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq(":TOL,10"));
}

TEST_F(GereratorGrammarTest, ResidueNum) {

    Atom::Node node = make_shared<Atom::residue_name_nums>(
            Atom::select_ranges{{Atom::numItemType(10, boost::optional<pair<uint, int>>{})}});

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq(":10"));
}

TEST_F(GereratorGrammarTest, AtomNameNums) {

    Atom::Node node = make_shared<Atom::atom_name_nums>(Atom::select_ranges{
            Atom::numItemType(10, make_pair<uint>(20, 2)),
            string{"OW"},
            string{"HW"}
    });

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("@10-20#2,OW,HW"));
}

TEST_F(GereratorGrammarTest, AtomType) {

    Atom::Node node = make_shared<Atom::atom_types>(Atom::select_ranges{
            Atom::numItemType(10, make_pair<uint>(20, 2)),
            string{"OW"},
            string{"HW"}
    });

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("@%10-20#2,OW,HW"));
}

TEST_F(GereratorGrammarTest, AtomTypeWithStepOne) {

    Atom::Node node = make_shared<Atom::atom_types>(Atom::select_ranges{
            Atom::numItemType(10, make_pair<uint>(20, 1)),
            string{"OW"},
            string{"HW"}
    });


    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("@%10-20,OW,HW"));
}

TEST_F(GereratorGrammarTest, NotOperator) {

    Atom::Node node = make_shared<Atom::Operator>(
            Atom::Op::NOT,
            make_shared<Atom::atom_name_nums>(Atom::select_ranges{string{"OW"}, string{"HW"}}));

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("!@OW,HW"));
}

TEST_F(GereratorGrammarTest, ANDOperator) {

    selections = {string{"OW"}, string{"HW"}};

    Atom::Node node = make_shared<Atom::Operator>(
            Atom::Op::AND,
            make_shared<Atom::residue_name_nums>(selections),
            make_shared<Atom::atom_name_nums>(selections));

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq(":OW,HW&@OW,HW"));
}

TEST_F(GereratorGrammarTest, NOT_AND_Operator) {

    selections = {string{"OW"}, string{"HW"}};

    Atom::Node node = make_shared<Atom::Operator>(
            Atom::Op::NOT, make_shared<Atom::Operator>(
                    Atom::Op::AND,
                    make_shared<Atom::residue_name_nums>(selections),
                    make_shared<Atom::atom_name_nums>(selections)));

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("!(:OW,HW&@OW,HW)"));
}

TEST_F(GereratorGrammarTest, OROperator) {

    selections = {string{"OW"}, string{"HW"}};

    Atom::Node node = make_shared<Atom::Operator>(
            Atom::Op::OR, make_shared<Atom::residue_name_nums>(selections),
            make_shared<Atom::atom_name_nums>(selections));

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq(":OW,HW|@OW,HW"));
}

TEST_F(GereratorGrammarTest, NOT_OR_Operator) {

    selections = {string{"OW"}, string{"HW"}};

    Atom::Node node = make_shared<Atom::Operator>(
            Atom::Op::NOT, make_shared<Atom::Operator>(
                    Atom::Op::OR,
                    make_shared<Atom::residue_name_nums>(selections),
                    make_shared<Atom::atom_name_nums>(selections)));

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("!(:OW,HW|@OW,HW)"));
}

TEST_F(GereratorGrammarTest, AND_OR_Operator) {

    selections = {string{"OW"}, string{"HW"}};

    Atom::Node node = make_shared<Atom::Operator>(
            Atom::Op::AND,
            make_shared<Atom::Operator>(
                    Atom::Op::OR,
                    make_shared<Atom::atom_types>(selections),
                    make_shared<Atom::atom_name_nums>(selections)),
            make_shared<Atom::Operator>(
                    Atom::Op::OR,
                    make_shared<Atom::residue_name_nums>(selections),
                    make_shared<Atom::atom_name_nums>(selections)));

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("(@%OW,HW|@OW,HW)&(:OW,HW|@OW,HW)"));
}

TEST_F(GereratorGrammarTest, OR_AND_Operator) {

    selections = {string{"OW"}, string{"HW"}};

    Atom::Node node = make_shared<Atom::Operator>(
            Atom::Op::OR,
            make_shared<Atom::Operator>(
                    Atom::Op::AND,
                    make_shared<Atom::atom_types>(selections),
                    make_shared<Atom::atom_name_nums>(selections)),
            make_shared<Atom::Operator>(
                    Atom::Op::AND,
                    make_shared<Atom::residue_name_nums>(selections),
                    make_shared<Atom::atom_name_nums>(selections)));

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("@%OW,HW&@OW,HW|:OW,HW&@OW,HW"));
}

TEST_F(GereratorGrammarTest, NOT_OR_AND_Operator) {

    selections = {string{"OW"}, string{"HW"}};

    Atom::Node node =
            make_shared<Atom::Operator>(
                    Atom::Op::NOT,
                    make_shared<Atom::Operator>
                            (Atom::Op::OR,
                             make_shared<Atom::Operator>(
                                     Atom::Op::AND,
                                     make_shared<Atom::atom_types>(
                                             selections),
                                     make_shared<Atom::atom_name_nums>(
                                             selections)),
                             make_shared<Atom::Operator>(
                                     Atom::Op::AND,
                                     make_shared<Atom::residue_name_nums>(
                                             selections),
                                     make_shared<Atom::atom_name_nums>(
                                             selections))));

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("!(@%OW,HW&@OW,HW|:OW,HW&@OW,HW)"));
}

TEST_F(GereratorGrammarTest, OR_OR_Operator) {

    selections = {string{"OW"}, string{"HW"}};

    Atom::Node node = make_shared<Atom::Operator>(
            Atom::Op::OR,
            make_shared<Atom::Operator>(
                    Atom::Op::OR, make_shared<Atom::atom_types>(selections),
                    make_shared<Atom::atom_name_nums>(selections)),
            make_shared<Atom::Operator>(
                    Atom::Op::OR,
                    make_shared<Atom::residue_name_nums>(selections),
                    make_shared<Atom::atom_name_nums>(selections)));

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("@%OW,HW|@OW,HW|:OW,HW|@OW,HW"));
}

TEST_F(GereratorGrammarTest, AND_AND_Operator) {

    selections = {string{"OW"}, string{"HW"}};

    Atom::Node node = make_shared<Atom::Operator>(
            Atom::Op::AND,
            make_shared<Atom::Operator>(
                    Atom::Op::AND, make_shared<Atom::atom_types>(selections),
                    make_shared<Atom::atom_name_nums>(selections)),
            make_shared<Atom::Operator>(
                    Atom::Op::AND,
                    make_shared<Atom::residue_name_nums>(selections),
                    make_shared<Atom::atom_name_nums>(selections)));

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("@%OW,HW&@OW,HW&:OW,HW&@OW,HW"));
}

TEST_F(GereratorGrammarTest, NOT_NOT_Operator) {

    selections = {string{"OW"}, string{"HW"}};

    Atom::Node node = make_shared<Atom::Operator>(
            Atom::Op::NOT,
            make_shared<Atom::Operator>(
                    Atom::Op::NOT, make_shared<Atom::atom_types>(selections)));

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("!!@%OW,HW"));
}

