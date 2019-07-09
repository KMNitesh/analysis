//
// Created by xiamr on 7/8/19.
//

#include <gmock/gmock.h>

#include "GeneratorGrammar.hpp"

using namespace std;
using namespace testing;

class GereratorGrammarTest : public Test {

protected:
    std::string generated;
    std::vector<boost::variant<boost::fusion::vector<uint, boost::optional<std::pair<uint, int>>>, std::string>> selections;
};

TEST_F(GereratorGrammarTest, AtomElementNames) {

    Atom::Node node = make_shared<Atom::atom_element_names>(std::vector<string>{"N", "O"});

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("@/N,O"));
}

TEST_F(GereratorGrammarTest, ResidueNames) {

    selections.push_back(string{"TOL"});
    selections.push_back(
            fusion::vector<uint, boost::optional<std::pair<uint, int>>>(10, boost::optional<std::pair<uint, int>>{}));

    Atom::Node node = make_shared<Atom::residue_name_nums>(selections);

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq(":TOL,10"));
}

TEST_F(GereratorGrammarTest, ResidueNum) {

    selections.push_back(
            fusion::vector<uint, boost::optional<std::pair<uint, int>>>(10, boost::optional<std::pair<uint, int>>{}));

    Atom::Node node = make_shared<Atom::residue_name_nums>(selections);

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq(":10"));
}

TEST_F(GereratorGrammarTest, AtomNameNums) {

    selections.push_back(fusion::vector<uint, boost::optional<std::pair<uint, int>>>(10, make_pair<uint>(20, 2)));
    selections.push_back(string{"OW"});
    selections.push_back(string{"HW"});

    Atom::Node node = make_shared<Atom::atom_name_nums>(selections);

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("@10-20#2,OW,HW"));
}

TEST_F(GereratorGrammarTest, AtomType) {

    selections.push_back(fusion::vector<uint, boost::optional<std::pair<uint, int>>>(10, make_pair<uint>(20, 2)));
    selections.push_back(string{"OW"});
    selections.push_back(string{"HW"});

    Atom::Node node = make_shared<Atom::atom_types>(selections);

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("@%10-20#2,OW,HW"));
}

TEST_F(GereratorGrammarTest, AtomTypeWithStepOne) {

    selections.push_back(fusion::vector<uint, boost::optional<std::pair<uint, int>>>(10, make_pair<uint>(20, 1)));
    selections.push_back(string{"OW"});
    selections.push_back(string{"HW"});

    Atom::Node node = make_shared<Atom::atom_types>(selections);

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("@%10-20,OW,HW"));
}

TEST_F(GereratorGrammarTest, NotOperator) {

    selections.push_back(string{"OW"});
    selections.push_back(string{"HW"});

    Atom::Node node = make_shared<Atom::Operator>(
            Atom::Op::NOT,
            make_shared<Atom::atom_name_nums>(selections));

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("!@OW,HW"));
}

TEST_F(GereratorGrammarTest, ANDOperator) {

    selections.push_back(string{"OW"});
    selections.push_back(string{"HW"});

    Atom::Node node = make_shared<Atom::Operator>(
            Atom::Op::AND,
            make_shared<Atom::residue_name_nums>(selections),
            make_shared<Atom::atom_name_nums>(selections));

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq(":OW,HW&@OW,HW"));
}

TEST_F(GereratorGrammarTest, NOT_AND_Operator) {

    selections.push_back(string{"OW"});
    selections.push_back(string{"HW"});

    Atom::Node node = make_shared<Atom::Operator>(
            Atom::Op::NOT, make_shared<Atom::Operator>(
                    Atom::Op::AND,
                    make_shared<Atom::residue_name_nums>(selections),
                    make_shared<Atom::atom_name_nums>(selections)));

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("!(:OW,HW&@OW,HW)"));
}

TEST_F(GereratorGrammarTest, OROperator) {

    selections.push_back(string{"OW"});
    selections.push_back(string{"HW"});

    Atom::Node node = make_shared<Atom::Operator>(
            Atom::Op::OR, make_shared<Atom::residue_name_nums>(selections),
            make_shared<Atom::atom_name_nums>(selections));

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq(":OW,HW|@OW,HW"));
}

TEST_F(GereratorGrammarTest, NOT_OR_Operator) {

    selections.push_back(string{"OW"});
    selections.push_back(string{"HW"});

    Atom::Node node = make_shared<Atom::Operator>(
            Atom::Op::NOT, make_shared<Atom::Operator>(
                    Atom::Op::OR,
                    make_shared<Atom::residue_name_nums>(selections),
                    make_shared<Atom::atom_name_nums>(selections)));

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("!(:OW,HW|@OW,HW)"));
}

TEST_F(GereratorGrammarTest, AND_OR_Operator) {

    selections.push_back(string{"OW"});
    selections.push_back(string{"HW"});

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

    selections.push_back(string{"OW"});
    selections.push_back(string{"HW"});

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

    selections.push_back(string{"OW"});
    selections.push_back(string{"HW"});

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

    selections.push_back(string{"OW"});
    selections.push_back(string{"HW"});

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

    selections.push_back(string{"OW"});
    selections.push_back(string{"HW"});

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

    selections.push_back(string{"OW"});
    selections.push_back(string{"HW"});

    Atom::Node node = make_shared<Atom::Operator>(
            Atom::Op::NOT,
            make_shared<Atom::Operator>(
                    Atom::Op::NOT, make_shared<Atom::atom_types>(selections)));

    ASSERT_TRUE(format_node(node, generated));

    ASSERT_THAT(generated, StrEq("!!@%OW,HW"));
}

