//
// Created by xiamr on 7/13/19.
//


#include <gmock/gmock.h>

#include "NormalVectorSelectorGrammar.hpp"

using namespace std;
using namespace testing;


TEST(NormalVectorSelectorGrammarTest, Simple) {

    NormalVectorGrammar<std::string::iterator, qi::ascii::space_type> grammar;

    std::string input_string = "normalVector([:1], [@O], [@1500])";
    NormalVectorSelectorNode ast;
    auto it = input_string.begin();
    bool status = qi::phrase_parse(it, input_string.end(), grammar, qi::ascii::space, ast);

    ASSERT_TRUE(status && it == input_string.end());

    ASSERT_EQ(ast, make_shared<NormalVectorSelectorNodeStruct>(
            make_shared<Atom::residue_name_nums>(
                    Atom::select_ranges{Atom::numItemType(1, boost::optional<pair<uint, int>>{})}),
            make_shared<Atom::atom_name_nums>(Atom::select_ranges{string{"O"}}),
            make_shared<Atom::atom_name_nums>(
                    Atom::select_ranges{Atom::numItemType(1500, boost::optional<pair<uint, int>>{})})
    ));

}