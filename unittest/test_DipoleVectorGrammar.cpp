//
// Created by xiamr on 7/13/19.
//

#include <gmock/gmock.h>

#include "DipoleVectorSelectorGrammar.hpp"

using namespace std;
using namespace testing;

TEST(DipoleVectorGrammarTest, Simple) {

    DipoleVectorGrammar<std::string::iterator, qi::ascii::space_type> grammar;

    std::string input_string = "dipoleVector([:1])";
    DipoleVectorSelectorNode ast;
    auto it = input_string.begin();
    bool status = qi::phrase_parse(it, input_string.end(), grammar, qi::ascii::space, ast);

    ASSERT_TRUE(status && it == input_string.end());

    ASSERT_EQ(ast, make_shared<DipoleVectorSelectorNodeStruct>(
            make_shared<Atom::residue_name_nums>(
                    Atom::select_ranges{Atom::numItemType(1, boost::optional<pair<uint, int>>{})}))
    );
}

