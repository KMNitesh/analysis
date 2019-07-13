//
// Created by xiamr on 7/13/19.
//

#include <gmock/gmock.h>

#include "TwoAtomSelectorGrammar.hpp"

using namespace std;
using namespace testing;


TEST(TwoAtomSelectorGrammarTest, Simple) {

    TwoAtomVectorGrammar<std::string::iterator, qi::ascii::space_type> grammar;

    std::string input_string = "twoAtomVector([:1], [@O])";
    TwoAtomVectorSelectorNode ast;
    auto it = input_string.begin();
    bool status = qi::phrase_parse(it, input_string.end(), grammar, qi::ascii::space, ast);

    ASSERT_TRUE(status && it == input_string.end());

    ASSERT_EQ(ast, make_shared<TwoAtomVectorSelectorNodeStruct>(
            make_shared<Atom::residue_name_nums>(
                    Atom::select_ranges{Atom::numItemType(1, boost::optional<pair<uint, int>>{})}),
            make_shared<Atom::atom_name_nums>(Atom::select_ranges{string{"O"}}))
    );

}