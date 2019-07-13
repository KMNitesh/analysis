//
// Created by xiamr on 7/13/19.
//

#include <gmock/gmock.h>

#include "RotACFGrammar.hpp"

using namespace std;
using namespace testing;

TEST(RotACFGrammarTest, Simple) {

    RotAcfGrammar<std::string::iterator, qi::ascii::space_type> grammar;

    std::string input_string = "rotacf(vector = normalVector([:1], [@O], [@1500]), "
                               "P=1, time_increment_ps = 0.1, max_time_grap_ps = 100,"
                               "out = a.xvg)";
    RotAcfNode ast;
    auto it = input_string.begin();
    bool status = qi::phrase_parse(it, input_string.end(), grammar, qi::ascii::space, ast);

    ASSERT_TRUE(status && it == input_string.end());
    auto node = make_shared<RotAcfNodeStruct>();
    node->vectorSelctor = vectorSelectorNode{make_shared<NormalVectorSelectorNodeStruct>(
            make_shared<Atom::residue_name_nums>(
                    Atom::select_ranges{Atom::numItemType(1, boost::optional<pair<uint, int>>{})}),
            make_shared<Atom::atom_name_nums>(Atom::select_ranges{string{"O"}}),
            make_shared<Atom::atom_name_nums>(
                    Atom::select_ranges{Atom::numItemType(1500, boost::optional<pair<uint, int>>{})}))};
    node->Legendre = 1;
    node->time_increment_ps = 0.1;
    node->max_time_grap_ps = 100;
    node->outfilename = "a.xvg";

    ASSERT_EQ(ast, node);
}