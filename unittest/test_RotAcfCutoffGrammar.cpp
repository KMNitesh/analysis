//
// Created by xiamr on 7/13/19.
//

#include <gmock/gmock.h>

#include "RotAcfCutoffGrammar.hpp"
#include "atom.hpp"

using namespace std;
using namespace testing;

TEST(RotAcfCutoffGrammarTest, Simple) {

    RotAcfCutoffGrammar<std::string::iterator, qi::ascii::space_type> grammar;

    std::string input_string = "rotacfcutoff(M = [@U], L= [@OW], vector = "
                               "normalVector([:1], [@O], [@1500]),"
                               " P=1, cutoff = 3.0, time_increment_ps = 0.1, out = my-raf.xvg  , max_time_grap_ps = 100)";
    RotAcfCutoffNode ast;
    auto it = input_string.begin();
    bool status = qi::phrase_parse(it, input_string.end(), grammar, qi::ascii::space, ast);

    ASSERT_TRUE(status && it == input_string.end());
    auto node = make_shared<RotAcfCutoffNodeStruct>();
    node->vectorSelctor = vectorSelectorNode{make_shared<NormalVectorSelectorNodeStruct>(
            make_shared<Atom::residue_name_nums>(
                    Atom::select_ranges{Atom::numItemType(1, boost::optional<pair<uint, int>>{})}),
            make_shared<Atom::atom_name_nums>(Atom::select_ranges{string{"O"}}),
            make_shared<Atom::atom_name_nums>(
                    Atom::select_ranges{Atom::numItemType(1500, boost::optional<pair<uint, int>>{})}))};
    node->M = make_shared<Atom::atom_name_nums>(Atom::select_ranges{string{"U"}});
    node->L = make_shared<Atom::atom_name_nums>(Atom::select_ranges{string{"OW"}});
    node->Legendre = 1;
    node->cutoff = 3.0;
    node->time_increment_ps = 0.1;
    node->max_time_grap_ps = 100;
    node->outfilename = "my-raf.xvg";

    ASSERT_EQ(ast, node);
}