//
// Created by xiamr on 7/13/19.
//


#include <gmock/gmock.h>

#include "LanguageGrammar.hpp"

using namespace std;
using namespace testing;

TEST(LanguageGrammarTest, Simple) {

    LanguageGrammar<std::string::iterator, qi::ascii::space_type> grammar;

    std::string input_string = "rotacf(vector = normalVector([:1], [@O], [@1500]), out = x-x-x.x, P=1,"
                               " time_increment_ps = 0.1, max_time_grap_ps = 100);"
                               "go(start = 1, end = 1000, step = 1, nthreads =12)";
    Language ast;
    auto it = input_string.begin();
    bool status = qi::phrase_parse(it, input_string.end(), grammar, qi::ascii::space, ast);

    ASSERT_TRUE(status && it == input_string.end());

    auto rotRafNode = make_shared<RotAcfNodeStruct>();
    rotRafNode->vectorSelctor = vectorSelectorNode{make_shared<NormalVectorSelectorNodeStruct>(
            make_shared<Atom::residue_name_nums>(
                    Atom::select_ranges{Atom::numItemType(1, boost::optional<pair<uint, int>>{})}),
            make_shared<Atom::atom_name_nums>(Atom::select_ranges{string{"O"}}),
            make_shared<Atom::atom_name_nums>(
                    Atom::select_ranges{Atom::numItemType(1500, boost::optional<pair<uint, int>>{})}))};
    rotRafNode->Legendre = 1;
    rotRafNode->time_increment_ps = 0.1;
    rotRafNode->max_time_grap_ps = 100;
    rotRafNode->outfilename = "x-x-x.x";

    auto goNode = make_shared<GoNodeStruct>();
    goNode->start = 1;
    goNode->end = 1000;
    goNode->step = 1;
    goNode->nthreads = 12;
    vector<shared_ptr<LanguageNodeVariant>> node = {make_shared<LanguageNodeVariant>(rotRafNode),
                                                    make_shared<LanguageNodeVariant>(goNode)};

    ASSERT_EQ(ast, node);
}
