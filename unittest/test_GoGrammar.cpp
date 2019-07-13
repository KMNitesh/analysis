//
// Created by xiamr on 7/13/19.
//

#include <gmock/gmock.h>

#include "GoGrammar.hpp"

using namespace std;
using namespace testing;

TEST(GoGrammarTest, Simple) {

    GoGrammar<std::string::iterator, qi::ascii::space_type> grammar;

    std::string input_string = "go(start = 1, end = 1000, step = 1, nthreads =12)";
    GoNode ast;
    auto it = input_string.begin();
    bool status = qi::phrase_parse(it, input_string.end(), grammar, qi::ascii::space, ast);

    ASSERT_TRUE(status && it == input_string.end());

    auto goNode = make_shared<GoNodeStruct>();

    goNode->start = 1;
    goNode->end = 1000;
    goNode->step = 1;
    goNode->nthreads = 12;

    ASSERT_EQ(ast, goNode);
}


