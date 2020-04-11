
#include "dsl/selection_grammar.hpp"
#include <gmock/gmock.h>

using namespace testing;

TEST(test_center_selection_grammer, test_com) {
    CenterGrammar<std::string::iterator, qi::ascii::space_type> grammar;

    std::string input_string = "com :1 ";
    SelectionAST ast;
    auto it = input_string.begin();
    bool status = qi::phrase_parse(it, input_string.end(), grammar > qi::eoi, qi::ascii::space, ast);

    auto ret = boost::get<MassCenterRuleNode>(&ast);

    ASSERT_TRUE(status && it == input_string.end() && ret);
}

TEST(test_center_selection_grammer, test_not_com) {
    CenterGrammar<std::string::iterator, qi::ascii::space_type> grammar;

    std::string input_string = "comdd :1 ";
    SelectionAST ast;
    auto it = input_string.begin();
    ASSERT_ANY_THROW(qi::phrase_parse(it, input_string.end(), grammar > qi::eoi, qi::ascii::space, ast));
}

TEST(test_center_selection_grammer, test_geom) {
    CenterGrammar<std::string::iterator, qi::ascii::space_type> grammar;

    std::string input_string = "geom  :1&@O ";
    SelectionAST ast;
    auto it = input_string.begin();
    bool status = qi::phrase_parse(it, input_string.end(), grammar > qi::eoi, qi::ascii::space, ast);

    auto ret = boost::get<GeomCenterRuleNode>(&ast);

    ASSERT_TRUE(status && it == input_string.end() && ret);
}

TEST(test_center_selection_grammer, test_not_geom) {
    CenterGrammar<std::string::iterator, qi::ascii::space_type> grammar;

    std::string input_string = "com [:1&@O] ";
    SelectionAST ast;
    auto it = input_string.begin();
    try {
        qi::phrase_parse(it, input_string.end(), grammar > qi::eoi, qi::ascii::space, ast);
    } catch (qi::expectation_failure<std::string::iterator> &x) {
        std::cerr << "Grammar Parse Failure ! Expecting : " << x.what_ << '\n';
        auto column = boost::spirit::get_column(std::begin(input_string), x.first);
        std::string pos = " (column: " + std::to_string(column) + ")";
        std::cerr << pos << ">>>>" << input_string << "<<<<\n";
        std::cerr << std::string(column + pos.size() + 3, ' ') << "^~~~ here\n";
    }
}

TEST(test_center_selection_grammer, test_noop) {
    CenterGrammar<std::string::iterator, qi::ascii::space_type> grammar;

    std::string input_string = ":1& @O ";
    SelectionAST ast;
    auto it = input_string.begin();
    bool status = qi::phrase_parse(it, input_string.end(), grammar > qi::eoi, qi::ascii::space, ast);

    auto ret = boost::get<NoopRuleNode>(&ast);

    ASSERT_TRUE(status && it == input_string.end() && ret);
}
