
#include <gmock/gmock.h>
#include "dsl/center_selection_grammar.hpp"

using namespace testing;

TEST(test_center_selection_grammer, test_com) {
    CenterGrammar<std::string::iterator, qi::ascii::space_type> grammar;

    std::string input_string = "com of :1 ";
    CenterRuleNode ast;
    auto it = input_string.begin();
    bool status = qi::phrase_parse(it, input_string.end(), grammar, qi::ascii::space, ast);

    auto ret = boost::get<std::shared_ptr<MassCenterRuleNode>>(&ast);

    ASSERT_TRUE(status && it == input_string.end() && ret);
}

TEST(test_center_selection_grammer, test_not_com) {
    CenterGrammar<std::string::iterator, qi::ascii::space_type> grammar;

    std::string input_string = "comdd of :1 ";
    CenterRuleNode ast;
    auto it = input_string.begin();
    bool status = qi::phrase_parse(it, input_string.end(), grammar, qi::ascii::space, ast);

    auto ret = boost::get<std::shared_ptr<MassCenterRuleNode>>(&ast);

    ASSERT_TRUE(!(status && it == input_string.end() && ret));
}

TEST(test_center_selection_grammer, test_geom) {
    CenterGrammar<std::string::iterator, qi::ascii::space_type> grammar;

    std::string input_string = "geom of :1& @O ";
    CenterRuleNode ast;
    auto it = input_string.begin();
    bool status = qi::phrase_parse(it, input_string.end(), grammar, qi::ascii::space, ast);

    auto ret = boost::get<std::shared_ptr<GeomCenterRuleNode>>(&ast);

    ASSERT_TRUE(status && it == input_string.end() && ret);
}

TEST(test_center_selection_grammer, test_not_geom) {
    CenterGrammar<std::string::iterator, qi::ascii::space_type> grammar;

    std::string input_string = "com of :1& @O ";
    CenterRuleNode ast;
    auto it = input_string.begin();
    bool status = qi::phrase_parse(it, input_string.end(), grammar, qi::ascii::space, ast);

    auto ret = boost::get<std::shared_ptr<GeomCenterRuleNode>>(&ast);

    ASSERT_TRUE(!(status && it == input_string.end() && ret));
}

TEST(test_center_selection_grammer, test_noop) {
    CenterGrammar<std::string::iterator, qi::ascii::space_type> grammar;

    std::string input_string = ":1& @O ";
    CenterRuleNode ast;
    auto it = input_string.begin();
    bool status = qi::phrase_parse(it, input_string.end(), grammar, qi::ascii::space, ast);

    auto ret = boost::get<std::shared_ptr<NoopRuleNode>>(&ast);

    ASSERT_TRUE(status && it == input_string.end() && ret);
}
