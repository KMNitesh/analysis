//
// Created by xiamr on 6/2/19.
//

#define BOOST_TEST_MODULE "Gromacs Style Center Selection Grammer Test"

#include <string>
#include <memory>

#include <boost/test/included/unit_test.hpp>

#include "center_selection_grammar.hpp"

BOOST_AUTO_TEST_SUITE(test_center_selection_grammer)

BOOST_AUTO_TEST_CASE(test_com) {
    CenterGrammar<std::string::iterator,qi::ascii::space_type> grammar;

    std::string input_string = "com of :1 ";
    CenterRuleNode  ast;
    auto it = input_string.begin();
    bool status = qi::phrase_parse(it, input_string.end(), grammar, qi::ascii::space, ast);

    auto ret = boost::get<std::shared_ptr<MassCenterRuleNode>>(&ast);

    BOOST_CHECK_MESSAGE(status && it == input_string.end() && ret, "test com");
    if (ret){
         boost::apply_visitor(print(), (*ret)->SelectionMask);
    }
}

BOOST_AUTO_TEST_CASE(test_not_com) {
    CenterGrammar<std::string::iterator,qi::ascii::space_type> grammar;

    std::string input_string = "comdd of :1 ";
    CenterRuleNode  ast;
    auto it = input_string.begin();
    bool status = qi::phrase_parse(it, input_string.end(), grammar, qi::ascii::space, ast);

    auto ret = boost::get<std::shared_ptr<MassCenterRuleNode>>(&ast);

    BOOST_CHECK_MESSAGE(!(status && it == input_string.end() && ret), "test not com");
}

BOOST_AUTO_TEST_CASE(test_geom) {
    CenterGrammar<std::string::iterator,qi::ascii::space_type> grammar;

    std::string input_string = "geom of :1& @O ";
    CenterRuleNode  ast;
    auto it = input_string.begin();
    bool status = qi::phrase_parse(it, input_string.end(), grammar, qi::ascii::space, ast);

    auto ret = boost::get<std::shared_ptr<GeomCenterRuleNode>>(&ast);

    BOOST_CHECK_MESSAGE(status && it == input_string.end() && ret, "test geom");
    if (ret){
         boost::apply_visitor(print(), (*ret)->SelectionMask);
    }
}

BOOST_AUTO_TEST_CASE(test_not_geom) {
    CenterGrammar<std::string::iterator,qi::ascii::space_type> grammar;

    std::string input_string = "com of :1& @O ";
    CenterRuleNode  ast;
    auto it = input_string.begin();
    bool status = qi::phrase_parse(it, input_string.end(), grammar, qi::ascii::space, ast);

    auto ret = boost::get<std::shared_ptr<GeomCenterRuleNode>>(&ast);

    BOOST_CHECK_MESSAGE(!(status && it == input_string.end() && ret), "test not geom");
    if (ret){
         boost::apply_visitor(print(), (*ret)->SelectionMask);
    }
}


BOOST_AUTO_TEST_CASE(test_noop) {
    CenterGrammar<std::string::iterator,qi::ascii::space_type> grammar;

    std::string input_string = ":1& @O ";
    CenterRuleNode  ast;
    auto it = input_string.begin();
    bool status = qi::phrase_parse(it, input_string.end(), grammar, qi::ascii::space, ast);

    auto ret = boost::get<std::shared_ptr<NoopRuleNode>>(&ast);

    BOOST_CHECK_MESSAGE(status && it == input_string.end() && ret, "test noop");
    if (ret){
         boost::apply_visitor(print(), (*ret)->SelectionMask);
    }
}


BOOST_AUTO_TEST_SUITE_END()