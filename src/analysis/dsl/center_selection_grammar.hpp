//
// Created by xiamr on 6/1/19.
//

#ifndef TINKER_CENTER_SELECTION_GRAMMAR_HPP
#define TINKER_CENTER_SELECTION_GRAMMAR_HPP

#include <boost/algorithm/string.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/phoenix.hpp>
#include <boost/regex.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/repository/include/qi_distinct.hpp>
#include <boost/spirit/repository/include/qi_kwd.hpp>
#include <boost/variant.hpp>
#include <tuple>

#include "center_selection_grammar_ast.hpp"
#include "dsl/grammar.hpp"
#include "utils/common.hpp"

namespace qi = boost::spirit::qi;
namespace fusion = boost::fusion;
namespace phoenix = boost::phoenix;

template <typename Iterator, typename Skipper>
struct CenterGrammar : qi::grammar<Iterator, CenterRuleNode(), Skipper> {
    CenterGrammar();

    qi::rule<Iterator, CenterRuleNode(), Skipper> expr;
    qi::rule<Iterator, Atom::Node(), Skipper> mask;
    Grammar<Iterator, Skipper> maskParser;
};

template <typename Iterator, typename Skipper>
CenterGrammar<Iterator, Skipper>::CenterGrammar() : CenterGrammar::base_type(expr) {
    using boost::spirit::repository::qi::distinct;
    using phoenix::bind;
    using phoenix::new_;
    using phoenix::ref;
    using phoenix::val;
    using qi::_1;
    using qi::_2;
    using qi::_3;
    using qi::_4;
    using qi::_val;
    using qi::alnum;
    using qi::alpha;
    using qi::as_string;
    using qi::eol;
    using qi::eps;
    using qi::fail;
    using qi::lexeme;
    using qi::lit;
    using qi::on_error;
    using qi::uint_;
    using qi::ascii::char_;

    mask = (("{" >> maskParser >> "}") | maskParser)[_val = _1];

    expr = (distinct(char_("a-zA-Z_0-9"))["com"] >> distinct(char_("a-zA-Z_0-9"))["of"] >>
            mask[_val = make_shared_<MassCenterRuleNode>(_1)]) |
           (distinct(char_("a-zA-Z_0-9"))["geom"] >> distinct(char_("a-zA-Z_0-9"))["of"] >>
            mask[_val = make_shared_<GeomCenterRuleNode>(_1)]) |
           mask[_val = make_shared_<NoopRuleNode>(_1)] | lit("quit")[_val = make_shared_<QuitRuleNode>()] |
           lit("help")[_val = make_shared_<HelpRuleNode>()];

    on_error<fail>(expr, std::cout << val("Error! Expecting") << _4 << val(" here: \"")
                                   << phoenix::construct<std::string>(_3, _2) << val("\"") << std::endl);
}

template <typename Iterator, typename Skipper>
std::tuple<CenterRuleNode, std::string> input_atom_selection(const CenterGrammar<Iterator, Skipper> &grammar,
                                                             const std::string &prompt) {
    for (;;) {
        CenterRuleNode mask;
        std::string input_string = input(prompt);
        boost::trim(input_string);
        if (input_string.empty()) continue;
        auto it = input_string.begin();

        bool status = qi::phrase_parse(it, input_string.end(), grammar, qi::ascii::space, mask);

        if (!(status and (it == input_string.end()))) {
            std::cout << "error-pos : " << std::endl;
            std::cout << input_string << std::endl;
            for (auto iter = input_string.begin(); iter != it; ++iter) std::cout << " ";
            std::cout << "^" << std::endl;

            continue;
        }
        return std::make_tuple(mask, input_string);
    }
}

inline std::string selectCentergroup(CenterRuleNode &ids, const std::string &prompt) {
    namespace qi = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;
    using ascii::char_;

    CenterGrammar<std::string::iterator, qi::ascii::space_type> grammar;

    std::string input;
    std::tie(ids, input) = input_atom_selection(grammar, prompt);
    return input;
}

#endif  // TINKER_CENTER_SELECTION_GRAMMAR_HPP
