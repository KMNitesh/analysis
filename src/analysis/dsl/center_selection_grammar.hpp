//
// Created by xiamr on 6/1/19.
//

#ifndef TINKER_CENTER_SELECTION_GRAMMAR_HPP
#define TINKER_CENTER_SELECTION_GRAMMAR_HPP

#include "center_selection_grammar_ast.hpp"
#include "dsl/grammar.hpp"
#include "utils/common.hpp"
#include <boost/algorithm/string.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/phoenix.hpp>
#include <boost/regex.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/support_line_pos_iterator.hpp>
#include <boost/spirit/repository/include/qi_distinct.hpp>
#include <boost/spirit/repository/include/qi_iter_pos.hpp>
#include <boost/spirit/repository/include/qi_kwd.hpp>
#include <boost/variant.hpp>
#include <tuple>

namespace qi = boost::spirit::qi;
namespace fusion = boost::fusion;
namespace phoenix = boost::phoenix;

template <typename Iterator, typename Skipper> struct CenterGrammar : qi::grammar<Iterator, CenterRuleNode(), Skipper> {
    CenterGrammar();

    qi::rule<Iterator, CenterRuleNode(), Skipper> expr, root;
    qi::rule<Iterator, AmberMask(), Skipper> mask;
    Grammar<Iterator, Skipper> maskParser;
};

template <typename Iterator, typename Skipper>
CenterGrammar<Iterator, Skipper>::CenterGrammar() : CenterGrammar::base_type(root) {
    using boost::spirit::repository::qi::distinct;
    using phoenix::bind;
    using phoenix::construct;
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
    using qi::eoi;
    using qi::eps;
    using qi::lexeme;
    using qi::lit;
    using qi::uint_;
    using qi::ascii::char_;

    mask = (("[" >> maskParser >> "]") | maskParser)[_val = _1];
    mask.name("mask");

    expr = (DISTINCT("com") > mask[_val = construct<MassCenterRuleNode>(_1)]) |
           (DISTINCT("geom") > mask[_val = construct<GeomCenterRuleNode>(_1)]) |
           ((DISTINCT("eda") > mask > mask)[_val = construct<EDARuleNode>(_1, _2)]) |
           ((DISTINCT("bond") > mask > mask)[_val = construct<BondRuleNode>(_1, _2)]) |
           ((DISTINCT("angle") > mask > mask > mask)[_val = construct<AngleRuleNode>(_1, _2, _3)]) |
           ((DISTINCT("dihedral") > mask > mask > mask > mask)[_val = construct<DihedralRuleNode>(_1, _2, _3, _4)]) |
           (DISTINCT("energy") > mask[_val = construct<BondedEnergyRuleNode>(_1)]) |
           (DISTINCT("help") >> -(as_string[lexeme[+(char_ - qi::ascii::space)]]))[_val = construct<HelpRuleNode>(_1)] |
           DISTINCT("quit")[_val = construct<QuitRuleNode>()] | mask[_val = construct<NoopRuleNode>(_1)];
    root = eps > expr;
}

std::tuple<CenterRuleNode, std::string>
input_atom_selection(const CenterGrammar<std::string::iterator, qi::ascii::space_type> &grammar,
                     const std::string &prompt) {

    for (;;) {
        CenterRuleNode mask;
        std::string input_string;
        if (isatty(STDIN_FILENO) == 1) {
            char *buf = readline::readline(prompt.c_str());
            input_string = buf;
            free(buf);
        } else {
            input_string = input(prompt);
        }
        boost::trim(input_string);
        if (input_string.empty())
            continue;

        std::string::iterator begin{std::begin(input_string)}, it{begin}, end{std::end(input_string)};
        try {
            qi::phrase_parse(it, end, grammar > qi::eoi, boost::spirit::ascii::space, mask);
            return {mask, input_string};
        } catch (const qi::expectation_failure<std::string::iterator> &x) {
            std::cerr << "Grammar Parse Failure ! Expecting : " << x.what_ << '\n';
            auto column = boost::spirit::get_column(begin, x.first);
            std::string pos = " (column: " + std::to_string(column) + ")";
            std::cerr << pos << ">>>>" << input_string << "<<<<\n";
            std::cerr << std::string(column + pos.size() + 3, ' ') << "^~~~ here\n";
        }
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

#endif // TINKER_CENTER_SELECTION_GRAMMAR_HPP
