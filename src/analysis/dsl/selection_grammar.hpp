//
// Created by xiamr on 6/1/19.
//

#ifndef TINKER_SELECTION_GRAMMAR_HPP
#define TINKER_SELECTION_GRAMMAR_HPP

#include "dsl/grammar.hpp"
#include "selection_grammar_ast.hpp"
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

template <typename Iterator, typename Skipper> struct CenterGrammar : qi::grammar<Iterator, SelectionAST(), Skipper> {
    CenterGrammar();

    qi::rule<Iterator, std::string(), Skipper> id, id_without_star, string_;
    qi::rule<Iterator, boost::any(), Skipper> value;
    qi::rule<Iterator, Condition(), qi::locals<LogicalOperator::Op>, Skipper> primary_condition;
    qi::rule<Iterator, Condition(), Skipper> like_condition, regex_condition, between_condition, in_condition;
    qi::rule<Iterator, Condition(), Skipper> combine_condtion_factor, combine_condtion_term, combine_condtion_expr,
        combine_condtion_root;

    qi::rule<Iterator, Condition(), Skipper> where_clause;
    qi::rule<Iterator, OrderBy(), qi::locals<OrderBy::Order>, Skipper> orderby_clause;
    qi::rule<Iterator, uint(), Skipper> limit_clause;
    qi::rule<Iterator, SelectStmt(),
             qi::locals<bool, boost::optional<Condition>, boost::optional<OrderBy>, boost::optional<uint>>, Skipper>
        select_stmt;

    qi::rule<Iterator, SelectionAST(), Skipper> expr, root;
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
    using qi::_a;
    using qi::_b;
    using qi::_c;
    using qi::_d;
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
    using qi::ascii::no_case;

    id_without_star %= as_string[lexeme[(alpha | char_("_")) >> *(alnum | char_("_"))]];
    id %= id_without_star | char_('*');
    id.name("id");
    string_ %= as_string[lexeme["'" >> *(char_ - "'") >> "'"]];
    value %= lexeme[qi::int_ >> !char_('.')] | qi::double_ | string_;
    value.name("value");

    primary_condition =
        (id_without_star >> (lit("<")[_a = LogicalOperator::Op::Less] | lit("<=")[_a = LogicalOperator::Op::LessEqual] |
                             lit("=")[_a = LogicalOperator::Op::Equal] | lit(">")[_a = LogicalOperator::Op::Great] |
                             lit(">=")[_a = LogicalOperator::Op::GreatEqual]) >>
         value)[_val = make_shared_<LogicalOperator>(_1, _a, _2)];

    like_condition =
        (id_without_star >> DISTINCT(no_case["like"]) >> string_)[_val = make_shared_<LikeOperator>(_1, _2)];

    regex_condition =
        (id_without_star >> DISTINCT(no_case["regex"]) >> string_)[_val = make_shared_<RegexOperator>(_1, _2)];

    between_condition = (id_without_star >> DISTINCT(no_case["between"]) >> qi::double_ >> DISTINCT(no_case["and"]) >>
                         qi::double_)[_val = make_shared_<BetweenOperator>(_1, _2, _3)];

    in_condition = (id_without_star >> DISTINCT(no_case["in"]) >> '(' >> (string_ % ',') >>
                    ')')[_val = make_shared_<InOperator>(_1, _2)];

    combine_condtion_factor %= ('(' > combine_condtion_root > ')') | primary_condition | like_condition |
                               regex_condition | between_condition | in_condition;

    combine_condtion_term =
        (DISTINCT(no_case["not"]) >
         combine_condtion_factor[_val = make_shared_<CombineOperator>(CombineOperator::Op::NOT, _1)]) |
        combine_condtion_factor[_val = _1];

    combine_condtion_expr =
        combine_condtion_term[_val = _1] >
        *(DISTINCT(no_case["and"]) >
          combine_condtion_term[_val = make_shared_<CombineOperator>(CombineOperator::Op::AND, _val, _1)]);

    combine_condtion_root =
        combine_condtion_expr[_val = _1] >
        *(DISTINCT(no_case["or"]) >
          combine_condtion_expr[_val = make_shared_<CombineOperator>(CombineOperator::Op::OR, _val, _1)]);

    primary_condition.name("condition");
    like_condition.name("condition");
    regex_conditon.name("condition");
    between_condition.name("condition");
    in_condition.name("condition");
    combine_condtion_factor.name("condition");
    combine_condtion_term.name("condition");
    combine_condtion_expr.name("condition");
    combine_condtion_root.name("condition");

    orderby_clause =
        (DISTINCT(no_case["order"]) > DISTINCT(no_case["by"])[_a = OrderBy::Order::ASC] > id_without_star >>
         -(DISTINCT(no_case["ASC"])[_a = OrderBy::Order::ASC] |
           DISTINCT(no_case["DESC"])[_a = OrderBy::Order::DESC]))[_val = construct<OrderBy>(_1, _a)];

    orderby_clause.name("order by clause");

    where_clause %= DISTINCT(no_case["where"]) > combine_condtion_root;

    where_clause.name("where_clause");

    limit_clause %= DISTINCT(no_case["limit"]) > uint_;

    select_stmt = (DISTINCT(no_case["select"])[_a = false] > -(DISTINCT(no_case["distinct"])[_a = true]) > (id % ',') >
                   DISTINCT(no_case["from"]) > mask >
                   *(where_clause[_b = _1] | orderby_clause[_c = _1] |
                     limit_clause[_d = _1]))[_val = construct<SelectStmt>(_a, _1, _2, _b, _c, _d)];

    mask = (("[" >> maskParser >> "]") | maskParser)[_val = _1];
    mask.name("mask");

    expr = (DISTINCT(no_case["com"]) > mask[_val = construct<MassCenterRuleNode>(_1)]) |
           (DISTINCT(no_case["geom"]) > mask[_val = construct<GeomCenterRuleNode>(_1)]) |
           ((DISTINCT(no_case["eda"]) > mask > mask)[_val = construct<EDARuleNode>(_1, _2)]) |
           ((DISTINCT(no_case["bond"]) > mask > mask)[_val = construct<BondRuleNode>(_1, _2)]) |
           ((DISTINCT(no_case["angle"]) > mask > mask > mask)[_val = construct<AngleRuleNode>(_1, _2, _3)]) |
           ((DISTINCT(no_case["dihedral"]) > mask > mask > mask >
             mask)[_val = construct<DihedralRuleNode>(_1, _2, _3, _4)]) |
           (DISTINCT(no_case["energy"]) > mask[_val = construct<BondedEnergyRuleNode>(_1)]) | select_stmt[_val = _1] |
           (DISTINCT(no_case["help"]) >>
            -(as_string[lexeme[+(char_ - qi::ascii::space)]]))[_val = construct<HelpRuleNode>(_1)] |
           DISTINCT(no_case["quit"])[_val = construct<QuitRuleNode>()] | mask[_val = construct<NoopRuleNode>(_1)];
    root = eps > expr;
}

std::tuple<SelectionAST, std::string>
input_atom_selection(const CenterGrammar<std::string::iterator, qi::ascii::space_type> &grammar,
                     const std::string &prompt) {

    for (;;) {
        SelectionAST mask;
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

inline std::string selectCentergroup(SelectionAST &ids, const std::string &prompt) {
    namespace qi = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;
    using ascii::char_;

    CenterGrammar<std::string::iterator, qi::ascii::space_type> grammar;

    std::string input;
    std::tie(ids, input) = input_atom_selection(grammar, prompt);
    return input;
}

#endif // TINKER_SELECTION_GRAMMAR_HPP
