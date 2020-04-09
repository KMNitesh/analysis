//
// Created by xiamr on 3/17/19.
//

#ifndef TINKER_GRAMMER_HPP
#define TINKER_GRAMMER_HPP

#include <boost/algorithm/string.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>
#include <boost/phoenix.hpp>
#include <boost/phoenix/function/adapt_function.hpp>
#include <boost/regex.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/repository/include/qi_distinct.hpp>
#include <boost/variant.hpp>

#include "data_structure/atom.hpp"
#include "dsl/AmberMask.hpp"
#include "dsl/MacroRules.hpp"
#include "utils/ProgramConfiguration.hpp"
#include "utils/common.hpp"

namespace qi = boost::spirit::qi;
namespace fusion = boost::fusion;
namespace phoenix = boost::phoenix;

template <typename Iterator, typename Skipper> struct Grammar : qi::grammar<Iterator, AmberMask(), Skipper> {
    Grammar();

    qi::rule<Iterator,
             boost::variant<fusion::vector<uint, boost::optional<std::pair<uint, int>>>, AmberMaskAST::Name>(), Skipper>
        select_item_rule;
    qi::rule<Iterator, boost::optional<std::pair<uint, int>>(), Skipper> num_range;

    qi::rule<Iterator, AmberMask(), Skipper> residue_select_rule;
    qi::rule<Iterator, AmberMask(), Skipper> molecule_select_rule;
    qi::rule<Iterator, AmberMask(), Skipper> nametype_select_rule;
    qi::rule<Iterator, AmberMask(), Skipper> select_rule;
    qi::rule<Iterator, AmberMask(), Skipper> term, factor;
    qi::rule<Iterator, AmberMask(), Skipper> expr;
    qi::rule<Iterator, AmberMask(), Skipper> maskParser, root;

    qi::rule<Iterator, AmberMaskAST::Name(), Skipper> str_with_wildcard;

    qi::rule<Iterator, AmberMask(), Skipper> macro_rule, custum_macro_rule;

    std::vector<qi::rule<Iterator, AmberMask(), Skipper>> r;
};

BOOST_PHOENIX_ADAPT_FUNCTION(std::string, replace_all_copy, boost::replace_all_copy, 3)

template <typename Iterator, typename Skipper>
Grammar<Iterator, Skipper>::Grammar() : Grammar::base_type(root, "mask") {
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
    using qi::eps;
    using qi::int_;
    using qi::lexeme;
    using qi::uint_;
    using qi::ascii::char_;

    str_with_wildcard = as_string[lexeme[+(alnum | char_("*?="))]][_val = replace_all_copy(_1, "=", "*")];

    select_item_rule =
        (str_with_wildcard >> -("-" > uint_ >> -("#" > int_)))[([](auto &attr, auto &context, bool &pass) {
            try {
                uint num = boost::lexical_cast<uint>(fusion::at_c<0>(attr));

                if (fusion::at_c<1>(attr)) {
                    int step = 1;
                    if (fusion::at_c<1>(fusion::at_c<1>(attr).get())) {
                        step = fusion::at_c<1>(fusion::at_c<1>(attr).get()).get();
                        if (step == 0) {
                            pass = false;
                            return;
                        }
                    }

                    fusion::at_c<0>(context.attributes) = fusion::vector<uint, boost::optional<std::pair<uint, int>>>(
                        num, std::make_pair(fusion::at_c<0>(fusion::at_c<1>(attr).get()), step));
                } else {
                    fusion::at_c<0>(context.attributes) = fusion::vector<uint, boost::optional<std::pair<uint, int>>>(
                        num, boost::optional<std::pair<uint, int>>());
                }

            } catch (boost::bad_lexical_cast &) {
                if (fusion::at_c<1>(attr)) {
                    pass = false;
                    return;
                } else {
                    fusion::at_c<0>(context.attributes) = fusion::at_c<0>(attr);
                }
            }
        })];

    num_range = ("-" > (uint_ >> -("#" > int_)))[([](auto &attr, auto &context, bool &pass) {
        auto step = fusion::at_c<1>(attr) ? fusion::at_c<1>(attr).get() : 1;
        if (step == 0) {
            pass = false;
            return;
        }
        fusion::at_c<0>(context.attributes) = std::make_pair(fusion::at_c<0>(attr), step);
    })];

    residue_select_rule =
        ':' > ('^' >> ((uint_ >> -num_range) % ",")[_val = make_shared_<AmberMaskAST::real_residue_nums>(_1)] |
               (select_item_rule % ',')[_val = make_shared_<AmberMaskAST::residue_name_nums>(_1)]);

    molecule_select_rule =
        '$' > ((uint_ >> -('-' > uint_ >> -('#' > int_))) % ',')[_val = make_shared_<AmberMaskAST::molecule_nums>(_1)];

    nametype_select_rule =
        '@' > (('%' > (select_item_rule % ',')[_val = make_shared_<AmberMaskAST::atom_types>(_1)]) |
               ('/' > (str_with_wildcard % ',')[_val = make_shared_<AmberMaskAST::atom_element_names>(_1)]) |
               (select_item_rule % ',')[_val = make_shared_<AmberMaskAST::atom_name_nums>(_1)]);

#define DISTINCT(x) distinct(char_("a-zA-Z_0-9") | char_("-+"))[x]

    macro_rule = (DISTINCT("System") | DISTINCT("All"))[_val = AMBERMASK::System] |
                 DISTINCT("Protein")[_val = AMBERMASK::Protein] | DISTINCT("Protein-H")[_val = AMBERMASK::Protein_H] |
                 DISTINCT("Backbone")[_val = AMBERMASK::Backbone] | DISTINCT("MainChain")[_val = AMBERMASK::MainChain] |
                 DISTINCT("MainChain+Cb")[_val = AMBERMASK::MainChain_plus_Cb] |
                 DISTINCT("MainChain+H")[_val = AMBERMASK::MainChain_plus_H] |
                 DISTINCT("C-alpha")[_val = AMBERMASK::C_alpha] | DISTINCT("SideChain")[_val = AMBERMASK::SideChain] |
                 DISTINCT("SideChain-H")[_val = AMBERMASK::SideChain_H] | DISTINCT("DNA")[_val = AMBERMASK::DNA] |
                 DISTINCT("DNA-H")[_val = AMBERMASK::DNA_H] | DISTINCT("RNA")[_val = AMBERMASK::RNA] |
                 DISTINCT("RNA-H")[_val = AMBERMASK::RNA_H] | DISTINCT("Water")[_val = AMBERMASK::Water];

    if (program_configuration) {
        const auto &all_macros = program_configuration->get_macro_mask();
        r.resize(all_macros.size() + 1);
        r[0] = eps(false);

        for (std::size_t i = 0; i < all_macros.size(); ++i) {
            const auto &[name, macro] = all_macros[i];
            r[i + 1] = DISTINCT(name)[_val = macro] | r[i][_val = _1];
        }
        custum_macro_rule = r[all_macros.size()];
    } else {
        custum_macro_rule = eps(false);
    }

    select_rule = custum_macro_rule | macro_rule | residue_select_rule | molecule_select_rule | nametype_select_rule;

    factor = ("(" > maskParser > ")") | select_rule;

    term = ("!" > factor[_val = make_shared_<AmberMaskAST::Operator>(AmberMaskAST::Op::NOT, _1)]) | factor[_val = _1];

    expr =
        term[_val = _1] > *("&" > term[_val = make_shared_<AmberMaskAST::Operator>(AmberMaskAST::Op::AND, _val, _1)]);

    maskParser =
        expr[_val = _1] > *("|" > expr[_val = make_shared_<AmberMaskAST::Operator>(AmberMaskAST::Op::OR, _val, _1)]);

    root = eps > maskParser;

    str_with_wildcard.name("str_with_wildcard");
    select_item_rule.name("select_item_rule");
    residue_select_rule.name("residue_select_rule");
    molecule_select_rule.name("molecule_select_rule");
    nametype_select_rule.name("nametype_select_rule");
    select_rule.name("mask");
    factor.name("mask");
    term.name("mask");
    expr.name("mask");
    maskParser.name("mask");
}

#endif // TINKER_GRAMMER_HPP
