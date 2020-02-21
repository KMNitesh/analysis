//
// Created by xiamr on 3/17/19.
//

#ifndef TINKER_GRAMMER_HPP
#define TINKER_GRAMMER_HPP


#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/repository/include/qi_distinct.hpp>
#include <boost/phoenix.hpp>
#include <boost/variant.hpp>
#include <boost/optional.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/phoenix/function/adapt_function.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>

#include "utils/common.hpp"
#include "data_structure/atom.hpp"

namespace qi = boost::spirit::qi;
namespace fusion = boost::fusion;
namespace phoenix = boost::phoenix;


template<typename Iterator, typename Skipper>
struct Grammar : qi::grammar<Iterator, Atom::Node(), Skipper> {
    Grammar();


    qi::rule<Iterator, boost::variant<fusion::vector<uint, boost::optional<std::pair<uint, int>>>,
            std::string>(), Skipper> select_item_rule;

    qi::rule<Iterator, Atom::Node(), Skipper> residue_select_rule;
    qi::rule<Iterator, Atom::Node(), Skipper> molecule_select_rule;
    qi::rule<Iterator, Atom::Node(), Skipper> nametype_select_rule;
    qi::rule<Iterator, Atom::Node(), Skipper> select_rule;
    qi::rule<Iterator, Atom::Node(), Skipper> term, factor;
    qi::rule<Iterator, Atom::Node(), Skipper> expr;
    qi::rule<Iterator, Atom::Node(), Skipper> maskParser;

    qi::rule<Iterator, std::string(), Skipper> str_with_wildcard;

    qi::rule<Iterator, Atom::Node(), Skipper> macro_rule;

};


BOOST_PHOENIX_ADAPT_FUNCTION(std::string, replace_all_copy, boost::replace_all_copy, 3)


template<typename Iterator, typename Skipper>
Grammar<Iterator, Skipper>::Grammar() : Grammar::base_type(maskParser, "mask") {
    using qi::uint_;
    using qi::int_;
    using qi::eps;
    using qi::as_string;
    using qi::lexeme;
    using qi::alpha;
    using qi::alnum;
    using qi::ascii::char_;
    using qi::_val;
    using qi::_1;
    using qi::_2;
    using qi::_3;
    using qi::_4;
    using phoenix::new_;
    using phoenix::val;
    using phoenix::bind;
    using phoenix::ref;
    using qi::on_error;
    using qi::fail;
    using boost::spirit::repository::qi::distinct;

    str_with_wildcard = as_string[lexeme[+(alnum | char_("*?="))]][_val = replace_all_copy(_1, "=", "*")];

    select_item_rule = (str_with_wildcard >> -("-" >> uint_ >> -("#" >> int_)))[(
            [](auto &attr, auto &context, bool &pass) {
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

                        fusion::at_c<0>(
                                context.attributes) = fusion::vector<uint, boost::optional<std::pair<uint, int>>>(
                                num, std::make_pair(fusion::at_c<0>(fusion::at_c<1>(attr).get()), step));
                    } else {
                        fusion::at_c<0>(
                                context.attributes) = fusion::vector<uint, boost::optional<std::pair<uint, int>>>(
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

    residue_select_rule = ':' >> (select_item_rule % ',')[_val = make_shared_<Atom::residue_name_nums>(_1)];

    molecule_select_rule =
            '$' >> ((uint_ >> -('-' >> uint_ >> -('#' >> int_))) % ',')[_val = make_shared_<Atom::molecule_nums>(_1)];

    nametype_select_rule =
            '@' >> ('%' >> (select_item_rule % ',')[_val = make_shared_<Atom::atom_types>(_1)]
                    | '/' >> (str_with_wildcard % ',')[_val = make_shared_<Atom::atom_element_names>(_1)]
                    | (select_item_rule % ',')[_val = make_shared_<Atom::atom_name_nums>(_1)]);

#define DISTINCT(x) distinct(char_("a-zA-Z_0-9") | char_("-"))[x]

    static auto protein = std::vector<boost::variant<Atom::numItemType, std::string>>{
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "HYP", "ILE", "LLE",
            "LEU", "LYS", "MET", "PHE", "PRO", "GLP", "SER", "THR", "TRP", "TYR", "VAL"};

    static auto dna = std::vector<boost::variant<Atom::numItemType, std::string>>{
            "DA5", "DA", "DA3", "DAN",
            "DT5", "DT", "DT3", "DTN",
            "DG5", "DG", "DG3", "DGN",
            "DC5", "DC", "DC3", "DCN"};

    static auto rna = std::vector<boost::variant<Atom::numItemType, std::string>>{
            "RA5", "RA", "RA3", "RAN",
            "RU5", "RU", "RU3", "RUN",
            "RG5", "RG", "RG3", "RGN",
            "RC5", "RC", "RC3", "RCN"};

    macro_rule = (DISTINCT("System") | DISTINCT("All"))
                 [_val = make_shared_<Atom::atom_element_names>(std::vector<std::string>{"*"})]

                 | DISTINCT("Protein")[_val = make_shared_<Atom::residue_name_nums>(protein)]

                 | DISTINCT("Protein-H")[_val = make_shared_<Atom::Operator>(
            Atom::Op::AND,
            make_shared_<Atom::residue_name_nums>(protein),
            make_shared_<Atom::Operator>(
                    Atom::Op::NOT,
                    make_shared_<Atom::atom_element_names>(std::vector<std::string>{"H*"})))]
                 | DISTINCT("Backbone")[_val = make_shared_<Atom::Operator>(
            Atom::Op::AND,
            make_shared_<Atom::residue_name_nums>(protein),
            make_shared_<Atom::atom_element_names>(std::vector<std::string>{"CA", "C", "O", "N", "H"}))]

                 | DISTINCT("MainChain")[_val = make_shared_<Atom::Operator>(
            Atom::Op::AND,
            make_shared_<Atom::residue_name_nums>(protein),
            make_shared_<Atom::atom_element_names>(std::vector<std::string>{"CA", "C", "N"}))]

                 | DISTINCT("C-alpha")[_val = make_shared_<Atom::Operator>(
            Atom::Op::AND,
            make_shared_<Atom::residue_name_nums>(protein),
            make_shared_<Atom::atom_element_names>(std::vector<std::string>{"CA"}))]

                 | DISTINCT("SideChain")[_val = make_shared_<Atom::Operator>(
            Atom::Op::AND,
            make_shared_<Atom::residue_name_nums>(protein),
            make_shared_<Atom::Operator>(
                    Atom::Op::NOT,
                    make_shared_<Atom::atom_element_names>(std::vector<std::string>{"CA", "C", "O", "N", "H"})))]

                 | DISTINCT("SideChain-H")[_val = make_shared_<Atom::Operator>(
            Atom::Op::AND,
            make_shared_<Atom::Operator>(
                    Atom::Op::AND,
                    make_shared_<Atom::residue_name_nums>(protein),
                    make_shared_<Atom::Operator>(
                            Atom::Op::NOT,
                            make_shared_<Atom::atom_element_names>(
                                    std::vector<std::string>{"CA", "C", "O", "N", "H"}))),
            make_shared_<Atom::Operator>(
                    Atom::Op::NOT, make_shared_<Atom::atom_element_names>(std::vector<std::string>{"H*"})))]

                 | DISTINCT("DNA")[_val = make_shared_<Atom::residue_name_nums>(dna)]

                 | DISTINCT("DNA-H")[_val = make_shared_<Atom::Operator>(
            Atom::Op::AND,
            make_shared_<Atom::residue_name_nums>(dna),
            make_shared_<Atom::Operator>(
                    Atom::Op::NOT, make_shared_<Atom::atom_element_names>(std::vector<std::string>{"H*"})))]

                 | DISTINCT("RNA")[_val = make_shared_<Atom::residue_name_nums>(rna)]

                 | DISTINCT("RNA-H")[_val = make_shared_<Atom::Operator>(
            Atom::Op::AND,
            make_shared_<Atom::residue_name_nums>(rna),
            make_shared_<Atom::Operator>(
                    Atom::Op::NOT,
                    make_shared_<Atom::atom_element_names>(std::vector<std::string>{"H*"})))]

                 | DISTINCT("Water")[_val = make_shared_<Atom::residue_name_nums>(
            std::vector<boost::variant<Atom::numItemType, std::string>>{"WAT", "SOL"})];

    select_rule = macro_rule | residue_select_rule | molecule_select_rule | nametype_select_rule;

    factor = "(" >> maskParser >> ")" | select_rule;

    term = "!" >> factor[_val = make_shared_<Atom::Operator>(Atom::Op::NOT, _1)] | factor[_val = _1];

    expr = term[_val = _1] >> *("&" >> term[_val = make_shared_<Atom::Operator>(Atom::Op::AND, _val, _1)]);

    maskParser = expr[_val = _1] >> *("|" >> expr[_val = make_shared_<Atom::Operator>(Atom::Op::OR, _val, _1)]);

    str_with_wildcard.name("str_with_wildcard");
    select_item_rule.name("select_item_rule");
    residue_select_rule.name("residue_select_rule");
    molecule_select_rule.name("molecule_select_rule");
    nametype_select_rule.name("nametype_select_rule");
    select_rule.name("select_rule");
    factor.name("factor");
    term.name("term");
    expr.name("expr");
    maskParser.name("mask");

    on_error<fail>(
            maskParser,
            std::cout
                    << val("Error! Expecting")
                    << _4
                    << val(" here: \"")
                    << phoenix::construct<std::string>(_3, _2)
                    << val("\"")
                    << std::endl
    );

}


#endif //TINKER_GRAMMER_HPP
