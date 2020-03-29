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
#include "utils/ProgramConfiguration.hpp"
#include "utils/common.hpp"

namespace qi = boost::spirit::qi;
namespace fusion = boost::fusion;
namespace phoenix = boost::phoenix;

template <typename Iterator, typename Skipper> struct Grammar : qi::grammar<Iterator, Atom::Node(), Skipper> {
    Grammar();

    qi::rule<Iterator, boost::variant<fusion::vector<uint, boost::optional<std::pair<uint, int>>>, std::string>(),
             Skipper>
        select_item_rule;

    qi::rule<Iterator, Atom::Node(), Skipper> residue_select_rule;
    qi::rule<Iterator, Atom::Node(), Skipper> molecule_select_rule;
    qi::rule<Iterator, Atom::Node(), Skipper> nametype_select_rule;
    qi::rule<Iterator, Atom::Node(), Skipper> select_rule;
    qi::rule<Iterator, Atom::Node(), Skipper> term, factor;
    qi::rule<Iterator, Atom::Node(), Skipper> expr;
    qi::rule<Iterator, Atom::Node(), Skipper> maskParser, root;

    qi::rule<Iterator, std::string(), Skipper> str_with_wildcard;

    qi::rule<Iterator, Atom::Node(), Skipper> macro_rule, custum_macro_rule;

    std::vector<qi::rule<Iterator, Atom::Node(), Skipper>> r;

    inline static auto protein = Atom::select_ranges{
        "ABU",   "ACE",   "AIB",   "ALA",   "ARG",   "ARGN",  "ASN",   "ASN1", "ASP",   "ASP1",  "ASPH",  "ASPP",
        "ASH",   "CT3",   "CYS",   "CYS1",  "CYS2",  "CYSH",  "DALA",  "GLN",  "GLU",   "GLUH",  "GLUP",  "GLH",
        "GLY",   "HIS",   "HIS1",  "HISA",  "HISB",  "HISH",  "HISD",  "HISE", "HISP",  "HSD",   "HSE",   "HSP",
        "HYP",   "ILE",   "LEU",   "LSN",   "LYS",   "LYSH",  "MELEU", "MET",  "MEVAL", "NAC",   "NME",   "NHE",
        "NH2",   "PHE",   "PHEH",  "PHEU",  "PHL",   "PRO",   "SER",   "THR",  "TRP",   "TRPH",  "TRPU",  "TYR",
        "TYRH",  "TYRU",  "VAL",   "PGLU",  "HID",   "HIE",   "HIP",   "LYP",  "LYN",   "CYN",   "CYM",   "CYX",
        "DAB",   "ORN",   "NALA",  "NGLY",  "NSER",  "NTHR",  "NLEU",  "NILE", "NVAL",  "NASN",  "NGLN",  "NARG",
        "NHID",  "NHIE",  "NHIP",  "NHISD", "NHISE", "NHISH", "NTRP",  "NPHE", "NTYR",  "NGLU",  "NASP",  "NLYS",
        "NORN",  "NDAB",  "NLYSN", "NPRO",  "NHYP",  "NCYS",  "NCYS2", "NMET", "NASPH", "NGLUH", "CALA",  "CGLY",
        "CSER",  "CTHR",  "CLEU",  "CILE",  "CVAL",  "CASN",  "CGLN",  "CARG", "CHID",  "CHIE",  "CHIP",  "CHISD",
        "CHISE", "CHISH", "CTRP",  "CPHE",  "CTYR",  "CGLU",  "CASP",  "CLYS", "CORN",  "CDAB",  "CLYSN", "CPRO",
        "CHYP",  "CCYS",  "CCYS2", "CMET",  "CASPH", "CGLUH"};

    inline static auto dna = Atom::select_ranges{"DA5", "DA", "DA3", "DAN", "DT5", "DT", "DT3", "DTN",
                                                 "DG5", "DG", "DG3", "DGN", "DC5", "DC", "DC3", "DCN"};

    inline static auto rna =
        Atom::select_ranges{"A",   "U",  "C",   "G",   "RA5", "RA", "RA3", "RAN", "RU5", "RU",  "RU3", "RUN",
                            "RG5", "RG", "RG3", "RGN", "RC5", "RC", "RC3", "RCN", "RT5", "RT3", "RTN"};

    inline static auto water = Atom::select_ranges{"SOL", "WAT", "HOH", "OHH", "TIP", "T3P", "T4P", "T5P", "T3H"};
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

    residue_select_rule = ':' > (select_item_rule % ',')[_val = make_shared_<Atom::residue_name_nums>(_1)];

    molecule_select_rule =
        '$' > ((uint_ >> -('-' > uint_ >> -('#' > int_))) % ',')[_val = make_shared_<Atom::molecule_nums>(_1)];

    nametype_select_rule = '@' > (('%' > (select_item_rule % ',')[_val = make_shared_<Atom::atom_types>(_1)]) |
                                  ('/' > (str_with_wildcard % ',')[_val = make_shared_<Atom::atom_element_names>(_1)]) |
                                  (select_item_rule % ',')[_val = make_shared_<Atom::atom_name_nums>(_1)]);

#define DISTINCT(x) distinct(char_("a-zA-Z_0-9") | char_("-+"))[x]

    macro_rule =
        (DISTINCT("System") | DISTINCT("All"))[_val = make_shared_<Atom::atom_name_nums>(Atom::select_ranges{"*"})]

        | DISTINCT("Protein")[_val = make_shared_<Atom::residue_name_nums>(protein)]

        | DISTINCT("Protein-H")[_val = make_shared_<Atom::Operator>(
                                    Atom::Op::AND, make_shared_<Atom::residue_name_nums>(protein),
                                    make_shared_<Atom::Operator>(Atom::Op::NOT, make_shared_<Atom::atom_name_nums>(
                                                                                    Atom::select_ranges{"H*"})))] |
        DISTINCT("Backbone")[_val = make_shared_<Atom::Operator>(
                                 Atom::Op::AND, make_shared_<Atom::residue_name_nums>(protein),
                                 make_shared_<Atom::atom_name_nums>(Atom::select_ranges{"CA", "C", "N"}))]

        | DISTINCT("MainChain")[_val = make_shared_<Atom::Operator>(
                                    Atom::Op::AND, make_shared_<Atom::residue_name_nums>(protein),
                                    make_shared_<Atom::atom_name_nums>(Atom::select_ranges{"CA", "C", "N", "O"}))]

        | DISTINCT(
              "MainChain+Cb")[_val = make_shared_<Atom::Operator>(
                                  Atom::Op::AND, make_shared_<Atom::residue_name_nums>(protein),
                                  make_shared_<Atom::atom_name_nums>(Atom::select_ranges{"CA", "C", "N", "O", "CB"}))]

        |
        DISTINCT("MainChain+H")[_val = make_shared_<Atom::Operator>(
                                    Atom::Op::AND, make_shared_<Atom::residue_name_nums>(protein),
                                    make_shared_<Atom::atom_name_nums>(Atom::select_ranges{"CA", "C", "N", "O", "H"}))]

        | DISTINCT("C-alpha")[_val = make_shared_<Atom::Operator>(
                                  Atom::Op::AND, make_shared_<Atom::residue_name_nums>(protein),
                                  make_shared_<Atom::atom_name_nums>(Atom::select_ranges{"CA"}))]

        | DISTINCT("SideChain")[_val = make_shared_<Atom::Operator>(
                                    Atom::Op::AND, make_shared_<Atom::residue_name_nums>(protein),
                                    make_shared_<Atom::Operator>(Atom::Op::NOT,
                                                                 make_shared_<Atom::atom_name_nums>(
                                                                     Atom::select_ranges{"CA", "C", "O", "N", "H"})))]

        | DISTINCT("SideChain-H")[_val = make_shared_<Atom::Operator>(
                                      Atom::Op::AND,
                                      make_shared_<Atom::Operator>(
                                          Atom::Op::AND, make_shared_<Atom::residue_name_nums>(protein),
                                          make_shared_<Atom::Operator>(
                                              Atom::Op::NOT, make_shared_<Atom::atom_name_nums>(
                                                                 Atom::select_ranges{"CA", "C", "O", "N", "H"}))),
                                      make_shared_<Atom::Operator>(Atom::Op::NOT, make_shared_<Atom::atom_name_nums>(
                                                                                      Atom::select_ranges{"H*"})))]

        | DISTINCT("DNA")[_val = make_shared_<Atom::residue_name_nums>(dna)]

        | DISTINCT("DNA-H")[_val = make_shared_<Atom::Operator>(
                                Atom::Op::AND, make_shared_<Atom::residue_name_nums>(dna),
                                make_shared_<Atom::Operator>(
                                    Atom::Op::NOT, make_shared_<Atom::atom_name_nums>(Atom::select_ranges{"H*"})))]

        | DISTINCT("RNA")[_val = make_shared_<Atom::residue_name_nums>(rna)]

        | DISTINCT("RNA-H")[_val = make_shared_<Atom::Operator>(
                                Atom::Op::AND, make_shared_<Atom::residue_name_nums>(rna),
                                make_shared_<Atom::Operator>(
                                    Atom::Op::NOT, make_shared_<Atom::atom_name_nums>(Atom::select_ranges{"H*"})))]

        | DISTINCT("Water")[_val = make_shared_<Atom::residue_name_nums>(water)];

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

    term = ("!" > factor[_val = make_shared_<Atom::Operator>(Atom::Op::NOT, _1)]) | factor[_val = _1];

    expr = term[_val = _1] > *("&" > term[_val = make_shared_<Atom::Operator>(Atom::Op::AND, _val, _1)]);

    maskParser = expr[_val = _1] > *("|" > expr[_val = make_shared_<Atom::Operator>(Atom::Op::OR, _val, _1)]);

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
