//
// Created by xiamr on 3/17/19.
//

#ifndef TINKER_GRAMMER_HPP
#define TINKER_GRAMMER_HPP

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/phoenix.hpp>
#include <boost/variant.hpp>
#include <boost/optional.hpp>
#include <boost/fusion/sequence/intrinsic/at_c.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/phoenix/function/adapt_function.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>

#include "atom.hpp"

namespace qi = boost::spirit::qi;
namespace fusion = boost::fusion;
namespace phoenix = boost::phoenix;


template<typename Iterator, typename Skipper>
struct Grammar : qi::grammar<Iterator, Atom::Node(), Skipper> {
    Grammar();


    qi::rule<Iterator, boost::variant<fusion::vector<uint, boost::optional<uint>>, std::string>(), Skipper> select_item_rule;
    qi::rule<Iterator, Atom::Node(), Skipper> residue_select_rule;
    qi::rule<Iterator, Atom::Node(), Skipper> nametype_select_rule;
    qi::rule<Iterator, Atom::Node(), Skipper> select_rule;
    qi::rule<Iterator, Atom::Node(), Skipper> factor, factor2;
    qi::rule<Iterator, Atom::Node(), Skipper> term;
    qi::rule<Iterator, Atom::Node(), Skipper> expr;

    qi::rule<Iterator, std::string(), Skipper> string_with_wildcard;




};

template<typename T>
struct make_shared_f {
    template<typename... A>
    struct result {
        typedef std::shared_ptr<T> type;
    };

    template<typename... A>
    typename result<A...>::type operator()(A &&... a) const {
        return std::make_shared<T>(std::forward<A>(a)...);
    }
};


template<typename T, typename... _Args>
inline auto make_shared_(_Args &&... __args) {
    return boost::phoenix::function<make_shared_f<T>>()(std::forward<_Args>(__args)...);
}

BOOST_PHOENIX_ADAPT_FUNCTION(std::string, replace_all_copy, boost::replace_all_copy, 3)


struct fa {
    template<typename Attrib, typename Context>
    void operator()(Attrib &attr, Context &context, bool &pass) const {
       pass = true;
        try {
            uint num = boost::lexical_cast<uint>(fusion::at_c<0>(attr));
            fusion::at_c<0>(context.attributes) = fusion::vector<uint, boost::optional<uint>>(num, fusion::at_c<1>(attr));
        } catch (boost::bad_lexical_cast &) {
            if (fusion::at_c<1>(attr)) {
                pass = false;
            } else {
                fusion::at_c<0>(context.attributes) = fusion::at_c<0>(attr);
            }
        }
    }

};

template<typename Iterator, typename Skipper>
Grammar<Iterator,Skipper>::Grammar() : Grammar::base_type(expr) {
    using qi::uint_;
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


    string_with_wildcard = as_string[lexeme[+(alnum | char_("*?="))]][_val = replace_all_copy(_1, "=", "*") ];

    select_item_rule = (string_with_wildcard >> ("-" >> uint_ | eps))[fa()];

    residue_select_rule = ":" >> (select_item_rule % ",")[_val = make_shared_<Atom::residue_name_nums>(_1)];

    nametype_select_rule = "@" >>
                               ((select_item_rule % ",")[_val = make_shared_<Atom::atom_name_nums>(_1)]
                                | "%" >> (select_item_rule % ",")[_val = make_shared_<Atom::atom_types>(_1)]
                                |
                                "/" >> (string_with_wildcard % ",")[_val = make_shared_<Atom::atom_element_names>(_1)]);

    select_rule = (residue_select_rule | nametype_select_rule)[_val = _1];

    factor2 = "(" >> expr[_val = _1] >> ")" | select_rule[_val = _1];

    factor = "!" >> factor2[_val = make_shared_<Atom::Operator>(Atom::Op::NOT, _1)] | factor2[_val = _1];

    term = factor[_val = _1] >> *("&" >> factor[_val = make_shared_<Atom::Operator>(Atom::Op::AND, _val, _1)]);

    expr = term[_val = _1] >> *("|" >> term[_val = make_shared_<Atom::Operator>(Atom::Op::OR, _val, _1)]);

    string_with_wildcard.name("string_with_wildcard");
    select_item_rule.name("select_item_rule");
    residue_select_rule.name("residue_select_rule");
    nametype_select_rule.name("nametype_select_rule");
    select_rule.name("select_rule");
    factor2.name("factor2");
    factor.name("factor");
    term.name("term");
    expr.name("expr");

    on_error<fail>(
            expr,
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
