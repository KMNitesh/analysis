//
// Created by xiamr on 7/13/19.
//

#ifndef TINKER_ROTACFCUTOFFGRAMMAR_HPP
#define TINKER_ROTACFCUTOFFGRAMMAR_HPP

#include <memory>

#include <boost/fusion/adapted/adt/adapt_adt.hpp>
#include <boost/fusion/include/adapt_adt.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/optional.hpp>
#include <boost/phoenix.hpp>
#include <boost/phoenix/function/adapt_function.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/variant.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>

#include "VectorSelectorFactoryGrammar.hpp"
#include "atom.hpp"
#include "common.hpp"

namespace qi = boost::spirit::qi;
namespace fusion = boost::fusion;
namespace phoenix = boost::phoenix;

struct RotAcfCutoffNodeStruct {
    Atom::Node M;
    Atom::Node L;
    boost::optional<vectorSelectorNode> vectorSelctor;
    boost::optional<int> Legendre;
    boost::optional<double> cutoff;
    boost::optional<double> time_increment_ps;
    boost::optional<double> max_time_grap_ps;
    std::string outfilename;
};

using RotAcfCutoffNode = std::shared_ptr<RotAcfCutoffNodeStruct>;

BOOST_FUSION_ADAPT_ADT(
        RotAcfCutoffNode,
        (Atom::Node, Atom::Node, obj->M, obj->M = val)
                (Atom::Node, Atom::Node, obj->L, obj->L = val)
                (boost::optional<vectorSelectorNode>, boost::optional<vectorSelectorNode>, obj->vectorSelctor,
                 obj->vectorSelctor = val)
                (boost::optional<int>, boost::optional<int>, obj->Legendre, obj->Legendre = val)
                (boost::optional<double>, boost::optional<double>, obj->cutoff, obj->cutoff = val)
                (boost::optional<double>, boost::optional<double>, obj->time_increment_ps, obj->time_increment_ps = val)
                (boost::optional<double>, boost::optional<double>, obj->max_time_grap_ps, obj->max_time_grap_ps = val)
                (std::string, std::string, obj->outfilename, obj->outfilename = val)

)

template<typename Iterator, typename Skipper>
struct RotAcfCutoffGrammar : qi::grammar<Iterator, RotAcfCutoffNode(), Skipper> {

    qi::rule<Iterator, RotAcfCutoffNode(), Skipper> RotAcfCutoffRule;
    qi::rule<Iterator, void(RotAcfCutoffNode), Skipper> option;
    VectorGrammar<Iterator, Skipper> vectorRule;
    Grammar<Iterator, Skipper> maskParser;
    qi::rule<Iterator, Atom::Node(), Skipper> mask;

    RotAcfCutoffGrammar() : RotAcfCutoffGrammar::base_type(RotAcfCutoffRule) {
        using phoenix::at_c;
        using qi::_1;
        using qi::_r1;
        using qi::_val;
        using qi::double_;
        using qi::int_;
        using qi::lit;
        using qi::as_string;
        using qi::lexeme;
        using qi::alnum;
        using qi::char_;

        mask %= "[" >> maskParser >> "]";

        option = lit("M") >> "=" >> mask[at_c<0>(_r1) = _1] |
                 lit("L") >> "=" >> mask[at_c<1>(_r1) = _1] |
                 lit("vector") >> "=" >> vectorRule[at_c<2>(_r1) = _1] |
                 lit("P") >> "=" >> int_[at_c<3>(_r1) = _1] |
                 lit("cutoff") >> "=" >> double_[at_c<4>(_r1) = _1] |
                 lit("time_increment_ps") >> "=" >> double_[at_c<5>(_r1) = _1] |
                 lit("max_time_grap_ps") >> "=" >> double_[at_c<6>(_r1) = _1] |
                 lit("out") >> "=" >> as_string[lexeme[+(alnum | char_("_.-"))]][at_c<7>(_r1) = _1];

        RotAcfCutoffRule =
                "rotacfcutoff" >> lit("(")[_val = make_shared_<RotAcfCutoffNodeStruct>()] >> (option(_val) % ",")
                               >> ")";
    }
};

inline bool operator==(const RotAcfCutoffNode &node1, const RotAcfCutoffNode &node2) {
    if (node1 and node2) {
        return node1->M == node2->M and node1->L == node2->L and node1->vectorSelctor == node2->vectorSelctor and
               node1->Legendre == node2->Legendre and
               (node1->cutoff and node2->cutoff
                ? std::abs(node1->cutoff.value() - node2->cutoff.value()) < std::numeric_limits<double>::epsilon()
                : false) and
               (node1->time_increment_ps and node2->time_increment_ps
                ? std::abs(node1->time_increment_ps.value() - node2->time_increment_ps.value()) <
                  std::numeric_limits<double>::epsilon()
                : false) and
               (node1->max_time_grap_ps and node2->max_time_grap_ps
                ? std::abs(node1->max_time_grap_ps.value() - node2->max_time_grap_ps.value()) <
                  std::numeric_limits<double>::epsilon()
                : false);
    }
    return false;
}

#endif // TINKER_ROTACFCUTOFFGRAMMAR_HPP
