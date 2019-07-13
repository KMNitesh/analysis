//
// Created by xiamr on 7/8/19.
//

#ifndef TINKER_ROTACFGRAMMAR_HPP
#define TINKER_ROTACFGRAMMAR_HPP

#include <memory>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/phoenix.hpp>
#include <boost/variant.hpp>
#include <boost/optional.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/phoenix/function/adapt_function.hpp>
#include <boost/fusion/adapted/adt/adapt_adt.hpp>
#include <boost/fusion/include/adapt_adt.hpp>


#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>

#include "common.hpp"
#include "atom.hpp"
#include "VectorSelectorFactoryGrammar.hpp"

namespace qi = boost::spirit::qi;
namespace fusion = boost::fusion;
namespace phoenix = boost::phoenix;

struct RotAcfNodeStruct {
    boost::optional<vectorSelectorNode> vectorSelctor;
    boost::optional<int> Legendre;
    boost::optional<double> time_increment_ps;
    boost::optional<double> max_time_grap_ps;
};


using RotAcfNode = std::shared_ptr<RotAcfNodeStruct>;


BOOST_FUSION_ADAPT_ADT(
        std::shared_ptr<RotAcfNodeStruct>,
        (boost::optional<vectorSelectorNode>, boost::optional<vectorSelectorNode>, obj->vectorSelctor, obj->vectorSelctor = val)
                (boost::optional<int>, boost::optional<int>, obj->Legendre, obj->Legendre = val)
                (boost::optional<double>, boost::optional<double>, obj->time_increment_ps, obj->time_increment_ps = val)
                (boost::optional<double>, boost::optional<double>, obj->max_time_grap_ps, obj->max_time_grap_ps = val)

)


template<typename Iterator, typename Skipper>
struct RotAcfGrammar : qi::grammar<Iterator, RotAcfNode(), Skipper> {

    qi::rule<Iterator, RotAcfNode(), Skipper> RotAcfRule;
    qi::rule<Iterator, void(RotAcfNode), Skipper> option;
    VectorGrammar<Iterator, Skipper> vectorRule;

    RotAcfGrammar() : RotAcfGrammar::base_type(RotAcfRule) {
        using qi::lit;
        using qi::_val;
        using qi::_1;
        using qi::_2;
        using qi::_3;
        using qi::_4;
        using qi::int_;
        using qi::double_;
        using phoenix::at_c;
        using qi::_r1;

        option = lit("vector") >> "=" >> vectorRule[at_c<0>(_r1) = _1] |
                 lit("P") >> "=" >> int_[at_c<1>(_r1) = _1] |
                 lit("time_increment_ps") >> "=" >> double_[at_c<2>(_r1) = _1] |
                 lit("max_time_grap_ps") >> "=" >> double_[at_c<3>(_r1) = _1];

        RotAcfRule = "rotacf" >> lit("(")[_val = make_shared_<RotAcfNodeStruct>()] >> (option(_val) % ",") >> ")";
    }
};


inline bool operator==(const RotAcfNode &node1, const RotAcfNode &node2) {
    if (node1 and node2) {
        return node1->vectorSelctor == node2->vectorSelctor and node1->Legendre == node2->Legendre and (
                node1->time_increment_ps and node2->time_increment_ps ?
                std::abs(node1->time_increment_ps.value() - node2->time_increment_ps.value()) <
                std::numeric_limits<double>::epsilon() : false) and (
                       node1->max_time_grap_ps and node2->max_time_grap_ps ?
                       std::abs(node1->max_time_grap_ps.value() - node2->max_time_grap_ps.value()) <
                       std::numeric_limits<double>::epsilon() : false);
    }
    return false;
}


#endif //TINKER_ROTACFGRAMMAR_HPP
