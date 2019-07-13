//
// Created by xiamr on 7/13/19.
//

#ifndef TINKER_GOGRAMMAR_HPP
#define TINKER_GOGRAMMAR_HPP

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

struct GoNodeStruct {
    unsigned int start = 1;
    unsigned int end = 0;
    unsigned int step = 1;
    unsigned int nthreads = 0;
};


using GoNode = std::shared_ptr<GoNodeStruct>;


BOOST_FUSION_ADAPT_ADT(
        GoNode,
        (unsigned int, unsigned int, obj->start, obj->start = val)
                (unsigned int, unsigned int, obj->end, obj->end = val)
                (unsigned int, unsigned int, obj->step, obj->step = val)
                (unsigned int, unsigned int, obj->nthreads, obj->nthreads = val)
)


template<typename Iterator, typename Skipper>
struct GoGrammar : qi::grammar<Iterator, GoNode(), Skipper> {

    qi::rule<Iterator, GoNode(), Skipper> GoRule;
    qi::rule<Iterator, void(GoNode), Skipper> option;

    GoGrammar() : GoGrammar::base_type(GoRule) {
        using qi::lit;
        using qi::_val;
        using qi::_1;
        using qi::uint_;
        using phoenix::at_c;
        using qi::_r1;

        option = lit("start") >> "=" >> uint_[at_c<0>(_r1) = _1] |
                 lit("end") >> "=" >> uint_[at_c<1>(_r1) = _1] |
                 lit("step") >> "=" >> uint_[at_c<2>(_r1) = _1] |
                 lit("nthreads") >> "=" >> uint_[at_c<3>(_r1) = _1];
        GoRule = "go" >> lit("(")[_val = make_shared_<GoNodeStruct>()] >> -(option(_val) % ",") >> ")";
    }
};


inline bool operator==(const GoNode &node1, const GoNode &node2) {
    if (node1 and node2) {
        return node1->start == node2->start and node1->end == node2->end and
               node1->step == node2->step and node1->nthreads == node2->nthreads;

    }
    return false;
}


#endif //TINKER_GOGRAMMAR_HPP
