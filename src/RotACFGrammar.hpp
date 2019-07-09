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
    vectorSelectorNode vectorSelctor;
    int Legendre;
    double time_increment_ps;
    double max_time_grap_ps;

    RotAcfNodeStruct(const vectorSelectorNode &vectorSelctor,
                     int legendre, double timeIncrementPs, double maxTimeGrapPs)
            : vectorSelctor(vectorSelctor), Legendre(legendre),
              time_increment_ps(timeIncrementPs), max_time_grap_ps(maxTimeGrapPs) {}

};

using RotAcfNode = std::shared_ptr<RotAcfNodeStruct>;

template<typename Iterator, typename Skipper>
struct RotAcfGrammar : qi::grammar<Iterator, RotAcfNode(), Skipper> {

    qi::rule<Iterator, DipoleVectorSelectorNode(), Skipper> RotAcfRule;
    VectorGrammar<Iterator, Skipper> vectorRule;

    RotAcfGrammar() : RotAcfGrammar::base_type(RotAcfRule) {
        using qi::lit;
        using qi::_val;
        using qi::_1;
        using qi::_2;
        using qi::_3;
        using qi::_4;
        using qi::uint_;
        using qi::double_;

        RotAcfRule =
                (lit("raf") >> "(" >> vectorRule >> "," >> uint_ >> "," >> double_ >> "," >> double_ >> ")")
                [_val = make_shared_<RotAcfNode>(_1, _2, _3, _4)];
    }
};


#endif //TINKER_ROTACFGRAMMAR_HPP
