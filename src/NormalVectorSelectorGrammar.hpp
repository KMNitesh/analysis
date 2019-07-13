//
// Created by xiamr on 7/8/19.
//

#ifndef TINKER_NORMALVECTORSELECTORGRAMMAR_HPP
#define TINKER_NORMALVECTORSELECTORGRAMMAR_HPP

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
#include "grammar.hpp"

namespace qi = boost::spirit::qi;
namespace fusion = boost::fusion;
namespace phoenix = boost::phoenix;

struct NormalVectorSelectorNodeStruct {
    Atom::Node id1;
    Atom::Node id2;
    Atom::Node id3;

    NormalVectorSelectorNodeStruct(const Atom::Node &id1,
                                   const Atom::Node &id2,
                                   const Atom::Node &id3)
            : id1(id1), id2(id2), id3(id3) {}
};

using NormalVectorSelectorNode = std::shared_ptr<NormalVectorSelectorNodeStruct>;

template<typename Iterator, typename Skipper>
struct NormalVectorGrammar : qi::grammar<Iterator, NormalVectorSelectorNode(), Skipper> {

    qi::rule<Iterator, NormalVectorSelectorNode(), Skipper> normal_vector_selector_rule;
    qi::rule<Iterator, Atom::Node(), Skipper> mask;
    Grammar<Iterator, Skipper> maskParser;

    NormalVectorGrammar() : NormalVectorGrammar::base_type(normal_vector_selector_rule) {
        using qi::lit;
        using qi::_val;
        using qi::_1;
        using qi::_2;
        using qi::_3;

        mask %= "[" >> maskParser >> "]";

        normal_vector_selector_rule =
                (lit("normalVector") >> "(" >> mask >> "," >> mask >> "," >> mask >> ")")
                [_val = make_shared_<NormalVectorSelectorNodeStruct>(_1, _2, _3)];
    }
};

inline bool operator==(const NormalVectorSelectorNode &node1, const NormalVectorSelectorNode &node2) {
    if (node1 && node2) {
        return node1->id1 == node2->id1 && node1->id2 == node2->id2 && node1->id3 == node2->id3;
    }
    return false;
}

#endif //TINKER_NORMALVECTORSELECTORGRAMMAR_HPP
