//
// Created by xiamr on 7/8/19.
//

#ifndef TINKER_VECTORSELECTORFACTORYGRAMMAR_HPP
#define TINKER_VECTORSELECTORFACTORYGRAMMAR_HPP

#include "NormalVectorSelectorGrammar.hpp"
#include "TwoAtomSelectorGrammar.hpp"
#include "DipoleVectorSelectorGrammar.hpp"


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

using vectorSelectorNode = boost::variant<NormalVectorSelectorNode, TwoAtomVectorSelectorNode, DipoleVectorSelectorNode>;

template<typename Iterator, typename Skipper>
struct VectorGrammar : qi::grammar<Iterator, vectorSelectorNode(), Skipper> {

    qi::rule<Iterator, vectorSelectorNode(), Skipper> vector_selector_rule;
    NormalVectorGrammar<Iterator, Skipper> normal_vector_selector_rule;
    TwoAtomVectorGrammar<Iterator, Skipper> two_atom_vector_selector_rule;
    DipoleVectorGrammar<Iterator, Skipper> dipole_vector_selector_rule;

    VectorGrammar() : VectorGrammar::base_type(vector_selector_rule) {
        vector_selector_rule %=
                normal_vector_selector_rule | two_atom_vector_selector_rule | dipole_vector_selector_rule;
    }
};


#endif //TINKER_VECTORSELECTORFACTORYGRAMMAR_HPP
