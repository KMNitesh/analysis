#ifndef TINKER_AMBERMASKLIST_HPP
#define TINKER_AMBERMASKLIST_HPP

#include "dsl/grammar.hpp"

template <typename Iterator, typename Skipper>
struct AmberMaskListGrammar : qi::grammar<Iterator, std::vector<std::vector<AmberMask>>(), Skipper> {
    Grammar<Iterator, Skipper> maskParser;

    qi::rule<Iterator, std::vector<AmberMask>(), Skipper> plane_parser;

    qi::rule<Iterator, std::vector<std::vector<AmberMask>>(), Skipper> root;

    AmberMaskListGrammar() : AmberMaskListGrammar::base_type{root, "mask_list"} {
        plane_parser = '[' >> maskParser >> ';' >> maskParser >> ';' >> maskParser >> ']';
        root = +plane_parser;
    }
};

#endif  // TINKER_AMBERMASKLIST_HPP