//
// Created by xiamr on 7/13/19.
//

#ifndef TINKER_LANGUAGEGRAMMAR_HPP
#define TINKER_LANGUAGEGRAMMAR_HPP

#include "GoGrammar.hpp"
#include "RotACFGrammar.hpp"

using LanguageNodeVariant = boost::variant<RotAcfNode, GoNode>;

using LanguageNode = std::shared_ptr<LanguageNodeVariant>;

using Language = std::vector<LanguageNode>;

template<typename Iterator, typename Skipper>
struct LanguageGrammar : qi::grammar<Iterator, Language(), Skipper> {

    qi::rule<Iterator, Language(), Skipper> languageRule;
    qi::rule<Iterator, LanguageNode(), Skipper> stmtRule;

    qi::rule<Iterator, Skipper> semicolons;

    GoGrammar<Iterator, Skipper> goRule;
    RotAcfGrammar<Iterator, Skipper> rotAcfRule;

    LanguageGrammar() : LanguageGrammar::base_type(languageRule) {
        using qi::lit;
        using qi::_val;
        using qi::_1;
        using phoenix::push_back;

        stmtRule = rotAcfRule[_val = make_shared_<LanguageNodeVariant>(_1)];


        semicolons = +(lit(";"));

        languageRule = (stmtRule[push_back(_val, _1)] % semicolons)
                >> -(semicolons >> goRule[push_back(_val, make_shared_<LanguageNodeVariant>(_1))]) >> -(semicolons);
    }
};

inline bool operator==(const LanguageNode &node1, const LanguageNode &node2) {
    if (node1 and node2) {
        return *node1 == *node2;
    }
    return false;
}


#endif //TINKER_LANGUAGEGRAMMAR_HPP
