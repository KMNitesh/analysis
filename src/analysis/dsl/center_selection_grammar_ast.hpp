//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_CENTER_SELECTION_GRAMMAR_AST_HPP
#define TINKER_CENTER_SELECTION_GRAMMAR_AST_HPP

#include <boost/blank.hpp>
#include <boost/optional.hpp>
#include <memory>

#include "dsl/AmberMask.hpp"

struct MassCenterRuleNode {
    AmberMask SelectionMask;

    MassCenterRuleNode(AmberMask mask) : SelectionMask(std::move(mask)) {}
};

struct GeomCenterRuleNode {
    AmberMask SelectionMask;

    GeomCenterRuleNode(AmberMask mask) : SelectionMask(std::move(mask)) {}
};

struct NoopRuleNode {
    AmberMask SelectionMask;

    NoopRuleNode(AmberMask mask) : SelectionMask(std::move(mask)) {}
};

struct EDARuleNode {
    AmberMask mask1, mask2;
    EDARuleNode(AmberMask mask1, AmberMask mask2) : mask1(std::move(mask1)), mask2(std::move(mask2)) {}
};

struct BondRuleNode {
    AmberMask mask1, mask2;
    BondRuleNode(AmberMask mask1, AmberMask mask2) : mask1(std::move(mask1)), mask2(std::move(mask2)) {}
};

struct AngleRuleNode {
    AmberMask mask1, mask2, mask3;
    AngleRuleNode(AmberMask mask1, AmberMask mask2, AmberMask mask3)
        : mask1(std::move(mask1)), mask2(std::move(mask2)), mask3(std::move(mask3)) {}
};

struct DihedralRuleNode {
    AmberMask mask1, mask2, mask3, mask4;
    DihedralRuleNode(AmberMask mask1, AmberMask mask2, AmberMask mask3, AmberMask mask4)
        : mask1(std::move(mask1)), mask2(std::move(mask2)), mask3(std::move(mask3)), mask4(std::move(mask4)) {}
};

struct BondedEnergyRuleNode {
    AmberMask mask;
    BondedEnergyRuleNode(AmberMask mask) : mask(std::move(mask)) {}
};

struct QuitRuleNode {};

struct HelpRuleNode {
    boost::optional<std::string> keyword;
    HelpRuleNode(boost::optional<std::string> keyword) : keyword(std::move(keyword)) {}
};

using CenterRuleNode =
    boost::variant<boost::blank, MassCenterRuleNode, GeomCenterRuleNode, NoopRuleNode, EDARuleNode, BondRuleNode,
                   AngleRuleNode, DihedralRuleNode, BondedEnergyRuleNode, QuitRuleNode, HelpRuleNode>;

#endif // TINKER_CENTER_SELECTION_GRAMMAR_AST_HPP
