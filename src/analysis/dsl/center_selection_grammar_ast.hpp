//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_CENTER_SELECTION_GRAMMAR_AST_HPP
#define TINKER_CENTER_SELECTION_GRAMMAR_AST_HPP

#include <boost/optional.hpp>
#include <memory>

#include "data_structure/atom.hpp"

struct MassCenterRuleNode {
    Atom::Node SelectionMask;

    MassCenterRuleNode(const Atom::Node &selectionMask) : SelectionMask(selectionMask) {}
};

struct GeomCenterRuleNode {
    Atom::Node SelectionMask;

    GeomCenterRuleNode(const Atom::Node &selectionMask) : SelectionMask(selectionMask) {}
};

struct NoopRuleNode {
    Atom::Node SelectionMask;

    NoopRuleNode(const Atom::Node &selectionMask) : SelectionMask(selectionMask) {}
};

struct EDARuleNode {
    Atom::Node mask1, mask2;
    EDARuleNode(const Atom::Node &mask1, const Atom::Node &mask2) : mask1(mask1), mask2(mask2) {}
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

struct HelpRuleNode {};

using CenterRuleNode =
    boost::variant<std::shared_ptr<MassCenterRuleNode>, std::shared_ptr<GeomCenterRuleNode>,
                   std::shared_ptr<NoopRuleNode>, std::shared_ptr<EDARuleNode>, std::shared_ptr<BondRuleNode>,
                   std::shared_ptr<AngleRuleNode>, std::shared_ptr<DihedralRuleNode>,
                   std::shared_ptr<BondedEnergyRuleNode>, std::shared_ptr<QuitRuleNode>, std::shared_ptr<HelpRuleNode>>;

#endif // TINKER_CENTER_SELECTION_GRAMMAR_AST_HPP
