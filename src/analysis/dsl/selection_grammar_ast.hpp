//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_SELECTION_GRAMMAR_AST_HPP
#define TINKER_SELECTION_GRAMMAR_AST_HPP

#include <boost/any.hpp>
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

struct LogicalOperator {
    enum class Op { Less, LessEqual, Equal, Great, GreatEqual };
    std::string id;
    Op op;
    boost::any value;

    LogicalOperator(std::string id, Op op, boost::any value) : id(std::move(id)), op(op), value(std::move(value)) {}
};

struct LikeOperator {
    std::string id, pattern;

    LikeOperator(std::string id, std::string pattern) : id(std::move(id)), pattern(std::move(pattern)) {}
};

struct RegexOperator {
    std::string id, regex;

    RegexOperator(std::string id, std::string regex) : id(std::move(id)), regex(std::move(regex)) {}
};

struct BetweenOperator {
    std::string id;
    double low_bound, up_bound;

    BetweenOperator(std::string id, double low_bound, double up_bound)
        : id(std::move(id)), low_bound(low_bound), up_bound(up_bound) {}
};

struct InOperator {
    std::string id;
    std::vector<std::string> fs;

    InOperator(std::string id, std::vector<std::string> fs) : id(std::move(id)), fs(std::move(fs)) {}
};

struct CombineOperator;

using Condition = boost::variant<boost::blank, std::shared_ptr<CombineOperator>, std::shared_ptr<LogicalOperator>,
                                 std::shared_ptr<LikeOperator>, std::shared_ptr<RegexOperator>,
                                 std::shared_ptr<BetweenOperator>, std::shared_ptr<InOperator>>;

struct CombineOperator {
    enum class Op { AND, OR, NOT };

    Op op;
    Condition cond1, cond2;

    CombineOperator(Op op, Condition cond1, Condition cond2 = {})
        : op(op), cond1(std::move(cond1)), cond2(std::move(cond2)) {}
};

struct OrderBy {
    enum class Order { ASC, DESC };

    std::string field;
    Order order;

    OrderBy() = default;
    OrderBy(std::string field, Order order) : field(std::move(field)), order(order) {}
};

struct SelectStmt {
    bool DISTINC;
    std::vector<std::string> fields;
    AmberMask mask;
    boost::optional<Condition> where;
    boost::optional<OrderBy> orderby;
    boost::optional<uint> limit;

    SelectStmt() = default;
    SelectStmt(bool DISTINC, std::vector<std::string> fields, AmberMask mask, boost::optional<Condition> where,
               boost::optional<OrderBy> orderby, boost::optional<uint> limit)
        : DISTINC(DISTINC), fields(std::move(fields)), mask(std::move(mask)), where(std::move(where)),
          orderby(std::move(orderby)), limit(std::move(limit)) {}
};

struct QuitRuleNode {};

struct HelpRuleNode {
    boost::optional<std::string> keyword;
    HelpRuleNode(boost::optional<std::string> keyword) : keyword(std::move(keyword)) {}
};

using SelectionAST =
    boost::variant<boost::blank, MassCenterRuleNode, GeomCenterRuleNode, NoopRuleNode, EDARuleNode, BondRuleNode,
                   AngleRuleNode, DihedralRuleNode, BondedEnergyRuleNode, SelectStmt, QuitRuleNode, HelpRuleNode>;

#endif // TINKER_SELECTION_GRAMMAR_AST_HPP
