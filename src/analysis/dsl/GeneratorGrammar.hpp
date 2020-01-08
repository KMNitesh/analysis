//
// Created by xiamr on 7/8/19.
//

#ifndef TINKER_GENERATORGRAMMAR_HPP
#define TINKER_GENERATORGRAMMAR_HPP

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/fusion/include/std_pair.hpp>
#include <boost/fusion/adapted/adt/adapt_adt.hpp>
#include <boost/fusion/include/adapt_adt.hpp>


#include "utils/common.hpp"
#include "data_structure/atom.hpp"


namespace karma = boost::spirit::karma;
namespace fusion = boost::fusion;
namespace phoenix = boost::phoenix;


BOOST_FUSION_ADAPT_ADT(
        std::shared_ptr<Atom::residue_name_nums>,
        (Atom::select_ranges, Atom::select_ranges, obj->val, /**/)
)

BOOST_FUSION_ADAPT_ADT(
        std::shared_ptr<Atom::molecule_nums>,
        (std::vector<Atom::molecule_nums::Attr>, std::vector<Atom::molecule_nums::Attr>, obj->val, /**/)
)

BOOST_FUSION_ADAPT_ADT(
        std::shared_ptr<Atom::atom_name_nums>,
        (Atom::select_ranges, Atom::select_ranges, obj->val, /**/)
)

BOOST_FUSION_ADAPT_ADT(
        std::shared_ptr<Atom::atom_types>,
        (Atom::select_ranges, Atom::select_ranges, obj->val, /**/)
)

BOOST_FUSION_ADAPT_ADT(
        std::shared_ptr<Atom::atom_element_names>,
        (std::vector<std::string>,
                const std::vector<std::string> &, obj->val, /**/)
)

BOOST_FUSION_ADAPT_ADT(
        std::shared_ptr<Atom::Operator>,
        (Atom::Op, Atom::Op, obj->op, /**/)
                (Atom::Node, Atom::Node, obj->node1, /**/)
                (Atom::Node, Atom::Node, obj->node2, /**/)
)


template<typename Iterator>
struct GeneratorGrammar : karma::grammar<Iterator, Atom::Node()> {
    karma::rule<Iterator, boost::variant<fusion::vector<uint, boost::optional<std::pair<uint, int>>>,
            std::string>()> select_item_rule;
    karma::rule<Iterator, fusion::vector<uint, boost::optional<std::pair<uint, int>>>()> num_range;
    karma::rule<Iterator, int()> step_num;
    karma::rule<Iterator, std::shared_ptr<Atom::residue_name_nums>()> residue_select_rule;
    karma::rule<Iterator, std::shared_ptr<Atom::molecule_nums>()> molecule_select_rule;
    karma::rule<Iterator, std::shared_ptr<Atom::atom_name_nums>()> atom_name_select_rule;
    karma::rule<Iterator, std::shared_ptr<Atom::atom_types>()> atom_type_rule;
    karma::rule<Iterator, std::shared_ptr<Atom::atom_element_names>()> atom_element_rule;
    karma::rule<Iterator, std::shared_ptr<Atom::Operator>(int)> Operator;

    karma::rule<Iterator, Atom::Node()> start;
    karma::rule<Iterator, Atom::Node(int)> expr;


    GeneratorGrammar() : GeneratorGrammar::base_type(start) {
        using karma::int_;
        using karma::double_;
        using karma::omit;
        using karma::eps;
        using karma::_val;
        using karma::_a;
        using karma::_1;
        using karma::_r1;
        using karma::string;
        using karma::lit;
        using karma::uint_;
        using phoenix::at_c;


        residue_select_rule = ":" << (select_item_rule % ",")[_1 = at_c<0>(_val)];
        molecule_select_rule = '$' << ((uint_ << -('-' << uint_ << -('#' << int_))) % ',')[_1 = at_c<0>(_val)];
        atom_name_select_rule = "@" << (select_item_rule % ",")[_1 = at_c<0>(_val)];
        atom_type_rule = "@%" << (select_item_rule % ",")[_1 = at_c<0>(_val)];
        atom_element_rule = "@/" << (string % ",")[_1 = at_c<0>(_val)];


        Operator = eps(at_c<0>(_val) == Atom::Op::NOT) << "!" << expr(3)[_1 = at_c<1>(_val)] |
                   eps(at_c<0>(_val) == Atom::Op::AND)
                           << (eps(_r1 > 2) << "(" | eps) << expr(2)[_1 = at_c<1>(_val)]
                           << "&" << expr(2)[_1 = at_c<2>(_val)] << (eps(_r1 > 2) << ")" | eps) |
                   eps(at_c<0>(_val) == Atom::Op::OR)
                           << (eps(_r1 > 1) << "(" | eps) << expr(1)[_1 = at_c<1>(_val)]
                           << "|" << expr(1)[_1 = at_c<2>(_val)] << (eps(_r1 > 1) << ")" | eps);

        expr = atom_element_rule | Operator(_r1) | atom_name_select_rule
               | atom_type_rule | residue_select_rule | molecule_select_rule;

        start = expr(0)[_1 = _val];

        step_num = eps(_val != 1) << "#" << int_[_1 = _val];
        num_range = uint_ << -("-" << uint_ << step_num);
        select_item_rule = num_range | string;
    }
};


inline bool format_node(const Atom::Node &node, std::string &generated) {
    std::back_insert_iterator<std::string> sink{generated};
    GeneratorGrammar<std::back_insert_iterator<std::string>> serializer;
    return karma::generate(sink, serializer, node);
}

#endif //TINKER_GENERATORGRAMMAR_HPP
