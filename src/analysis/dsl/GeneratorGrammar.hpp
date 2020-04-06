//
// Created by xiamr on 7/8/19.
//

#ifndef TINKER_GENERATORGRAMMAR_HPP
#define TINKER_GENERATORGRAMMAR_HPP

#include <boost/config/warning_disable.hpp>
#include <boost/fusion/adapted/adt/adapt_adt.hpp>
#include <boost/fusion/include/adapt_adt.hpp>
#include <boost/fusion/include/std_pair.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>

#include "dsl/AmberMask.hpp"
#include "dsl/MacroRules.hpp"
#include "utils/common.hpp"

namespace karma = boost::spirit::karma;
namespace fusion = boost::fusion;
namespace phoenix = boost::phoenix;

BOOST_FUSION_ADAPT_ADT(std::shared_ptr<AmberMaskAST::residue_name_nums>,
                       (AmberMaskAST::select_ranges, AmberMaskAST::select_ranges, obj->val, /**/))

BOOST_FUSION_ADAPT_ADT(std::shared_ptr<AmberMaskAST::molecule_nums>,
                       (std::vector<AmberMaskAST::molecule_nums::Attr>, std::vector<AmberMaskAST::molecule_nums::Attr>,
                        obj->val,
                        /**/))

BOOST_FUSION_ADAPT_ADT(std::shared_ptr<AmberMaskAST::atom_name_nums>,
                       (AmberMaskAST::select_ranges, AmberMaskAST::select_ranges, obj->val, /**/))

BOOST_FUSION_ADAPT_ADT(std::shared_ptr<AmberMaskAST::atom_types>,
                       (AmberMaskAST::select_ranges, AmberMaskAST::select_ranges, obj->val, /**/))

BOOST_FUSION_ADAPT_ADT(std::shared_ptr<AmberMaskAST::atom_element_names>,
                       (std::vector<AmberMaskAST::Name>, const std::vector<AmberMaskAST::Name> &, obj->val, /**/))

BOOST_FUSION_ADAPT_ADT(std::shared_ptr<AmberMaskAST::Operator>, (AmberMaskAST::Op, AmberMaskAST::Op, obj->op,
                                                                 /**/)(AmberMask, AmberMask, obj->node1,
                                                                       /**/)(AmberMask, AmberMask, obj->node2, /**/))

template <typename Iterator> struct GeneratorGrammar : karma::grammar<Iterator, AmberMask()> {
    karma::rule<Iterator,
                boost::variant<fusion::vector<uint, boost::optional<std::pair<uint, int>>>, AmberMaskAST::Name>()>
        select_item_rule;
    karma::rule<Iterator, AmberMaskAST::Name()> name;
    karma::rule<Iterator, fusion::vector<uint, boost::optional<std::pair<uint, int>>>()> num_range;
    karma::rule<Iterator, int()> step_num;
    karma::rule<Iterator, std::shared_ptr<AmberMaskAST::residue_name_nums>()> residue_select_rule;
    karma::rule<Iterator, std::shared_ptr<AmberMaskAST::molecule_nums>()> molecule_select_rule;
    karma::rule<Iterator, std::shared_ptr<AmberMaskAST::atom_name_nums>()> atom_name_select_rule;
    karma::rule<Iterator, std::shared_ptr<AmberMaskAST::atom_types>()> atom_type_rule;
    karma::rule<Iterator, std::shared_ptr<AmberMaskAST::atom_element_names>()> atom_element_rule;
    karma::rule<Iterator, std::shared_ptr<AmberMaskAST::Operator>(int)> Operator;

    karma::rule<Iterator, AmberMask()> start;
    karma::rule<Iterator, AmberMask(int)> expr;

    GeneratorGrammar() : GeneratorGrammar::base_type(start) {
        using karma::_1;
        using karma::_a;
        using karma::_r1;
        using karma::_val;
        using karma::double_;
        using karma::eps;
        using karma::int_;
        using karma::lit;
        using karma::omit;
        using karma::string;
        using karma::uint_;
        using phoenix::at_c;

        residue_select_rule = ":" << (select_item_rule % ",")[_1 = at_c<0>(_val)];
        molecule_select_rule = '$' << ((uint_ << -('-' << uint_ << -step_num)) % ',')[_1 = at_c<0>(_val)];
        atom_name_select_rule = "@" << (select_item_rule % ",")[_1 = at_c<0>(_val)];
        atom_type_rule = "@%" << (select_item_rule % ",")[_1 = at_c<0>(_val)];
        atom_element_rule = "@/" << (name % ",")[_1 = at_c<0>(_val)];

        Operator = eps(at_c<0>(_val) == AmberMaskAST::Op::NOT) << "!" << expr(3)[_1 = at_c<1>(_val)] |
                   eps(at_c<0>(_val) == AmberMaskAST::Op::AND)
                       << (eps(_r1 > 2) << "(" | eps) << expr(2)[_1 = at_c<1>(_val)] << "&"
                       << expr(2)[_1 = at_c<2>(_val)] << (eps(_r1 > 2) << ")" | eps) |
                   eps(at_c<0>(_val) == AmberMaskAST::Op::OR)
                       << (eps(_r1 > 1) << "(" | eps) << expr(1)[_1 = at_c<1>(_val)] << "|"
                       << expr(1)[_1 = at_c<2>(_val)] << (eps(_r1 > 1) << ")" | eps);

        expr = eps(_val == AMBERMASK::System) << "System" | eps(_val == AMBERMASK::Protein) << "Protein" |
               eps(_val == AMBERMASK::Protein_H) << "Protein-H" | eps(_val == AMBERMASK::Backbone) << "Backbone" |
               eps(_val == AMBERMASK::MainChain) << "MainChain" |
               eps(_val == AMBERMASK::MainChain_plus_Cb) << "MainChain+Cb" |
               eps(_val == AMBERMASK::MainChain_plus_H) << "MainChain+H" |
               eps(_val == AMBERMASK::C_alpha) << "C-alpha" | eps(_val == AMBERMASK::SideChain) << "SideChain" |
               eps(_val == AMBERMASK::SideChain_H) << "SideChain-H" | eps(_val == AMBERMASK::DNA) << "DNA" |
               eps(_val == AMBERMASK::DNA_H) << "DNA-H" | eps(_val == AMBERMASK::RNA) << "RNA" |
               eps(_val == AMBERMASK::RNA_H) << "RNA-H" | eps(_val == AMBERMASK::Water) << "Water" | atom_element_rule |
               Operator(_r1) | atom_name_select_rule | atom_type_rule | residue_select_rule | molecule_select_rule;

        start = expr(0)[_1 = _val];

        step_num = eps(_val != 1) << "#" << int_[_1 = _val];
        num_range = uint_ << -("-" << uint_ << step_num);
        select_item_rule = num_range | name;
        name = string[_1 = (&_val)->*&AmberMaskAST::Name::name];
    }
};

inline bool format_node(const AmberMask &node, std::string &generated) {
    std::back_insert_iterator<std::string> sink{generated};
    GeneratorGrammar<std::back_insert_iterator<std::string>> serializer;
    return karma::generate(sink, serializer, node);
}

#endif // TINKER_GENERATORGRAMMAR_HPP
