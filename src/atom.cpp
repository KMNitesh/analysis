//
// Created by xiamr on 3/17/19.
//


#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/phoenix.hpp>
#include <boost/variant.hpp>
#include <boost/optional.hpp>
#include <boost/fusion/sequence/intrinsic/at_c.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/phoenix/function/adapt_function.hpp>

#include <fnmatch.h>

#include "common.hpp"
#include "atom.hpp"
#include "grammar.hpp"

namespace qi = boost::spirit::qi;
namespace fusion = boost::fusion;
namespace phoenix = boost::phoenix;



void print::operator()(const std::shared_ptr<Atom::residue_name_nums> &residues) const {
    indent(space_num);
    bool first = true;
    struct print_res : boost::static_visitor<> {
        void operator()(fusion::vector<uint, boost::optional<uint>> &i) {
            std::cout << fusion::at_c<0>(i);
            auto op = fusion::at_c<1>(i);
            if (op) std::cout << "-" << op.get();
        }

        void operator()(std::string &i) {
            std::cout << i;
        }
    } p;
    for (auto &i : residues->val) {
        if (first) {
            std::cout << "residues : ";
            first = false;
        } else {
            std::cout << ",";
        }
        boost::apply_visitor(p, i);
    }
    std::cout << std::endl;
}

void print::operator()(const std::shared_ptr<Atom::atom_name_nums> &names) const {
    struct print_res : boost::static_visitor<> {
        void operator()(const fusion::vector<uint, boost::optional<uint>> &i) {
            std::cout << fusion::at_c<0>(i);
            auto op = fusion::at_c<1>(i);
            if (op) std::cout << "-" << op.get();
        }

        void operator()(const std::string &i) {
            std::cout << i;
        }
    } p;
    indent(space_num);
    bool first = true;
    for (auto &i : names->val) {
        if (first) {
            std::cout << "names : ";
            first = false;
        } else {
            std::cout << ",";
        }
        boost::apply_visitor(p, i);
    }
    std::cout << std::endl;
}

void print::operator()(const std::shared_ptr<Atom::atom_types> &types) const {
    indent(space_num);
    struct print_types : boost::static_visitor<> {
        void operator()(const fusion::vector<uint, boost::optional<uint>> &i) {
            std::cout << fusion::at_c<0>(i);
            auto op = fusion::at_c<1>(i);
            if (op) std::cout << "-" << op.get();
        }

        void operator()(const std::string &i) {
            std::cout << i;
        }
    } p;
    bool first = true;
    for (auto &i : types->val) {
        if (first) {
            std::cout << "types : ";
            first = false;
        } else {
            std::cout << ",";
        }
        boost::apply_visitor(p, i);
    }
    std::cout << std::endl;
}

void print::operator()(const std::shared_ptr<Atom::atom_element_names> &ele) const {
    indent(space_num);
    bool first = true;
    for (auto &i : ele->val) {
        if (first) {
            std::cout << "elements : ";
            first = false;
        } else {
            std::cout << ",";
        }
        std::cout << i;
    }
    std::cout << std::endl;
}

void print::operator()(const std::shared_ptr<Atom::Operator> &op) const {
    switch (op->op) {
        case Atom::Op::NOT:
            indent(space_num);
            std::cout << "!" << std::endl;
            boost::apply_visitor(print(space_num + 1), op->node1);
            break;
        case Atom::Op::AND:
            indent(space_num);
            std::cout << "&" << std::endl;
            boost::apply_visitor(print(space_num + 1), op->node1);
            boost::apply_visitor(print(space_num + 1), op->node2);
            break;
        case Atom::Op::OR:
            indent(space_num);
            std::cout << "|" << std::endl;
            boost::apply_visitor(print(space_num + 1), op->node1);
            boost::apply_visitor(print(space_num + 1), op->node2);
            break;
    }
}



template<typename Iterator, typename Skipper>
Atom::AtomIndenter input_atom_selection(const Grammar<Iterator,Skipper> &grammar, const std::string &promot) {

    for (;;) {
        Atom::AtomIndenter mask;
        mask.input_string = input(promot);
        boost::trim(mask.input_string);
        auto it = mask.input_string.begin();

        bool status = qi::phrase_parse(it, mask.input_string.end(), grammar, qi::ascii::space, mask.ast);

        if (status) {
            std::cout << "Parsed Abstract Syntax Tree :" << std::endl;
            boost::apply_visitor(print(), mask.ast);
        }

        if (!(status and (it == mask.input_string.end()))) {
            std::cout << "error-pos : " << std::endl;
            std::cout << mask.input_string << std::endl;
            for (auto iter = mask.input_string.begin(); iter != it; ++iter) std::cout << " ";
            std::cout << "^" << std::endl;

            continue;
        }
        return mask;

    }
}


void
Atom::select2group(Atom::AtomIndenter &ids1, Atom::AtomIndenter &ids2,
                   const std::string &prompt1, const std::string &prompt2) {

    Grammar<std::string::iterator, qi::ascii::space_type> grammar;

    ids1 = input_atom_selection(grammar, prompt1);
    ids2 = input_atom_selection(grammar, prompt2);

}

void Atom::select1group(AtomIndenter &ids, const std::string &prompt) {
    namespace qi = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;
    using ascii::char_;

    Grammar<std::string::iterator, qi::ascii::space_type> grammar;

    ids = input_atom_selection(grammar, prompt);
}


bool Atom::is_match(const std::shared_ptr<Atom> &atom, const Atom::AtomIndenter &id) {
    struct Equal : boost::static_visitor<bool> {

        explicit Equal(const std::shared_ptr<Atom> &atom) : atom(atom) {}


        bool operator()(const std::shared_ptr<residue_name_nums> &residues) const {
            if (!atom->residue_name or !atom->residue_num) {
                std::cerr << "residue selection syntax is invaild in current context" << std::endl;
                exit(1);
            }
            struct Equal_residue : boost::static_visitor<bool> {
                explicit Equal_residue(const std::shared_ptr<Atom> &atom) : atom(atom) {}

                bool operator()(const fusion::vector<uint, boost::optional<uint>> &i) {
                    auto op = fusion::at_c<1>(i);
                    if (op) {
                        if (op.get() >= fusion::at_c<0>(i))
                            return atom->residue_num >= fusion::at_c<0>(i) and atom->residue_num <= op.get();
                        else
                            return atom->residue_num >= op.get() and atom->residue_num <= fusion::at_c<0>(i);
                    } else {
                        return atom->residue_num == fusion::at_c<0>(i);
                    }
                }

                bool operator()(const std::string &pattern) {
                    if (fnmatch(pattern.c_str(), atom->residue_name.get().c_str(), FNM_CASEFOLD) == 0) return true;
                    std::string num_str = boost::lexical_cast<std::string>(atom->residue_num.get());
                    return fnmatch(pattern.c_str(), num_str.c_str(), FNM_CASEFOLD) == 0;
                }

            private:
                const std::shared_ptr<Atom> &atom;
            } equal(atom);

            for (auto &i : residues->val) {
                if (boost::apply_visitor(equal, i)) return true;
            }
            return false;
        }

        bool operator()(const std::shared_ptr<atom_name_nums> &names) const {
            struct Equal_atom : boost::static_visitor<bool> {
                explicit Equal_atom(const std::shared_ptr<Atom> &atom) : atom(atom) {}

                bool operator()(const fusion::vector<uint, boost::optional<uint>> &i) {
                    auto op = fusion::at_c<1>(i);
                    if (op) {
                        if (op.get() >= fusion::at_c<0>(i))
                            return atom->seq >= fusion::at_c<0>(i) and atom->seq <= op.get();
                        else
                            return atom->seq >= op.get() and atom->seq <= fusion::at_c<0>(i);
                    } else {
                        return atom->seq == fusion::at_c<0>(i);
                    }
                }

                bool operator()(const std::string &pattern) {
                    if (fnmatch(pattern.c_str(), atom->atom_name.c_str(), FNM_CASEFOLD) == 0) return true;
                    std::string num_str = boost::lexical_cast<std::string>(atom->seq);
                    return fnmatch(pattern.c_str(), num_str.c_str(), FNM_CASEFOLD) == 0;
                }

            private:
                const std::shared_ptr<Atom> &atom;
            } equal(atom);

            for (auto &i : names->val) {
                if (boost::apply_visitor(equal, i)) return true;
            }
            return false;
        }

        bool operator()(const std::shared_ptr<atom_types> types) const {
            struct Equal_types : boost::static_visitor<bool> {
                explicit Equal_types(const std::shared_ptr<Atom> &atom) : atom(atom) {}

                bool operator()(const fusion::vector<uint, boost::optional<uint>> &i) {
                    auto op = fusion::at_c<1>(i);
                    if (op) {
                        if (op.get() >= fusion::at_c<0>(i))
                            return atom->typ >= static_cast<int>(fusion::at_c<0>(i)) and atom->typ <= static_cast<int>(op.get());
                        else
                            return atom->typ >= static_cast<int>(op.get()) and atom->typ <= static_cast<int>(fusion::at_c<0>(i));
                    } else {
                        return atom->typ == static_cast<int>(fusion::at_c<0>(i));
                    }
                }

                bool operator()(const std::string &pattern) {
                    return fnmatch(pattern.c_str(), atom->type_name.c_str(), FNM_CASEFOLD) == 0;
                    std::string num_str = boost::lexical_cast<std::string>(atom->typ);
                    return fnmatch(pattern.c_str(), num_str.c_str(), FNM_CASEFOLD) == 0;
                }

            private:
                const std::shared_ptr<Atom> &atom;
            } equal(atom);

            for (auto &i : types->val) {
                if (boost::apply_visitor(equal, i)) return true;
            }
            return false;
        }

        bool operator()(const std::shared_ptr<atom_element_names> &ele) const {
            if (!atom->atom_symbol) {
                std::cerr << "atom element symbol selection syntax is invaild in current context" << std::endl;
                exit(1);
            }
            for (auto &pattern : ele->val) {
                if (fnmatch(pattern.c_str(), atom->atom_symbol.get().c_str(), FNM_CASEFOLD) == 0) return true;
            }
            return false;

        }

        bool operator()(const std::shared_ptr<Operator> &op) const {
            switch (op->op) {
                case Op::NOT:
                    return not boost::apply_visitor(Equal(atom), op->node1);
                    break;
                case Op::AND: {
                    Equal equal(atom);
                    return boost::apply_visitor(equal, op->node1) and boost::apply_visitor(equal, op->node2);
                }
                    break;
                case Op::OR: {
                    Equal equal(atom);
                    return boost::apply_visitor(equal, op->node1) or boost::apply_visitor(equal, op->node2);
                }
                    break;
                default:
                    std::abort();
            }
        }

    private:
        const std::shared_ptr<Atom> &atom;
    };


    return boost::apply_visitor(Equal(atom), id.ast);
}
