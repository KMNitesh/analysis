//
// Created by xiamr on 3/17/19.
//

#include "atom.hpp"

#include <fnmatch.h>

#include <boost/fusion/include/at_c.hpp>
#include <boost/optional.hpp>
#include <boost/spirit/include/support_line_pos_iterator.hpp>
#include <boost/variant.hpp>

#include "dsl/GeneratorGrammar.hpp"
#include "dsl/grammar.hpp"
#include "molecule.hpp"
#include "utils/ThrowAssert.hpp"
#include "utils/common.hpp"

namespace qi = boost::spirit::qi;
namespace fusion = boost::fusion;
namespace phoenix = boost::phoenix;

std::ostream &operator<<(std::ostream &os, const AmberMask &mask) {
    os << to_string(mask);
    return os;
}

void print::operator()(const std::shared_ptr<Atom::residue_name_nums> &residues) const {
    indent(space_num);
    bool first = true;
    struct print_res : boost::static_visitor<> {
        void operator()(const fusion::vector<uint, boost::optional<std::pair<uint, int>>> &i) {
            std::cout << fusion::at_c<0>(i);
            auto op = fusion::at_c<1>(i);
            if (op) {
                std::cout << "-" << op.get().first;
                if (op.get().second != 1) {
                    std::cout << '#' << op.get().second;
                }
            }
        }

        void operator()(Atom::Name &i) { std::cout << i; }
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

void print::operator()(const std::shared_ptr<Atom::molecule_nums> &molecule) const {
    indent(space_num);
    bool first = true;
    for (auto &i : molecule->val) {
        if (first) {
            std::cout << "molecule num : ";
            first = false;
        } else {
            std::cout << ",";
        }
        std::cout << fusion::at_c<0>(i);
        auto op = fusion::at_c<1>(i);
        if (op) {
            std::cout << "-" << fusion::at_c<0>(op.get());
            auto &step = fusion::at_c<1>(op.get());
            if (step and step.get() != 1) {
                std::cout << '#' << step.get();
            }
        }
    }
    std::cout << std::endl;
}

void print::operator()(const std::shared_ptr<Atom::atom_name_nums> &names) const {
    struct print_res : boost::static_visitor<> {
        void operator()(const fusion::vector<uint, boost::optional<std::pair<uint, int>>> &i) {
            std::cout << fusion::at_c<0>(i);
            auto op = fusion::at_c<1>(i);
            if (op) {
                std::cout << "-" << op.get().first;
                if (op.get().second != 1) {
                    std::cout << '#' << op.get().second;
                }
            }
        }

        void operator()(const std::string &i) { std::cout << i; }
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
        void operator()(const fusion::vector<uint, boost::optional<std::pair<uint, int>>> &i) {
            std::cout << fusion::at_c<0>(i);
            auto op = fusion::at_c<1>(i);
            if (op) {
                std::cout << "-" << op.get().first;
                if (op.get().second != 1) {
                    std::cout << '#' << op.get().second;
                }
            }
        }

        void operator()(const std::string &i) { std::cout << i; }
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
    if (op) {
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
}

void print::indent(int space_num) const { std::cout << std::string(3 * space_num, ' '); }

template <typename Iterator, typename Skipper>
Atom::AmberMask input_atom_selection(const Grammar<Iterator, Skipper> &grammar, const std::string &prompt,
                                     std::string input_string = "", bool allow_empty = false, bool quiet = false) {
    for (;;) {
        Atom::AmberMask mask;
        if (input_string.empty()) {
            if (isatty(STDIN_FILENO) == 1) {
                char *buf = readline::readline(prompt.c_str());
                input_string = buf;
                free(buf);
            } else {
                input_string = input(prompt);
            }
        }

        boost::trim(input_string);
        if (input_string.empty()) {
            if (allow_empty)
                return mask;
            else
                continue;
        }

        std::string::iterator begin{std::begin(input_string)}, it{begin}, end{std::end(input_string)};
        try {
            qi::phrase_parse(it, end, grammar > qi::eoi, boost::spirit::ascii::space, mask);
            if (!quiet) {
                std::cout << "Parsed Abstract Syntax Tree :\n";
                boost::apply_visitor(print(), mask);
            }
            return mask;
        } catch (const qi::expectation_failure<std::string::iterator> &x) {
            std::cerr << "Grammar Parse Failure ! Expecting : " << x.what_ << '\n';
            auto column = boost::spirit::get_column(begin, x.first);
            std::string pos = " (column: " + std::to_string(column) + ")";
            std::cerr << pos << ">>>>" << input_string << "<<<<\n";
            std::cerr << std::string(column + pos.size() + 3, ' ') << "^~~~ here\n";
        }
        input_string.clear();
    }
}

AmberMask parse_atoms(const std::string &input_string, bool quiet) {
    Grammar<std::string::iterator, qi::ascii::space_type> grammar;
    return input_atom_selection(grammar, "Enter mask for read trajectory > ", input_string, false, quiet);
}

void Atom::select2group(Atom::AmberMask &ids1, Atom::AmberMask &ids2, const std::string &prompt1,
                        const std::string &prompt2) {
    Grammar<std::string::iterator, qi::ascii::space_type> grammar;

    ids1 = input_atom_selection(grammar, prompt1);
    ids2 = input_atom_selection(grammar, prompt2);
}

void Atom::select1group(AmberMask &ids, const std::string &prompt, bool allow_empty) {
    namespace qi = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;
    using ascii::char_;

    Grammar<std::string::iterator, qi::ascii::space_type> grammar;

    ids = input_atom_selection(grammar, prompt, "", allow_empty);
}

bool is_match_impl(const std::shared_ptr<Atom> &atom, const Atom::Node &ast) {
    return boost::apply_visitor(AtomEqual(atom), ast);
}

bool Atom::is_match(const std::shared_ptr<Atom> &atom, const Atom::AmberMask &id) { return is_match_impl(atom, id); }

bool AtomEqual::operator()(const std::shared_ptr<Atom::residue_name_nums> &residues) const {
    if (residues) {
        if (!atom->residue_name or !atom->residue_num) {
            throw std::runtime_error("residue selection syntax is invaild in current context");
        }
        struct Equal_residue : boost::static_visitor<bool> {
            explicit Equal_residue(const std::shared_ptr<Atom> &atom) : atom(atom) {}

            bool operator()(const fusion::vector<uint, boost::optional<std::pair<uint, int>>> &i) {
                auto op = fusion::at_c<1>(i);
                auto res_num = atom->residue_num.get();
                if (op) {
                    if (op.get().first >= fusion::at_c<0>(i)) {
                        return res_num >= fusion::at_c<0>(i) and res_num <= op.get().first and
                               ((res_num - fusion::at_c<0>(i)) % op.get().second == 0);
                    } else {
                        return res_num >= op.get().first and res_num <= fusion::at_c<0>(i) and
                               ((res_num - fusion::at_c<0>(i)) % op.get().second == 0);
                    }
                } else {
                    return res_num == fusion::at_c<0>(i);
                }
                return false;
            }

            bool operator()(const Atom::Name &pattern) {
                if (pattern.has_GLOB) {
                    if (fnmatch(pattern.name.c_str(), atom->residue_name.get().c_str(), 0) == 0)
                        return true;
                    if (!pattern.has_alpha) {
                        std::string num_str = std::to_string(atom->residue_num.get());
                        return fnmatch(pattern.name.c_str(), num_str.c_str(), 0) == 0;
                    }
                } else if (pattern.has_alpha) {
                    return pattern.name == atom->residue_name.get();
                }
                return false;
            }

        private:
            const std::shared_ptr<Atom> &atom;
        } equal(atom);

        for (auto &i : residues->val) {
            if (boost::apply_visitor(equal, i))
                return true;
        }
    }
    return false;
}

bool AtomEqual::operator()(const std::shared_ptr<Atom::molecule_nums> &molecules) const {
    if (molecules) {
        auto Equal_molecule = [this](const auto &i) -> bool {
            auto seq = atom->molecule.lock()->sequence;
            auto op = fusion::at_c<1>(i);
            if (op) {
                auto step = fusion::at_c<1>(op.get()) ? fusion::at_c<1>(op.get()).get() : 1;
                if (fusion::at_c<0>(op.get()) >= fusion::at_c<0>(i)) {
                    return seq >= fusion::at_c<0>(i) and seq <= fusion::at_c<0>(op.get()) and
                           ((seq - fusion::at_c<0>(i)) % step == 0);
                } else {
                    return seq >= fusion::at_c<0>(op.get()) and seq <= fusion::at_c<0>(i) and
                           ((seq - fusion::at_c<0>(i)) % step == 0);
                }
            } else {
                return seq == fusion::at_c<0>(i);
            }
            return false;
        };

        for (auto &i : molecules->val) {
            if (Equal_molecule(i))
                return true;
        }
    }
    return false;
}

bool AtomEqual::operator()(const std::shared_ptr<Atom::atom_name_nums> &names) const {
    if (names) {
        struct Equal_atom : boost::static_visitor<bool> {
            explicit Equal_atom(const std::shared_ptr<Atom> &atom) : atom(atom) {}

            bool operator()(const fusion::vector<uint, boost::optional<std::pair<uint, int>>> &i) {
                auto op = fusion::at_c<1>(i);
                if (op) {
                    if (op.get().first >= fusion::at_c<0>(i))
                        return atom->seq >= fusion::at_c<0>(i) and atom->seq <= op.get().first and
                               ((atom->seq - fusion::at_c<0>(i)) % op.get().second) == 0;
                    else
                        return atom->seq >= op.get().first and atom->seq <= fusion::at_c<0>(i) and
                               ((atom->seq - fusion::at_c<0>(i)) % op.get().second) == 0;
                } else {
                    return atom->seq == fusion::at_c<0>(i);
                }
            }

            bool operator()(const Atom::Name &pattern) {
                if (pattern.has_GLOB) {
                    if (fnmatch(pattern.name.c_str(), atom->atom_name.c_str(), 0) == 0)
                        return true;
                    if (!pattern.has_alpha) {
                        std::string num_str = std::to_string(atom->seq);
                        return fnmatch(pattern.name.c_str(), num_str.c_str(), 0) == 0;
                    }
                } else if (pattern.has_alpha) {
                    return pattern.name == atom->atom_name;
                }
                return false;
            }

        private:
            const std::shared_ptr<Atom> &atom;
        } equal(atom);

        for (auto &i : names->val) {
            if (boost::apply_visitor(equal, i))
                return true;
        }
    }
    return false;
}

bool AtomEqual::operator()(const std::shared_ptr<Atom::atom_types> &types) const {
    if (types) {
        struct Equal_types : boost::static_visitor<bool> {
            explicit Equal_types(const std::shared_ptr<Atom> &atom) : atom(atom) {}

            bool operator()(const fusion::vector<uint, boost::optional<std::pair<uint, int>>> &i) {
                auto op = fusion::at_c<1>(i);
                if (op) {
                    if (op.get().first >= fusion::at_c<0>(i))
                        return atom->typ >= static_cast<int>(fusion::at_c<0>(i)) and
                               atom->typ <= static_cast<int>(op.get().first) and
                               ((atom->typ - fusion::at_c<0>(i)) % op.get().second) == 0;
                    else
                        return atom->typ >= static_cast<int>(op.get().first) and
                               atom->typ <= static_cast<int>(fusion::at_c<0>(i)) and
                               ((atom->typ - fusion::at_c<0>(i)) % op.get().second) == 0;
                } else {
                    return atom->typ == static_cast<int>(fusion::at_c<0>(i));
                }
            }

            bool operator()(const Atom::Name &pattern) {
                if (pattern.has_GLOB) {
                    if (fnmatch(pattern.name.c_str(), atom->type_name.c_str(), 0) == 0)
                        return true;
                    if (!pattern.has_alpha) {
                        std::string num_str = std::to_string(atom->typ);
                        return fnmatch(pattern.name.c_str(), num_str.c_str(), 0) == 0;
                    }
                } else if (pattern.has_alpha) {
                    return pattern.name == atom->type_name;
                }
                return false;
            }

        private:
            const std::shared_ptr<Atom> &atom;
        } equal(atom);

        for (auto &i : types->val) {
            if (boost::apply_visitor(equal, i))
                return true;
        }
    }
    return false;
}

bool AtomEqual::operator()(const std::shared_ptr<Atom::atom_element_names> &ele) const {
    if (ele) {
        if (!atom->atom_symbol) {
            throw std::runtime_error("atom element symbol selection syntax is invaild in current context");
        }
        for (auto &pattern : ele->val) {
            if (pattern.has_GLOB) {
                if (fnmatch(pattern.name.c_str(), atom->atom_symbol.get().c_str(), 0) == 0)
                    return true;
            } else {
                if (pattern.name == atom->atom_symbol.get())
                    return true;
            }
        }
    }
    return false;
}

bool AtomEqual::operator()(const std::shared_ptr<Atom::Operator> &op) const {
    if (op) {
        switch (op->op) {
        case Atom::Op::NOT:
            return not boost::apply_visitor(AtomEqual(atom), op->node1);
            break;
        case Atom::Op::AND: {
            AtomEqual equal(atom);
            return boost::apply_visitor(equal, op->node1) and boost::apply_visitor(equal, op->node2);
        } break;
        case Atom::Op::OR: {
            AtomEqual equal(atom);
            return boost::apply_visitor(equal, op->node1) or boost::apply_visitor(equal, op->node2);
        } break;
        default:
            throw std::runtime_error("invalid Operator");
        }
    }
    return false;
}

std::string to_string(const AmberMask &mask) {
    std::string generated;
    throw_assert(format_node(mask, generated), "amberMask format error");
    return generated;
}
