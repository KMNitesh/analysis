
#ifndef TINKER_AMBERMASK_HPP
#define TINKER_AMBERMASK_HPP

#include <boost/algorithm/string.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/sequence/intrinsic/at_c.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/optional.hpp>
#include <boost/optional/optional_io.hpp>
#include <boost/phoenix/function/adapt_function.hpp>
#include <boost/variant.hpp>

class Atom;

namespace AmberMaskAST {

using numItemType = boost::fusion::vector<uint, boost::optional<std::pair<uint, int>>>;

struct Name {
    std::string name;
    bool has_GLOB = false;
    bool has_alpha = false;

    Name() : Name("") {}
    Name(const char *str) : Name(std::string(str)) {}
    Name(std::string name) : name(std::move(name)) { check(); }

    operator std::string() { return name; }

private:
    void check();
};

using select_ranges = std::vector<boost::variant<numItemType, Name>>;

struct residue_name_nums {
    select_ranges val;

    explicit residue_name_nums(const select_ranges &val) : val(val) {}
};

struct real_residue_nums {
    std::vector<numItemType> val;

    explicit real_residue_nums(const std::vector<numItemType> &val) : val(val) {}
};

struct molecule_nums {
    using Attr = boost::fusion::vector<uint, boost::optional<boost::fusion::vector<uint, boost::optional<int>>>>;
    std::vector<Attr> val;

    explicit molecule_nums(const std::vector<Attr> &val) : val(val) {}
};

struct atom_name_nums {
    select_ranges val;

    explicit atom_name_nums(const select_ranges &val) : val(val) {}

    explicit atom_name_nums(const Name &name) { val.emplace_back(name); }
};

struct atom_types {
    select_ranges val;

    explicit atom_types(const select_ranges &val) : val(val) {}

    explicit atom_types(const Name &type) { val.push_back(type); }

    explicit atom_types(int typenum) {
        boost::fusion::vector<uint, boost::optional<std::pair<uint, int>>> t;
        boost::fusion::at_c<0>(t) = typenum;
        val.emplace_back(t);
    }
};

struct atom_element_names {
    std::vector<Name> val;

    explicit atom_element_names(const std::vector<Name> &val) : val(val) {}
};

struct Operator;

using AmberMask =
    boost::variant<boost::blank, std::shared_ptr<Operator>, std::shared_ptr<residue_name_nums>,
                   std::shared_ptr<real_residue_nums>, std::shared_ptr<molecule_nums>, std::shared_ptr<atom_name_nums>,
                   std::shared_ptr<atom_types>, std::shared_ptr<atom_element_names>>;

enum class Op { NOT, AND, OR };

struct Operator {
    explicit Operator(Op op, AmberMask node1 = AmberMask(), AmberMask node2 = AmberMask())
        : op(op), node1(node1), node2(node2) {}

    Op op;
    AmberMask node1;
    AmberMask node2;
};

bool is_match(const std::shared_ptr<Atom> &atom, const AmberMask &id);

void select2group(AmberMask &ids1, AmberMask &ids2, const std::string &prompt1 = "Enter mask for atom1 : ",
                  const std::string &prompt2 = "Enter mask for atom2 : ");

void select1group(AmberMask &ids, const std::string &prompt = "Enter mask for atom : ", bool allow_empty = false);

inline bool isBlank(const AmberMask &mask) { return mask.which() == 0; }

std::ostream &operator<<(std::ostream &os, const AmberMask &mask);

std::ostream &operator<<(std::ostream &os, const Name &name);

bool operator==(const Name &name1, const Name &name2);

bool operator==(const std::shared_ptr<residue_name_nums> &residues1,
                const std::shared_ptr<residue_name_nums> &residues2);

bool operator==(const std::shared_ptr<real_residue_nums> &residues1,
                const std::shared_ptr<real_residue_nums> &residues2);

bool operator==(const std::shared_ptr<molecule_nums> &molecules1, const std::shared_ptr<molecule_nums> &molecules2);

bool operator==(const std::shared_ptr<atom_name_nums> &names1, const std::shared_ptr<atom_name_nums> &names2);

bool operator==(const std::shared_ptr<atom_types> &types1, const std::shared_ptr<atom_types> &types2);

bool operator==(const std::shared_ptr<atom_element_names> &ele1, const std::shared_ptr<atom_element_names> &ele2);

bool operator==(const std::shared_ptr<Operator> &op1, const std::shared_ptr<Operator> &op2);

struct print : boost::static_visitor<> {
    int space_num;

    explicit print(int space_num = 0) : space_num(space_num) {}

    void indent(int space_num) const;

    void operator()(const boost::blank &) const {};

    void operator()(const std::shared_ptr<residue_name_nums> &residues) const;

    void operator()(const std::shared_ptr<real_residue_nums> &residues) const;

    void operator()(const std::shared_ptr<molecule_nums> &molecules) const;

    void operator()(const std::shared_ptr<atom_name_nums> &names) const;

    void operator()(const std::shared_ptr<atom_types> &types) const;

    void operator()(const std::shared_ptr<atom_element_names> &ele) const;

    void operator()(const std::shared_ptr<Operator> &op) const;
};

std::string to_string(const AmberMask &mask);

AmberMask parse_atoms(const std::string &input_string, bool quiet);

} // namespace AmberMaskAST

using AmberMask = AmberMaskAST::AmberMask;

#endif // TINKER_AMBERMASK_HPP