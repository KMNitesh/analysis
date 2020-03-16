//
// Created by xiamr on 3/17/19.
//

#ifndef TINKER_ATOM_HPP
#define TINKER_ATOM_HPP

#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/sequence/intrinsic/at_c.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/optional.hpp>
#include <boost/optional/optional_io.hpp>
#include <boost/phoenix/function/adapt_function.hpp>
#include <boost/variant.hpp>

#include "utils/common.hpp"
#include "utils/std.hpp"

class Molecule;

struct lj_t {
    double c6,c12;
};

class Atom {
    boost::optional<int> at_no;

public:
    boost::graph_traits<graph_t>::vertex_descriptor vertex_descriptor;
    std::size_t seq;
    std::string atom_name;

    const boost::optional<int> &getAtNo() const { return at_no; }

    void setAtNo(const boost::optional<int> &atNo) { at_no = atNo; }

    double x, y, z;  // position

    std::tuple<double, double, double> getCoordinate() const { return {x, y, z}; }

    std::tuple<double, double, double> getVelocities() const { return {vx, vy, vz}; }

    double vx = 0.0, vy = 0.0, vz = 0.0;  // velocity

    int typ;  // atom type
    std::string type_name;

    boost::optional<double> charge;

    boost::optional<double> mass;

    boost::optional<lj_t>  lj_param;

    std::list<std::size_t> con_list;  // atom num that connect to

    std::weak_ptr<Molecule> molecule;

    boost::optional<std::string> residue_name;
    boost::optional<uint> residue_num;

    boost::optional<std::string> atom_symbol;

    bool mark = false;  // used in NMR analysis

    bool adj(const std::shared_ptr<Atom> &atom) {
        for (auto i : con_list) {
            if (atom->seq == i) return true;
        }
        return false;
    }

    using numItemType = boost::fusion::vector<uint, boost::optional<std::pair<uint, int>>>;

    using select_ranges = std::vector<boost::variant<numItemType, std::string>>;

    struct residue_name_nums {
        select_ranges val;

        explicit residue_name_nums(const select_ranges &val) : val(val) {}
    };

    struct molecule_nums {
        using Attr = boost::fusion::vector<uint, boost::optional<boost::fusion::vector<uint, boost::optional<int>>>>;
        std::vector<Attr> val;

        explicit molecule_nums(const std::vector<Attr> &val) : val(val) {}
    };

    struct atom_name_nums {
        select_ranges val;

        explicit atom_name_nums(const select_ranges &val) : val(val) {}

        explicit atom_name_nums(const std::string &name) { val.emplace_back(name); }
    };

    struct atom_types {
        select_ranges val;

        explicit atom_types(const select_ranges &val) : val(val) {}

        explicit atom_types(const std::string &type) { val.push_back(type); }

        explicit atom_types(int typenum) {
            boost::fusion::vector<uint, boost::optional<std::pair<uint, int>>> t;
            boost::fusion::at_c<0>(t) = typenum;
            val.emplace_back(t);
        }
    };

    struct atom_element_names {
        std::vector<std::string> val;

        explicit atom_element_names(const std::vector<std::string> &val) : val(val) {}
    };

    struct Operator;

    using Node = boost::variant<std::shared_ptr<Operator>, std::shared_ptr<residue_name_nums>,
                                std::shared_ptr<molecule_nums>, std::shared_ptr<atom_name_nums>,
                                std::shared_ptr<atom_types>, std::shared_ptr<atom_element_names>>;

    enum class Op { NOT, AND, OR };

    struct Operator {
        explicit Operator(Op op, Node node1 = Node(), Node node2 = Node()) : op(op), node1(node1), node2(node2) {}

        Op op;
        Node node1;
        Node node2;
    };

    using AmberMask = Atom::Node;

    static bool is_match(const std::shared_ptr<Atom> &atom, const AmberMask &id);

    static void select2group(Atom::AmberMask &ids1, Atom::AmberMask &ids2,
                             const std::string &prompt1 = "Enter mask for atom1 : ",
                             const std::string &prompt2 = "Enter mask for atom2 : ");

    static void select1group(AmberMask &ids, const std::string &prompt = "Enter mask for atom : ");
};

using AmberMask = Atom::AmberMask;

struct print : boost::static_visitor<> {
    int space_num;

    explicit print(int space_num = 0) : space_num(space_num) {}

    void indent(int space_num) const;

    void operator()(const std::shared_ptr<Atom::residue_name_nums> &residues) const;

    void operator()(const std::shared_ptr<Atom::molecule_nums> &molecules) const;

    void operator()(const std::shared_ptr<Atom::atom_name_nums> &names) const;

    void operator()(const std::shared_ptr<Atom::atom_types> &types) const;

    void operator()(const std::shared_ptr<Atom::atom_element_names> &ele) const;

    void operator()(const std::shared_ptr<Atom::Operator> &op) const;
};

struct AtomEqual : boost::static_visitor<bool> {
    explicit AtomEqual(const std::shared_ptr<Atom> &atom) : atom(atom) {}

    bool operator()(const std::shared_ptr<Atom::residue_name_nums> &residues) const;

    bool operator()(const std::shared_ptr<Atom::molecule_nums> &molecules) const;

    bool operator()(const std::shared_ptr<Atom::atom_name_nums> &names) const;

    bool operator()(const std::shared_ptr<Atom::atom_types> &types) const;

    bool operator()(const std::shared_ptr<Atom::atom_element_names> &ele) const;

    bool operator()(const std::shared_ptr<Atom::Operator> &op) const;

private:
    const std::shared_ptr<Atom> &atom;
};

std::ostream &operator<<(std::ostream &os, const Atom::AmberMask &mask);

inline bool operator==(const std::shared_ptr<Atom::residue_name_nums> &residues1,
                       const std::shared_ptr<Atom::residue_name_nums> &residues2) {
    if (residues1 && residues2) {
        return residues1->val == residues2->val;
    }
    return false;
}

inline bool operator==(const std::shared_ptr<Atom::molecule_nums> &molecules1,
                       const std::shared_ptr<Atom::molecule_nums> &molecules2) {
    if (molecules1 && molecules2) {
        return molecules1->val == molecules2->val;
    }
    return false;
}

inline bool operator==(const std::shared_ptr<Atom::atom_name_nums> &names1,
                       const std::shared_ptr<Atom::atom_name_nums> &names2) {
    if (names1 && names2) {
        return names1->val == names2->val;
    }
    return false;
}

inline bool operator==(const std::shared_ptr<Atom::atom_types> &types1,
                       const std::shared_ptr<Atom::atom_types> &types2) {
    if (types1 && types2) {
        return types1->val == types2->val;
    }
    return false;
}

inline bool operator==(const std::shared_ptr<Atom::atom_element_names> &ele1,
                       const std::shared_ptr<Atom::atom_element_names> &ele2) {
    if (ele1 && ele2) {
        return ele1->val == ele2->val;
    }
    return false;
}

inline bool operator==(const std::shared_ptr<Atom::Operator> &op1, const std::shared_ptr<Atom::Operator> &op2) {
    if (op1 && op2) {
        if (op1->op == op2->op) {
            if (op1->op == Atom::Op::NOT) {
                return op1->node1 == op2->node1;
            } else {
                return op1->node1 == op2->node1 and op1->node2 == op2->node2;
            }
        }
    }
    return false;
}

std::string to_string(const AmberMask &mask);

#endif  // TINKER_ATOM_HPP
