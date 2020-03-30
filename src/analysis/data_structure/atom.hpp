//
// Created by xiamr on 3/17/19.
//

#ifndef TINKER_ATOM_HPP
#define TINKER_ATOM_HPP

#include <boost/algorithm/string.hpp>
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
    double c6, c12;
};

class Atom {
    boost::optional<int> at_no;

public:
    boost::graph_traits<graph_t>::vertex_descriptor vertex_descriptor;
    std::size_t seq;
    std::string atom_name;

    const boost::optional<int> &getAtNo() const { return at_no; }

    void setAtNo(const boost::optional<int> &atNo) { at_no = atNo; }

    double x, y, z; // position

    std::tuple<double, double, double> getCoordinate() const { return {x, y, z}; }

    std::tuple<double, double, double> getVelocities() const { return {vx, vy, vz}; }

    double vx = 0.0, vy = 0.0, vz = 0.0; // velocity

    int typ; // atom type
    std::string type_name;

    boost::optional<double> charge;

    boost::optional<double> mass;

    boost::optional<lj_t> lj_param;

    std::array<double, 2> get_lj_p() const {
        const auto c6 = (*lj_param).c6;
        const auto c12 = (*lj_param).c12;
        const auto sigma = c6 != 0.0 ? std::pow(c12 / c6, 1.0 / 6) : 0.0;
        const auto epsilon = c12 != 0.0 ? c6 * c6 / (4 * c12) : 0.0;
        return {sigma, epsilon};
    }

    std::list<std::size_t> con_list; // atom num that connect to

    std::weak_ptr<Molecule> molecule;

    boost::optional<std::string> residue_name;
    boost::optional<uint> residue_num;

    boost::optional<std::string> atom_symbol;

    bool mark = false; // used in NMR analysis

    bool adj(const std::shared_ptr<Atom> &atom) {
        for (auto i : con_list) {
            if (atom->seq == i)
                return true;
        }
        return false;
    }

    using numItemType = boost::fusion::vector<uint, boost::optional<std::pair<uint, int>>>;

    struct Name {
        std::string name;
        bool has_GLOB;
        bool has_alpha;

        Name() : Name("") {}
        Name(const char *str) : Name(std::string(str)) {}
        Name(std::string name) : name(std::move(name)) {
            for (char c : name) {
                if (boost::is_any_of("*?")(c))
                    has_GLOB = true;
                else if (boost::is_alpha()(c))
                    has_alpha = true;
                if (has_GLOB and has_alpha)
                    return;
            }
        }

        operator std::string() { return name; }
    };

    using select_ranges = std::vector<boost::variant<numItemType, Name>>;

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

    using Node = boost::variant<boost::blank, std::shared_ptr<Operator>, std::shared_ptr<residue_name_nums>,
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

    static void select1group(AmberMask &ids,
                             const std::string &prompt = "Enter mask for atom : ", bool allow_empty = false);

    static bool isBlank(const AmberMask &mask) { return mask.which() == 0; }
};

using AmberMask = Atom::AmberMask;

struct print : boost::static_visitor<> {
    int space_num;

    explicit print(int space_num = 0) : space_num(space_num) {}

    void indent(int space_num) const;

    void operator()(const boost::blank &) const {};

    void operator()(const std::shared_ptr<Atom::residue_name_nums> &residues) const;

    void operator()(const std::shared_ptr<Atom::molecule_nums> &molecules) const;

    void operator()(const std::shared_ptr<Atom::atom_name_nums> &names) const;

    void operator()(const std::shared_ptr<Atom::atom_types> &types) const;

    void operator()(const std::shared_ptr<Atom::atom_element_names> &ele) const;

    void operator()(const std::shared_ptr<Atom::Operator> &op) const;
};

struct AtomEqual : boost::static_visitor<bool> {
    explicit AtomEqual(const std::shared_ptr<Atom> &atom) : atom(atom) {}

    bool operator()(const boost::blank &) const { return false; };

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

inline std::ostream &operator<<(std::ostream &os, const Atom::Name &name) { return os << name.name; }

inline bool operator==(const Atom::Name &name1, const Atom::Name &name2) { return name1.name == name2.name; }

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

AmberMask parse_atoms(const std::string &input_string, bool quiet = false);

#endif // TINKER_ATOM_HPP
