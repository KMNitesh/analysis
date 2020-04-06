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

#include "dsl/AmberMask.hpp"
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

    std::array<double, 2> get_lj_p() const;

    std::list<std::size_t> con_list; // atom num that connect to

    std::weak_ptr<Molecule> molecule;

    boost::optional<std::string> residue_name;
    boost::optional<uint> residue_num;

    boost::optional<std::string> atom_symbol;

    bool mark = false; // used in NMR analysis

    bool adj(const std::shared_ptr<Atom> &atom);
};

#endif // TINKER_ATOM_HPP
