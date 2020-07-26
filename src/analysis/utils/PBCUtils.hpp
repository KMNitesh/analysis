

#ifndef TINKER_PBCUTILS_HPP
#define TINKER_PBCUTILS_HPP

#include <boost/graph/breadth_first_search.hpp>
#include <boost/range/numeric.hpp>

#include "ana_module/Trajconv.hpp"
#include "data_structure/atom.hpp"
#include "data_structure/frame.hpp"
#include "data_structure/molecule.hpp"

class Frame;

class Visitor : public boost::default_bfs_visitor {
public:
    explicit Visitor(std::vector<std::shared_ptr<Molecule>> &mols,
                     std::vector<std::pair<std::shared_ptr<Molecule>, std::shared_ptr<Molecule>>> &mole_seq)
        : mols(mols), mole_seq(mole_seq){};

    template <typename Edge, typename Graph>
    void tree_edge(Edge e, const Graph &g) {
        auto source = boost::source(e, g);
        auto target = boost::target(e, g);
        mole_seq.emplace_back(mols[source], mols[target]);
    }

private:
    std::vector<std::shared_ptr<Molecule>> &mols;
    std::vector<std::pair<std::shared_ptr<Molecule>, std::shared_ptr<Molecule>>> &mole_seq;
};

class AggregateVisitor : public boost::default_bfs_visitor {
public:
    explicit AggregateVisitor(const std::shared_ptr<Frame> &frame) : frame(frame){};

    template <typename Edge, typename Graph>
    void tree_edge(Edge e, const Graph &g) const {
        auto &source = g[boost::source(e, g)];
        auto &target = g[boost::target(e, g)];
        auto r = target->getCoordinate() - source->getCoordinate();
        frame->image(r);
        std::tie(target->x, target->y, target->z) = r + source->getCoordinate();
    }

private:
    const std::shared_ptr<Frame> &frame;
};

class PBCUtils {
public:
    void do_move_center_basedto_atom(const AmberMask &mask, const std::shared_ptr<Frame> &frame) const;

    void do_move_center_basedto_atom_group(const AmberMask &mask, const std::shared_ptr<Frame> &frame) const;

    void do_move_center_basedto_molecule(const AmberMask &mask, const std::shared_ptr<Frame> &frame) const;

    void do_move_all_atom_into_box(const std::shared_ptr<Frame> &frame) const;

    void do_molecule_aggregate(const std::shared_ptr<Frame> &frame) const;

    void doPBC(Trajconv::PBCType pbc_mode, const AmberMask &mask, const std::shared_ptr<Frame> &frame) const;

    static std::shared_ptr<Atom> find_atom(const AmberMask &mask, const std::shared_ptr<Frame> &frame);

    static std::vector<std::shared_ptr<Atom>> find_atoms(const AmberMask &mask, const std::shared_ptr<Frame> &frame);

    static std::shared_ptr<Molecule> find_molecule(const AmberMask &mask, const std::shared_ptr<Frame> &frame);

    using MolPair = std::pair<std::set<std::shared_ptr<Molecule>>,
                              std::vector<std::pair<std::shared_ptr<Molecule>, std::shared_ptr<Molecule>>>>;

    template <typename Iterator>
    static MolPair calculate_intermol(Iterator &&atoms, const std::shared_ptr<Frame> &frame) {
        std::set<std::shared_ptr<Molecule>> mols_set;
        for (auto &atom : atoms) {
            mols_set.insert(atom->molecule.lock());
        }
        return {mols_set, calculate_intermol_imp(mols_set, frame)};
    }

    static void move(std::set<std::shared_ptr<Molecule>> &mols_set,
                     std::vector<std::pair<std::shared_ptr<Molecule>, std::shared_ptr<Molecule>>> &mol_seq,
                     const std::shared_ptr<Frame> &frame);

    static void move(MolPair &mols, const std::shared_ptr<Frame> &frame) { move(mols.first, mols.second, frame); };

private:
    static std::vector<std::pair<std::shared_ptr<Molecule>, std::shared_ptr<Molecule>>> calculate_intermol_imp(
        const std::set<std::shared_ptr<Molecule>> &mols_set, const std::shared_ptr<Frame> &frame);

    static std::tuple<double, double, double> cal_box_center_shift(const std::shared_ptr<Frame> &frame);

    mutable std::vector<std::shared_ptr<Atom>> selected_atoms;
    mutable std::shared_ptr<Molecule> selected_molecule;
    mutable std::set<std::shared_ptr<Molecule>> mols_set;
};

#endif  // TINKER_PBCUTILS_HPP
