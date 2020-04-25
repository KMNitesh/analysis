//
// Created by xiamr on 8/9/19.
//

#ifndef TINKER_HBONDLIFETIMECUTOFF_HPP
#define TINKER_HBONDLIFETIMECUTOFF_HPP

#include <boost/container_hash/hash.hpp>

#include "AbstractAnalysis.hpp"
#include "data_structure/atom.hpp"
#include "dsl/AmberMask.hpp"
#include "utils/std.hpp"

class Frame;

class HBondLifeTimeCutoff : public AbstractAnalysis {
public:
    HBondLifeTimeCutoff();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    [[nodiscard]] static std::string_view title() {
        return "Hydrogen Bond LifeTime for Water in Solvation Shell(Ci, history independent)";
    }

    struct InnerAtom {
        int index;
        std::shared_ptr<std::deque<int>> list_ptr1, list_ptr2;

        InnerAtom(int index, std::shared_ptr<std::deque<int>> list_ptr1,  std::shared_ptr<std::deque<int>> list_ptr2)
            : index(index), list_ptr1(std::move(list_ptr1)), list_ptr2(std::move(list_ptr2)) {}
    };

    struct InnerAtomHasher {
        typedef InnerAtom argument_type;
        typedef std::size_t result_type;

        result_type operator()(argument_type const &s) const noexcept {
            size_t seed = 0;
            boost::hash_combine(seed, s.index);
            boost::hash_combine(seed, s.list_ptr1);
            boost::hash_combine(seed, s.list_ptr2);
            return seed;
        }
    };

protected:
    AmberMask center_Metal_atom_mask;
    AmberMask Ow_atom_mask;

    double dist_R_cutoff;
    double angle_HOO_cutoff;

    double time_increment_ps;
    double max_time_grap_ps;

    double cutoff2;

    std::deque<std::shared_ptr<std::deque<int>>> hb_histroy;

    std::unordered_set<std::shared_ptr<Atom>> metal;

    std::unordered_set<InnerAtom, InnerAtomHasher> inner_atoms;

    //  Ow, Hw
    std::vector<std::pair<std::shared_ptr<Atom>, std::deque<std::shared_ptr<Atom>>>> water_struct;

    [[nodiscard]] virtual std::vector<double> calculateAcf() const;

    void printData(std::ostream &os, const std::vector<double> &acf, std::string_view title) const;

    auto find_in(int seq);

    int find_hbond(const std::shared_ptr<Atom> &o1, const std::shared_ptr<Atom> &hydrogen,
                   std::shared_ptr<Frame> &frame) const;
};

#endif  // TINKER_HBONDLIFETIMECUTOFF_HPP
