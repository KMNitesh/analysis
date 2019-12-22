//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_DIFFUSECUTOFF_HPP
#define TINKER_DIFFUSECUTOFF_HPP

#include "utils/std.hpp"
#include <boost/container_hash/hash.hpp>
#include "AbstractAnalysis.hpp"
#include "data_structure/atom.hpp"


class Frame;

class DiffuseCutoff : public AbstractAnalysis {
public:

    DiffuseCutoff();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    [[nodiscard]] std::string description() override;

    void setParameters(const Atom::Node &M, const Atom::Node &L,
                       double cutoff, double time_increment_ps, const std::string &outfilename);

    [[nodiscard]] static std::string_view title() {
        return "Self-Diffusion Coefficient Calculation based on Einstein Equation within Solvation Shell";
    }

    struct InnerAtom {
        int index;
        std::deque<std::tuple<double, double, double>> *list_ptr = nullptr;

        InnerAtom(int index, std::deque<std::tuple<double, double, double>> *list_ptr)
                : index(index), list_ptr(list_ptr) {}
    };

    struct InnerAtomHasher {
        typedef InnerAtom argument_type;
        typedef std::size_t result_type;

        result_type operator()(argument_type const &s) const noexcept {
            size_t seed = 0;
            boost::hash_combine(seed, s.index);
            boost::hash_combine(seed, s.list_ptr);
            return seed;
        }
    };

    ~DiffuseCutoff() override;

private:

    double time_increment_ps = 0.1;
    double cutoff2;

    AmberMask mask1;
    AmberMask mask2;

    std::unordered_set<std::shared_ptr<Atom>> group1;
    std::unordered_set<std::shared_ptr<Atom>> group2;

    std::unordered_set<InnerAtom, InnerAtomHasher> inner_atoms;

    std::deque<std::deque<std::tuple<double, double, double>> *> rcm;

    auto find_in(int seq);

};


#endif //TINKER_DIFFUSECUTOFF_HPP
