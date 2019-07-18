//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_DIFFUSECUTOFF_HPP
#define TINKER_DIFFUSECUTOFF_HPP


#include <memory>
#include <unordered_set>
#include <string>
#include <map>
#include <list>
#include <utility>

#include "common.hpp"
#include "BasicAnalysis.hpp"
#include "atom.hpp"

class Frame;

class DiffuseCutoff : public BasicAnalysis {

public:

    DiffuseCutoff() {
        enable_outfile = true;
        enable_forcefield = true;
    }

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void readInfo() override;

    std::string getOutfileName() override { return outfilename; }

    void setParameters(const Atom::Node &M, const Atom::Node &L,
                       double cutoff, double time_increment_ps, const std::string &outfilename);

    static const std::string title() {
        return "Diffusion Cutoff Coefficient by Einstein equation";
    }

    struct InnerAtom {
        int index;
        std::list<std::tuple<double, double, double>> *list_ptr = nullptr;

        InnerAtom(int index, std::list<std::tuple<double, double, double>> *list_ptr)
                : index(index), list_ptr(list_ptr) {}
    };

    struct InnerAtomHasher {
        typedef InnerAtom argument_type;
        typedef std::size_t result_type;

        result_type operator()(argument_type const &s) const noexcept {
            result_type
            const h1(std::hash<int>()
            (s.index));
            result_type
            const h2(std::hash<std::list<std::tuple<double, double, double>> *>()
            (s.list_ptr));
            return h1 ^ (h2 << 1); // or use boost::hash_combine (see Discussion)
        }
    };

private:

    double time_increment_ps = 0.1;
    double cutoff2;

    Atom::AtomIndenter ids1;
    Atom::AtomIndenter ids2;

    std::unordered_set<std::shared_ptr<Atom>> group1;
    std::unordered_set<std::shared_ptr<Atom>> group2;

    std::unordered_set<InnerAtom, InnerAtomHasher> inner_atoms;

    std::list<std::list<std::tuple<double, double, double>> *> rcm;

    auto find_in(int seq);

    std::string outfilename;

};


#endif //TINKER_DIFFUSECUTOFF_HPP
