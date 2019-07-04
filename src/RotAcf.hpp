//
// Created by xiamr on 7/3/19.
//

#ifndef TINKER_ROTACF_HPP
#define TINKER_ROTACF_HPP

#include <memory>
#include <unordered_set>
#include <string>
#include <map>
#include <list>
#include <tuple>
#include <vector>

#include "common.hpp"
#include "BasicAnalysis.hpp"
#include "atom.hpp"

class Frame;

class Molecule;

class RotAcf : public BasicAnalysis {
public:

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    explicit RotAcf() {
        enable_outfile = true;
        enable_tbb = true;
    }

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    static const std::string title() { return "Rotational autocorrelation function"; }

protected:

    double time_increment_ps = 0.1;

    Atom::AtomIndenter ids1;
    Atom::AtomIndenter ids2;
    Atom::AtomIndenter ids3;

    std::list<std::tuple<std::shared_ptr<Atom>, std::shared_ptr<Atom>, std::shared_ptr<Atom>>> pairs;

    std::vector<std::vector<std::tuple<double, double, double>>> rots;

    std::tuple<double, double, double> calVector(
            std::tuple<std::shared_ptr<Atom>, std::shared_ptr<Atom>, std::shared_ptr<Atom>> &atoms,
            std::shared_ptr<Frame> &frame);

    std::vector<double> calculate() const;

    std::vector<double> integrate(const std::vector<double> &acf) const;
};


#endif //TINKER_ROTACF_HPP
