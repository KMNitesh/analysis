//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_RADICALDISTRIBTUIONFUNCTION_HPP
#define TINKER_RADICALDISTRIBTUIONFUNCTION_HPP

#include <unordered_set>
#include <map>
#include <memory>

#include "common.hpp"
#include "BasicAnalysis.hpp"
#include "atom.hpp"

class Frame;

class RadicalDistribtuionFunction : public BasicAnalysis {

    double rmax;
    double width;

    bool intramol = false;
    int nframe = 0;
    int numj = 0, numk = 0;

    int nbin;
    std::map<int, int> hist;
    std::map<int, double> gr, gs, integral;

    double xbox, ybox, zbox;

    Atom::AtomIndenter ids1;
    Atom::AtomIndenter ids2;

    std::unordered_set<std::shared_ptr<Atom>> group1;
    std::unordered_set<std::shared_ptr<Atom>> group2;

public:
    RadicalDistribtuionFunction() {
        enable_outfile = true;
    }

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    static const std::string title() {
        return "Radical Distribution Function";
    }
};


#endif //TINKER_RADICALDISTRIBTUIONFUNCTION_HPP
