//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_RADICALDISTRIBTUIONFUNCTION_HPP
#define TINKER_RADICALDISTRIBTUIONFUNCTION_HPP

#include <map>
#include <memory>
#include <unordered_set>

#include "AbstractAnalysis.hpp"
#include "data_structure/atom.hpp"
#include "dsl/AmberMask.hpp"

class Frame;

class RadicalDistribtuionFunction : public AbstractAnalysis {
public:
    RadicalDistribtuionFunction();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    void setParameters(const AmberMask &id1, const AmberMask &id2, double max_dist, double width, bool intramol,
                       const std::string outfilename);

    [[nodiscard]] static std::string_view title() { return "Radical Distribution Function"; }

private:
    double rmax;
    double width;

    bool intramol = false;
    int nframe = 0;
    int numj = 0, numk = 0;

    int nbin;
    std::map<int, int> hist;
    std::map<int, double> gr, gs, integral;

    double volume;

    AmberMask ids1;
    AmberMask ids2;

    std::unordered_set<std::shared_ptr<Atom>> group1;
    std::unordered_set<std::shared_ptr<Atom>> group2;
};

#endif // TINKER_RADICALDISTRIBTUIONFUNCTION_HPP
