//
// Created by xiamr on 9/24/19.
//

#ifndef TINKER_LOCALSTRUCTUREINDEXFORLIQUID_HPP
#define TINKER_LOCALSTRUCTUREINDEXFORLIQUID_HPP

#include "std.hpp"
#include "BasicAnalysis.hpp"
#include "atom.hpp"

class Frame;

class LocalStructureIndexForLiquid : public BasicAnalysis {
public:
    LocalStructureIndexForLiquid();

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    static std::string title() { return "Local Structure Index (LSI) for Liquid"; }

protected:

    double cutoff2;
    int r_index;

    std::deque<std::pair<double, double>> localStructureIndices; // (Ri, LSI)
};


#endif //TINKER_LOCALSTRUCTUREINDEXFORLIQUID_HPP
