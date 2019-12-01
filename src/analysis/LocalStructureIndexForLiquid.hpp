//
// Created by xiamr on 9/24/19.
//

#ifndef TINKER_LOCALSTRUCTUREINDEXFORLIQUID_HPP
#define TINKER_LOCALSTRUCTUREINDEXFORLIQUID_HPP

#include "std.hpp"
#include "AbstractAnalysis.hpp"
#include "atom.hpp"

class Frame;

class LocalStructureIndexForLiquid : public AbstractAnalysis {
public:

    LocalStructureIndexForLiquid();

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    [[nodiscard]] static std::string_view title() { return "Local Structure Index (LSI) for Liquid"; }

protected:

    double cutoff2;
    int r_index;

    std::deque<std::pair<double, double>> localStructureIndices; // (Ri, LSI)
};


#endif //TINKER_LOCALSTRUCTUREINDEXFORLIQUID_HPP
