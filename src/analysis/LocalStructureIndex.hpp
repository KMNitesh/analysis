//
// Created by xiamr on 9/22/19.
//

#ifndef TINKER_LOCALSTRUCTUREINDEX_HPP
#define TINKER_LOCALSTRUCTUREINDEX_HPP

#include "std.hpp"
#include "BasicAnalysis.hpp"
#include "atom.hpp"

class Frame;

class LocalStructureIndex : public BasicAnalysis {
public:
    LocalStructureIndex();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    static std::string title() { return "Local Structure Index (LSI) (PCCP 2011,13, 19918-19924)"; }

    template<typename RandomAccessRange>
    static double calculateLSI(RandomAccessRange &distance_within_cutoff_range);

protected:
    AmberMask metal_mask;
    AmberMask Ow_atom_mask;

    std::shared_ptr<Atom> metal;
    std::vector<std::shared_ptr<Atom>> Ow_atoms;
    double cutoff2;

    std::deque<std::pair<double, double>> localStructureIndices;
};


#endif //TINKER_LOCALSTRUCTUREINDEX_HPP
