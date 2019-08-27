//
// Created by xiamr on 8/26/19.
//

#ifndef TINKER_CLUSTERVOLUME_HPP
#define TINKER_CLUSTERVOLUME_HPP

#include "std.hpp"
#include "BasicAnalysis.hpp"
#include "atom.hpp"
#include "HBond.hpp"

class Frame;

class ClusterVolume : public BasicAnalysis {
public:
    ClusterVolume();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    static std::string title() { return "vdW Volume"; }

protected:

    static double getVdwRadii(const std::shared_ptr<Atom> &atom);

    AmberMask atom_mask;

    std::unordered_set<std::shared_ptr<Atom>> atom_group;

    int grid_x;
    int grid_y;
    int grid_z;

    std::deque<std::pair<double, double>> volumes;

    inline static std::unordered_map<Symbol, double> vdWRadiis{
            {Symbol::Hydrogen,   1.10},
            {Symbol::Carbon,     1.70},
            {Symbol::Nitrogen,   1.55},
            {Symbol::Oxygen,     1.52},
            {Symbol::Phosphorus, 1.80},
            {Symbol::Sulfur,     1.80}
    };
};


#endif //TINKER_CLUSTERVOLUME_HPP
