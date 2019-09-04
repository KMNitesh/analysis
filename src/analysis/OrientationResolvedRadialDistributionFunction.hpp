//
// Created by xiamr on 9/4/19.
//

#ifndef TINKER_ORIENTATIONRESOLVEDRADIALDISTRIBUTIONFUNCTION_HPP
#define TINKER_ORIENTATIONRESOLVEDRADIALDISTRIBUTIONFUNCTION_HPP

#include "std.hpp"
#include <boost/multi_array.hpp>
#include "BasicAnalysis.hpp"
#include "atom.hpp"
#include "VectorSelector.hpp"

class Frame;

class OrientationResolvedRadialDistributionFunction : public BasicAnalysis {
public:
    OrientationResolvedRadialDistributionFunction();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    static std::string title() {
        return "Orientation-Resolved Radial Distribution Functions (JCTC 2019, 15, 803−812)";
    }

protected:

    AmberMask reference_atom_mask;
    AmberMask water_Ow_atom_mask;

    std::shared_ptr<Atom> reference_atom;
    std::vector<std::shared_ptr<Atom>> water_OW_atoms;

    std::shared_ptr<VectorSelector> vectorSelector;


    double distance_width;
    double angle_width;

    double max_distance;

    boost::multi_array<double,2> hist;

    int distance_bins;
    int angle_bins;

    int nframes = 0;

    void normalize();

    void write(std::ostream &os) const;
};


#endif //TINKER_ORIENTATIONRESOLVEDRADIALDISTRIBUTIONFUNCTION_HPP
