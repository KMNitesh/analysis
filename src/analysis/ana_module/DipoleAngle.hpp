//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_DIPOLEANGLE_HPP
#define TINKER_DIPOLEANGLE_HPP


#include <memory>
#include <unordered_set>
#include <string>
#include <map>
#include <list>
#include <utility>
#include "AbstractAnalysis.hpp"
#include "data_structure/atom.hpp"

class Frame;

class DipoleAngle : public AbstractAnalysis {
public:

    DipoleAngle();

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void readInfo() override;

    [[nodiscard]] static std::string_view title() { return "Dipole Angle"; }

protected:

    Atom::AmberMask ids1;
    Atom::AmberMask ids2;

    std::unordered_set<std::shared_ptr<Atom>> group1;
    std::unordered_set<std::shared_ptr<Atom>> group2;


    double distance_width;
    double angle_width;

    int distance_bins;
    int angle_bins;

    std::map<std::pair<int, int>, size_t> hist;
};


#endif //TINKER_DIPOLEANGLE_HPP
