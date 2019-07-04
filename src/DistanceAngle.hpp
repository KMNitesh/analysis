//
// Created by xiamr on 7/4/19.
//

#ifndef TINKER_DISTANCEANGLE_HPP
#define TINKER_DISTANCEANGLE_HPP

#include <memory>
#include <unordered_set>
#include <string>
#include <map>
#include <unordered_map>
#include <list>
#include <tuple>
#include <vector>

#include "common.hpp"
#include "BasicAnalysis.hpp"
#include "atom.hpp"


class DistanceAngle : public BasicAnalysis {
public:

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    explicit DistanceAngle() {
        enable_outfile = true;
    }

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    static const std::string title() { return "Gibbs Free Energy of Distance Angle for AnO2"; }


protected:

    Atom::AtomIndenter id1, id2, id3;

    std::list<std::pair<std::shared_ptr<Atom>, std::shared_ptr<Atom>>> pairs;

    std::unordered_set<std::shared_ptr<Atom>> group3;


    double distance_width;
    double angle_width;

    int distance_bins;
    int angle_bins;

    double temperature;  // unit: K

    std::map<std::pair<int, int>, size_t> hist;

    void printData(std::ostream &os) const;

};


#endif //TINKER_DISTANCEANGLE_HPP