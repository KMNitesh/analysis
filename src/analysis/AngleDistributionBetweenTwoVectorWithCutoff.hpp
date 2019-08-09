//
// Created by xiamr on 8/6/19.
//

#ifndef TINKER_ANGLEDISTRIBUTIONBETWEENTWOVECTORWITHCUTOFF_HPP
#define TINKER_ANGLEDISTRIBUTIONBETWEENTWOVECTORWITHCUTOFF_HPP

#include "std.hpp"
#include "BasicAnalysis.hpp"
#include "atom.hpp"
#include "Histrogram.hpp"
#include "VectorSelector.hpp"

class Frame;

class AngleDistributionBetweenTwoVectorWithCutoff : public BasicAnalysis {
public:

    AngleDistributionBetweenTwoVectorWithCutoff();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    std::string description() override;

    static const std::string title() { return "Angle Distribution Between Two Vector with Cutoff"; }

    void setParameters(const Atom::Node &M,
                       const Atom::Node &L,
                       std::shared_ptr<VectorSelector> vector1,
                       std::shared_ptr<VectorSelector> vector2,
                       double angle_max,
                       double angle_width,
                       double cutoff1,
                       double cutoff2,
                       const std::string &outfilename);

protected:

    Atom::AmberMask ids1;
    Atom::AmberMask ids2;

    std::unordered_set<std::shared_ptr<Atom>> group1;
    std::unordered_set<std::shared_ptr<Atom>> group2;

    Histrogram angle_hist;
    Histrogram cos_hist;

    double cutoff1, cutoff2;

    std::shared_ptr<VectorSelector> vector1, vector2;

    void init_cos_hist(double angle_max, double angle_width);
};


#endif //TINKER_ANGLEDISTRIBUTIONBETWEENTWOVECTORWITHCUTOFF_HPP
