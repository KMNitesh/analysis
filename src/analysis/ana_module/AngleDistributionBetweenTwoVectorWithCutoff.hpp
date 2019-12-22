//
// Created by xiamr on 8/6/19.
//

#ifndef TINKER_ANGLEDISTRIBUTIONBETWEENTWOVECTORWITHCUTOFF_HPP
#define TINKER_ANGLEDISTRIBUTIONBETWEENTWOVECTORWITHCUTOFF_HPP

#include "utils/std.hpp"
#include "AbstractAnalysis.hpp"
#include "data_structure/atom.hpp"
#include "utils/Histogram.hpp"
#include "utils/VectorSelector.hpp"

class Frame;

class AngleDistributionBetweenTwoVectorWithCutoff : public AbstractAnalysis {
public:

    AngleDistributionBetweenTwoVectorWithCutoff();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    [[nodiscard]] std::string description() override;

    [[nodiscard]] static std::string_view title() { return "Angle Distribution Between Two Vector with Cutoff"; }

    void setParameters(const AmberMask &M,
                       const AmberMask &L,
                       std::shared_ptr<VectorSelector> vector1,
                       std::shared_ptr<VectorSelector> vector2,
                       double angle_max,
                       double angle_width,
                       double cutoff1,
                       double cutoff2,
                       const std::string &outfilename);

protected:

    Atom::AmberMask metal_mask;
    Atom::AmberMask ligand_mask;

    std::unordered_set<std::shared_ptr<Atom>> group1;
    std::unordered_set<std::shared_ptr<Atom>> group2;

    Histogram angle_hist;
    Histogram cos_hist;

    double cutoff1, cutoff2;

    std::shared_ptr<VectorSelector> vector1, vector2;

    std::multimap<int, std::pair<int, double>> angle_evolution;

    int nframe = 0;

    void init_cos_hist(double angle_max, double angle_width);

    void saveJson(std::ostream &os) const;
};


#endif //TINKER_ANGLEDISTRIBUTIONBETWEENTWOVECTORWITHCUTOFF_HPP
