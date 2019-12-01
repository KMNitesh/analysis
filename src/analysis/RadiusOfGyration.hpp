//
// Created by xiamr on 9/9/19.
//

#ifndef TINKER_RADIUSOFGYRATION_HPP
#define TINKER_RADIUSOFGYRATION_HPP

#include "std.hpp"
#include "AbstractAnalysis.hpp"
#include "atom.hpp"

class Frame;

class RadiusOfGyration : public AbstractAnalysis {
public:

    RadiusOfGyration();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    [[nodiscard]] static std::string_view title() { return "Radius of gyration (mass-weighted)"; }

protected:

    AmberMask atomMask;
    std::deque<double> series;
    std::unordered_set<std::shared_ptr<Molecule>> moles;

    bool bIncludeHydrogen;
};


#endif //TINKER_RADIUSOFGYRATION_HPP
