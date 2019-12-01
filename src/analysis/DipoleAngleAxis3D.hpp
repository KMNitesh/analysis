//
// Created by xiamr on 7/8/19.
//

#ifndef TINKER_DIPOLEANGLEAXIS3D_HPP
#define TINKER_DIPOLEANGLEAXIS3D_HPP

#include <unordered_set>
#include <memory>
#include <list>
#include "AbstractAnalysis.hpp"
#include "atom.hpp"

class Frame;

class DipoleAngleAxis3D : public AbstractAnalysis {
public:

    DipoleAngleAxis3D();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    [[nodiscard]] static std::string_view title() { return "Dipole Angle 3 Dimension Distribution with Axis"; }

protected:

    Atom::AmberMask amberMask;

    std::unordered_set<std::shared_ptr<Molecule>> group;

    void printData(std::ostream &os) const;

    std::list<std::tuple<double, double, double>> distributions;
};


#endif //TINKER_DIPOLEANGLEAXIS3D_HPP
