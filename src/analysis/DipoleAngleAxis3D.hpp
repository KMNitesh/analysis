//
// Created by xiamr on 7/8/19.
//

#ifndef TINKER_DIPOLEANGLEAXIS3D_HPP
#define TINKER_DIPOLEANGLEAXIS3D_HPP

#include <unordered_set>
#include <memory>
#include <list>
#include "AbstractAnalysis.hpp"
#include "common.hpp"
#include "atom.hpp"

class Frame;

class DipoleAngleAxis3D : public AbstractAnalysis {
public:
    DipoleAngleAxis3D() {
        enable_forcefield = true;
        enable_outfile = true;
    }

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    static const std::string title() { return "Dipole Angle 3 Dimension Distribution with Axis"; }

protected:

    Atom::AmberMask ids;

    std::unordered_set<std::shared_ptr<Molecule>> group;

    void printData(std::ostream &os) const;

    std::list<std::tuple<double, double, double>> distributions;
};


#endif //TINKER_DIPOLEANGLEAXIS3D_HPP
