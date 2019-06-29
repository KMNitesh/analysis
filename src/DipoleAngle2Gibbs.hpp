//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_DIPOLEANGLE2GIBBS_HPP
#define TINKER_DIPOLEANGLE2GIBBS_HPP

#include "DipoleAngle.hpp"

class DipoleAngle2Gibbs : public DipoleAngle {
public:
    void print(std::ostream &os) override;

    void readInfo() override;

    static const std::string title() {
        return "Dipole Angle to Gibbs Free Energy";
    };

protected:

    const double kb = 1.380649e-23; // unit: J/K
    double temperature;  // unit: K
    const double avogadro_constant = 6.022140857e23;

    void printData(std::ostream &os) const;
};

#endif //TINKER_DIPOLEANGLE2GIBBS_HPP
