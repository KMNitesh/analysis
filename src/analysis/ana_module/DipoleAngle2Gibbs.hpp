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

    [[nodiscard]] static std::string_view title() { return "Dipole Angle to Gibbs Free Energy"; };

protected:
    double temperature;  // unit: K

    void printData(std::ostream &os) const;
};

#endif  // TINKER_DIPOLEANGLE2GIBBS_HPP
