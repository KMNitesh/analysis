//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_DIPOLEANGLESINGLEDISTANCENORMAL_HPP
#define TINKER_DIPOLEANGLESINGLEDISTANCENORMAL_HPP

#include "DipoleAngle.hpp"

class DipoleAngleSingleDistanceNormal : public DipoleAngle {
public:
    void print(std::ostream &os) override;

    [[nodiscard]] static std::string_view title() { return "Dipole Angle of single distance normal"; };

};

#endif //TINKER_DIPOLEANGLESINGLEDISTANCENORMAL_HPP
