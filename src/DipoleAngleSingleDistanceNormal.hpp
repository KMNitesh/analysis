//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_DIPOLEANGLESINGLEDISTANCENORMAL_HPP
#define TINKER_DIPOLEANGLESINGLEDISTANCENORMAL_HPP

#include "DipoleAngle.hpp"

class DipoleAngleSingleDistanceNormal : public DipoleAngle {
public:
    void print() override;

    static const std::string title() {
        return "Dipole Angle of single distance normal";
    };

};

#endif //TINKER_DIPOLEANGLESINGLEDISTANCENORMAL_HPP
