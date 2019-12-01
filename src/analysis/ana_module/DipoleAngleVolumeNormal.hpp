//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_DIPOLEANGLEVOLUMENORMAL_HPP
#define TINKER_DIPOLEANGLEVOLUMENORMAL_HPP

#include "DipoleAngle.hpp"

class DipoleAngleVolumeNormal : public DipoleAngle {
public:

    void print(std::ostream &os) override;

    [[nodiscard]] static std::string_view title() { return "Dipole Angle of volume normal"; }
};


#endif //TINKER_DIPOLEANGLEVOLUMENORMAL_HPP
