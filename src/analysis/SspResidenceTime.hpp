//
// Created by xiamr on 9/16/19.
//

#ifndef TINKER_SSPRESIDENCETIME_HPP
#define TINKER_SSPRESIDENCETIME_HPP

#include "ResidenceTime.hpp"

class SspResidenceTime : public ResidenceTime {
public:
    void print(std::ostream &os) override;

    void readInfo() override;

    static std::string
    title() { return "Residence Time based on SSP approach (J. Phys. Chem. B, Vol. 112, No. 26, 2008)"; }

protected:

    void calculateSSP();
};


#endif //TINKER_SSPRESIDENCETIME_HPP
