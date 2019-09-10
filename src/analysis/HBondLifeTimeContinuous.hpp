//
// Created by xiamr on 9/10/19.
//

#ifndef TINKER_HBONDLIFETIMECONTINUOUS_HPP
#define TINKER_HBONDLIFETIMECONTINUOUS_HPP

#include "HBondLifeTime.hpp"

class HBondLifeTimeContinuous : public HBondLifeTime {
public:
    void print(std::ostream &os) override;

    static std::string title() { return "Hydrogen Bond LifeTime for Water(Si, history dependent)"; }

protected:
    [[nodiscard]] std::vector<double> calculateAcf() const override;
};


#endif //TINKER_HBONDLIFETIMECONTINUOUS_HPP
