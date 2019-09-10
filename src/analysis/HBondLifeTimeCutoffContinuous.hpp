//
// Created by xiamr on 9/10/19.
//

#ifndef TINKER_HBONDLIFETIMECUTOFFCONTINUOUS_HPP
#define TINKER_HBONDLIFETIMECUTOFFCONTINUOUS_HPP

#include "HBondLifeTimeCutoff.hpp"

class HBondLifeTimeCutoffContinuous : public HBondLifeTimeCutoff {
public:
    void print(std::ostream &os) override;

    static std::string
    title() { return "Hydrogen Bond LifeTime for Water in Solvation Shell(Si, history dependent)"; }

protected:
    [[nodiscard]] std::vector<double> calculateAcf() const override;
};


#endif //TINKER_HBONDLIFETIMECUTOFFCONTINUOUS_HPP
