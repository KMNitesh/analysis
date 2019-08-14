//
// Created by xiamr on 8/13/19.
//

#ifndef TINKER_IRSPECTRUMELECTRICALFLUX_HPP
#define TINKER_IRSPECTRUMELECTRICALFLUX_HPP

#include "std.hpp"
#include "BasicAnalysis.hpp"
#include "atom.hpp"

class Frame;

class IRSpectrumElectricalFlux : public BasicAnalysis {
public:
    IRSpectrumElectricalFlux();

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    static std::string title() { return "Infrared radiation (IR) Spectrum from Electrical Flux"; }

protected:

    double time_increment_ps;

    std::deque<std::tuple<double, double, double>> electricalFlux;


};


#endif //TINKER_IRSPECTRUMELECTRICALFLUX_HPP
