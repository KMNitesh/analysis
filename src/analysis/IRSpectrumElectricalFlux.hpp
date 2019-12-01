//
// Created by xiamr on 8/13/19.
//

#ifndef TINKER_IRSPECTRUMELECTRICALFLUX_HPP
#define TINKER_IRSPECTRUMELECTRICALFLUX_HPP

#include "std.hpp"
#include "AbstractAnalysis.hpp"
#include "atom.hpp"

class Frame;

class IRSpectrumElectricalFlux : public AbstractAnalysis {
public:

    IRSpectrumElectricalFlux();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    [[nodiscard]] static std::string_view title() { return "Infrared radiation (IR) Spectrum from Electrical Flux"; }

protected:

    double time_increment_ps;

    std::deque<std::tuple<double, double, double>> electricalFlux;

    AmberMask selected_mols_mask;

    std::unordered_set<std::shared_ptr<Molecule>> selected_mols;
};


#endif //TINKER_IRSPECTRUMELECTRICALFLUX_HPP
