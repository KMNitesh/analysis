
#ifndef TINKER_BONDEDENERGYCALC_HPP
#define TINKER_BONDEDENERGYCALC_HPP

#include "ana_module/AbstractAnalysis.hpp"
#include "data_structure/atom.hpp"
#include "dsl/AmberMask.hpp"
#include "utils/BondEnergyCalculator.hpp"

class Frame;

class BondedEnergyCalc : public AbstractAnalysis {
public:
    BondedEnergyCalc();
    [[nodiscard]] static std::string_view title() { return "Bonded Energy Calculator"; }

    virtual void processFirstFrame([[maybe_unused]] std::shared_ptr<Frame> &frame);

    virtual void process(std::shared_ptr<Frame> &frame);

    virtual void print(std::ostream &os);

    virtual void readInfo();

private:
    AmberMask mask;
    std::unique_ptr<BondEnergyCalculator> calculator;
    std::deque<std::map<int, BondEnergyCalculator::Term>> energy_terms;
};

#endif  // TINKER_BONDEDENERGYCALC_HPP