//
// Created by xiamr on 8/11/19.
//

#ifndef TINKER_IRSPECTRUM_HPP
#define TINKER_IRSPECTRUM_HPP

#include "std.hpp"
#include "AbstractAnalysis.hpp"
#include "atom.hpp"

class Frame;

class IRSpectrum : public AbstractAnalysis {
public:

    IRSpectrum();

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    [[nodiscard]] static std::string_view title() { return "Infrared radiation (IR) Spectrum"; }

    [[nodiscard]] static std::vector<double> calculateIntense(const std::vector<double> &acf, double time_increment_ps);

    template<typename Container>
    [[nodiscard]] static std::vector<double> calculateAcf(const Container &evolution);

    static void calculateSpectrum(const std::string &out);

    static void printData(std::ostream &os,
                          const std::deque<std::tuple<double, double, double>> &dipole_evolution,
                          double time_increment_ps,
                          boost::optional<AmberMask> mask = boost::optional<AmberMask>{});

protected:

    [[nodiscard]] std::tuple<double, double, double> getDipole(std::shared_ptr<Frame> &frame);

    double time_increment_ps;

    std::deque<std::tuple<double, double, double>> dipole_evolution;

    AmberMask selected_mols_mask;

    std::unordered_set<std::shared_ptr<Molecule>> selected_mols;
};


#endif //TINKER_IRSPECTRUM_HPP
