//
// Created by xiamr on 8/11/19.
//

#ifndef TINKER_IRSPECTRUM_HPP
#define TINKER_IRSPECTRUM_HPP

#include "std.hpp"
#include "BasicAnalysis.hpp"
#include "atom.hpp"

class Frame;

class IRSpectrum : public BasicAnalysis {
public:
    IRSpectrum();

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    static std::string title() { return "Infrared radiation (IR) Spectrum"; }

    static std::vector<double> calculateIntense(const std::vector<long double> &acf, double time_increment_ps);

    template<typename Container>
    static std::vector<long double> calculateAcf(const Container &evolution);

    static void calculateSpectrum(const std::string &out);

    static void printData(std::ostream &os,
                          const std::deque<std::tuple<double, double, double>> &dipole_evolution,
                          double time_increment_ps);

protected:

    double time_increment_ps;

    std::deque<std::tuple<double, double, double>> dipole_evolution;
};


#endif //TINKER_IRSPECTRUM_HPP
