//
// Created by xiamr on 7/3/19.
//

#ifndef TINKER_ROTACF_HPP
#define TINKER_ROTACF_HPP

#include <memory>
#include <unordered_set>
#include <string>
#include <map>
#include <list>
#include <tuple>
#include <vector>
#include <functional>

#include "common.hpp"
#include "BasicAnalysis.hpp"
#include "atom.hpp"
#include "VectorSelector.hpp"

class Frame;

class Molecule;

class RotAcf : public BasicAnalysis {
public:

    explicit RotAcf() {
        enable_outfile = true;
        enable_tbb = true;
    }

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    std::string description() override;

    void readInfo() override;

    void setParameters(const std::shared_ptr<VectorSelector> &vector, int LegendrePolynomial,
                       double time_increment_ps, double max_time_grap_ps, std::string outfilename);

    static const std::string title() { return "Rotational Autocorrelation Function"; }

protected:

    double time_increment_ps = 0.1;

    std::vector<std::vector<std::tuple<double, double, double>>> rots;

    std::shared_ptr<VectorSelector> vectorSelector;

    template<typename Function>
    std::vector<double> calculate(Function f) const;

    std::vector<double> integrate(const std::vector<double> &acf) const;

    int LegendrePolynomial;
    double max_time_grap;

    const std::unordered_map<int, std::string> LegendreStr{
            {1, "P1 = x"},
            {2, "P2 = (1/2)(3x^2 - 1)"},
            {3, "P3 = (1/2)(5x^3 - 3x)"},
            {4, "P4 = (1/8)(35x^4 - 30x^2 + 3)"},
    };
};


#endif //TINKER_ROTACF_HPP
