//
// Created by xiamr on 7/3/19.
//

#ifndef TINKER_ROTACF_HPP
#define TINKER_ROTACF_HPP

#include "AbstractAnalysis.hpp"
#include "utils/VectorSelector.hpp"

class Frame;

class Molecule;

class RotAcf : public AbstractAnalysis {
public:
    explicit RotAcf();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    [[nodiscard]] std::string description();

    void readInfo() override;

    void setParameters(const std::shared_ptr<VectorSelector> &vector, int LegendrePolynomial, double time_increment_ps,
                       double max_time_gap_ps, const std::string &outfilename);

    [[nodiscard]] static std::string_view title() { return "Rotational Autocorrelation Function"; }

protected:
    template <typename Function>
    [[nodiscard]] std::vector<double> calculate(Function f) const;

    [[nodiscard]] std::vector<double> integrate(const std::vector<double> &acf) const;

    std::vector<std::vector<std::tuple<double, double, double>>> rots;

    std::shared_ptr<VectorSelector> vectorSelector;

    int LegendrePolynomial;

    double time_increment_ps = 0.1;

    double max_time_gap;
};

#endif  // TINKER_ROTACF_HPP
