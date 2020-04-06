//
// Created by xiamr on 7/3/19.
//

#include "RotAcf.hpp"

#include <tbb/tbb.h>

#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm_ext.hpp>
#include <boost/range/irange.hpp>
#include <functional>

#include "data_structure/frame.hpp"
#include "data_structure/molecule.hpp"
#include "utils/LegendrePolynomial.hpp"
#include "utils/ThrowAssert.hpp"
#include "utils/VectorSelectorFactory.hpp"
#include "utils/common.hpp"
#include "data_structure/atom.hpp"

RotAcf::RotAcf() {
    enable_outfile = true;
    enable_tbb = true;
}

void RotAcf::processFirstFrame(std::shared_ptr<Frame> &frame) { rots.resize(vectorSelector->initialize(frame)); }

void RotAcf::process(std::shared_ptr<Frame> &frame) {
    auto vectors = vectorSelector->calculateVectors(frame);
    auto it1 = vectors.begin();
    auto it2 = rots.begin();
    assert(vectors.size() == rots.size());
    for (; it1 != vectors.end(); ++it1, ++it2) {
        it2->push_back(*it1);
    }
}

void RotAcf::print(std::ostream &os) {
    const std::unordered_map<int, std::function<std::vector<double>()>> func_mapping{
        {1, [this] { return calculate(LegendrePolynomialLevel1()); }},
        {2, [this] { return calculate(LegendrePolynomialLevel2()); }},
        {3, [this] { return calculate(LegendrePolynomialLevel3()); }},
        {4, [this] { return calculate(LegendrePolynomialLevel4()); }}};
    auto acf = func_mapping.at(LegendrePolynomial)();
    auto integration = integrate(acf);
    os << description() << '\n';
    os << "    Time Gap      ACF               integrate\n";
    os << "      (ps)                            (ps)\n";
    for (auto t : boost::irange(acf.size())) {
        os << boost::format("%12.2f%18.14f%15.5f\n") % (t * time_increment_ps) % acf[t] % integration[t];
    }
}

std::vector<double> RotAcf::integrate(const std::vector<double> &acf) const {
    std::vector<double> integrate(acf.size());
    integrate[0] = 0.0;
    for (auto i : boost::irange(std::size_t(1), integrate.size())) {
        integrate[i] = integrate[i - 1] + 0.5 * (acf[i - 1] + acf[i]) * time_increment_ps;
    }
    return integrate;
}

namespace {
template <typename Function>
class RotAcfParallelBody {
public:
    const std::vector<std::vector<std::tuple<double, double, double>>> &rots;
    std::vector<double, tbb::tbb_allocator<double>> acf;

    size_t array_length;
    size_t max_time_grap_step;
    Function f;

    explicit RotAcfParallelBody(const std::vector<std::vector<std::tuple<double, double, double>>> &rots,
                                size_t max_time_grap_step, size_t array_length, Function f)
        : rots(rots), acf(array_length), array_length(array_length), max_time_grap_step(max_time_grap_step), f(f) {}

    RotAcfParallelBody(const RotAcfParallelBody &rhs, tbb::split)
        : rots(rhs.rots),
          acf(rhs.array_length),
          array_length(rhs.array_length),
          max_time_grap_step(rhs.max_time_grap_step),
          f(rhs.f) {}

    void join(const RotAcfParallelBody &rhs) {
        for (size_t i = 1; i < array_length; i++) {
            acf[i] += rhs.acf[i];
        }
    }

    void operator()(const tbb::blocked_range<std::size_t> &range) {
        for (auto index = range.begin(); index != range.end(); index++) {
            auto &_vector = rots[index];
            auto total_size = _vector.size();

            for (size_t i = 0; i < total_size - 1; i++) {
                for (size_t j = i + 1; j < std::min(total_size, max_time_grap_step + i); j++) {
                    auto m = j - i;

                    assert(i < _vector.size());
                    assert(j < _vector.size());

                    double cos = dot_multiplication(_vector[i], _vector[j]);

                    acf[m] += f(cos);
                }
            }
        }
    }
};
}  // namespace

template <typename Function>
std::vector<double> RotAcf::calculate(Function f) const {
    size_t max_time_grap_step = std::ceil(max_time_grap / time_increment_ps) + 1;

    auto array_length = std::min(rots[0].size(), max_time_grap_step);
    RotAcfParallelBody parallelBody(rots, max_time_grap_step, array_length, f);

    tbb::parallel_reduce(tbb::blocked_range(std::size_t(0), rots.size()), parallelBody);

    std::vector<double> acf(array_length);
    auto total_mols = rots.size();
    auto total_size = rots[0].size();
    for (size_t i = 1; i < array_length; i++) {
        acf[i] = parallelBody.acf[i] / total_mols / (total_size - i);
    }

    acf[0] = 1.0;
    return acf;
}

void RotAcf::readInfo() {
    vectorSelector = VectorSelectorFactory::getVectorSelector();
    vectorSelector->readInfo();
    std::cout << "Legendre Polynomial\n";
    boost::range::for_each(boost::irange<int>(1, LegendreStr.size() + 1),
                           [](auto i) { std::cout << std::to_string(i) << ". " << LegendreStr.at(i) << '\n'; });

    LegendrePolynomial = choose<int>(1, LegendreStr.size(), "select > ");

    time_increment_ps =
        choose(0.0, std::numeric_limits<double>::max(), "Enter the Time Increment in Picoseconds [0.1]:", Default(0.1));

    max_time_grap = choose(0.0, std::numeric_limits<double>::max(), "Enter the Max Time Grap in Picoseconds :");
}

void RotAcf::setParameters(const std::shared_ptr<VectorSelector> &vector, int LegendrePolynomial,
                           double time_increment_ps, double max_time_grap_ps, const std::string &outfilename) {
    vectorSelector = vector;
    if (!(LegendrePolynomial >= 1 and LegendrePolynomial <= LegendreStr.size())) {
        throw std::runtime_error("legendre polynomial must be in the rang ef 1.." + std::to_string(LegendreStr.size()));
    }

    this->LegendrePolynomial = LegendrePolynomial;

    if (time_increment_ps <= 0) {
        throw std::runtime_error("`time_increment_ps` must be postive");
    } else {
        this->time_increment_ps = time_increment_ps;
    }

    if (max_time_grap_ps <= 0) {
        throw std::runtime_error("`max_time_grap_ps` must be postive");
    } else if (max_time_grap_ps <= this->time_increment_ps) {
        throw std::runtime_error("`max_time_grap_ps` must be larger than `time_increment_ps`");
    } else {
        this->max_time_grap = max_time_grap_ps;
    }
    this->outfilename = outfilename;
    boost::trim(this->outfilename);
    if (this->outfilename.empty()) {
        throw std::runtime_error("outfilename cannot empty");
    }
}

std::string RotAcf::description() {
    std::stringstream ss;
    std::string title_line = "------ " + std::string(title()) + " ------";
    ss << title_line << "\n";
    ss << " vector            = " << vectorSelector->description() << "\n";
    ss << " P                 = " << LegendrePolynomial << "  [ " << LegendreStr.at(LegendrePolynomial) << " ] \n";
    ss << " time_increment_ps = " << time_increment_ps << " (ps)\n";
    ss << " max_time_grap_ps  = " << max_time_grap << " (ps)\n";
    ss << " outfilename       = " << outfilename << "\n";
    ss << std::string(title_line.size(), '-') << '\n';
    return ss.str();
}
