//
// Created by xiamr on 7/3/19.
//

#include <functional>
#include <tbb/tbb.h>
#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm_ext.hpp>
#include <boost/range/irange.hpp>

#include "RotAcf.hpp"
#include "frame.hpp"
#include "molecule.hpp"
#include "ThrowAssert.hpp"
#include "VectorSelectorFactory.hpp"

using namespace std;

void RotAcf::processFirstFrame(std::shared_ptr<Frame> &frame) {
    rots.resize(vectorSelector->initialize(frame));
}

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
            {1, [this] { return calculate([](auto x) { return x; }); }},
            {2, [this] { return calculate([](auto x) { return 0.5 * (3 * x * x - 1); }); }},
            {3, [this] { return calculate([](auto x) { return 0.5 * (5 * x * x * x - 3 * x); }); }},
            {4, [this] { return calculate([](auto x) { return 1.0 / 8.0 * (35 * x * x * x * x - 30 * x * x + 3); }); }}
    };
    auto acf = func_mapping.at(LegendrePolynomial)();
    vector<double> integration = integrate(acf);
    os << description() << '\n';
    os << "    Time Gap      ACF               integrate\n";
    os << "      (ps)                            (ps)\n";
    for (std::size_t t : boost::irange(acf.size() + 1)) {
        os << boost::format("%12.2f%18.14f%15.5f\n") % (t * time_increment_ps) % acf[t] % integration[t];
    }
}

vector<double> RotAcf::integrate(const vector<double> &acf) const {
    vector<double> integrate(acf.size());
    integrate[0] = 0.0;
    for (size_t i = 1; i < integrate.size(); i++) {
        integrate[i] = integrate[i - 1] + 0.5 * (acf[i - 1] + acf[i]) * time_increment_ps;
    }
    return integrate;
}

template<typename Function>
class RotAcfParallelBody {
public:
    const vector<vector<tuple<double, double, double>>> &rots;
    double *acf = nullptr;
    unsigned long long *ntime = nullptr;

    size_t array_length;
    size_t max_time_grap_step;
    Function f;

    explicit RotAcfParallelBody(const vector<vector<tuple<double, double, double>>> &rots,
                                size_t max_time_grap_step, size_t array_length, Function f)
            : rots(rots), acf(new double[array_length]{}),
              ntime(new unsigned long long[array_length]{}),
              array_length(array_length),
              max_time_grap_step(max_time_grap_step), f(f) {}

    RotAcfParallelBody(const RotAcfParallelBody &rhs, tbb::split)
            : rots(rhs.rots), acf(new double[rhs.array_length]{}),
              ntime(new unsigned long long[rhs.array_length]{}),
              array_length(rhs.array_length),
              max_time_grap_step(rhs.max_time_grap_step), f(rhs.f) {}

    void join(const RotAcfParallelBody &rhs) {
        for (size_t i = 1; i < array_length; i++) {
            acf[i] += rhs.acf[i];
            ntime[i] += rhs.ntime[i];
        }
    }

    virtual ~RotAcfParallelBody() {
        delete[] acf;
        delete[] ntime;
    }

    void operator()(const tbb::blocked_range<int> &range) {
        for (auto index = range.begin(); index != range.end(); index++) {
            auto &_vector = rots[index];
            auto total_size = _vector.size();

            for (size_t i = 0; i < total_size - 1; i++) {
                for (size_t j = i + 1; j < total_size; j++) {
                    auto m = j - i;

                    assert(i < _vector.size());
                    assert(j < _vector.size());

                    if (m > max_time_grap_step) break;

                    double cos = dot_multiplication(_vector[i], _vector[j]);

                    acf[m] += f(cos);
                    ntime[m]++;
                }
            }
        }
    }
};

template<typename Function>
vector<double> RotAcf::calculate(Function f) const {
    size_t max_time_grap_step = std::ceil(max_time_grap / time_increment_ps);

    RotAcfParallelBody parallelBody(rots, max_time_grap_step, min(rots[0].size(), max_time_grap_step + 1), f);

    tbb::parallel_reduce(tbb::blocked_range<int>(0, rots.size()), parallelBody, tbb::auto_partitioner());

    std::vector<double> acf(max_time_grap_step + 1, 0);
    for (size_t i = 1; i <= max_time_grap_step; i++) {
        assert(parallelBody.ntime[i] > 0);
        acf[i] = parallelBody.acf[i] / parallelBody.ntime[i];
    }

    acf[0] = 1.0;
    return acf;
}

void RotAcf::readInfo() {

    vectorSelector = VectorSelectorFactory::getVectorSelector();
    vectorSelector->readInfo();
    std::cout << "Legendre Polynomial\n";
    boost::range::for_each(boost::irange<int>(1, LegendreStr.size() + 1),
                           [this](auto i) { cout << to_string(i) << ". " << LegendreStr.at(i) << '\n'; });

    LegendrePolynomial = choose<int>(1, LegendreStr.size(), "select > ");

    time_increment_ps = choose(0.0, std::numeric_limits<double>::max(),
                               "Enter the Time Increment in Picoseconds [0.1]:", true, 0.1);

    max_time_grap = choose(0.0, std::numeric_limits<double>::max(), "Enter the Max Time Grap in Picoseconds :");
}

void RotAcf::setParameters(const shared_ptr<VectorSelector> &vector, int LegendrePolynomial,
                           double time_increment_ps, double max_time_grap_ps, const string &outfilename) {

    vectorSelector = vector;
    if (!(LegendrePolynomial >= 1 and LegendrePolynomial <= LegendreStr.size())) {
        throw runtime_error("legendre polynomial must be in the rang ef 1.." + to_string(LegendreStr.size()));
    }

    this->LegendrePolynomial = LegendrePolynomial;

    if (time_increment_ps <= 0) {
        throw runtime_error("`time_increment_ps` must be postive");
    } else {
        this->time_increment_ps = time_increment_ps;
    }

    if (max_time_grap_ps <= 0) {
        throw runtime_error("`max_time_grap_ps` must be postive");
    } else if (max_time_grap_ps <= this->time_increment_ps) {
        throw runtime_error("`max_time_grap_ps` must be larger than `time_increment_ps`");
    } else {
        this->max_time_grap = max_time_grap_ps;
    }
    this->outfilename = outfilename;
    boost::trim(this->outfilename);
    if (this->outfilename.empty()) {
        throw runtime_error("outfilename cannot empty");
    }
}

string RotAcf::description() {
    stringstream ss;
    string title_line = "------ " + title() + " ------";
    ss << title_line << "\n";
    ss << " vector            = " << vectorSelector->description() << "\n";
    ss << " P                 = " << LegendrePolynomial << "  [ " << LegendreStr.at(LegendrePolynomial) << " ] \n";
    ss << " time_increment_ps = " << time_increment_ps << " (ps)\n";
    ss << " max_time_grap_ps  = " << max_time_grap << " (ps)\n";
    ss << " outfilename       = " << outfilename << "\n";
    ss << string(title_line.size(), '-') << '\n';
    return ss.str();
}
