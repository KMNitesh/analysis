//
// Created by xiamr on 7/3/19.
//

#include <functional>
#include <tbb/tbb.h>
#include "RotAcf.hpp"
#include "frame.hpp"
#include "molecule.hpp"
#include "ThrowAssert.hpp"
#include "VectorSelectorFactory.hpp"
#include "RotACFGrammar.hpp"


using namespace std;

void RotAcf::processFirstFrame(std::shared_ptr<Frame> &frame) {
    rots.resize(vectorSelector->initialize(frame));
}

void RotAcf::process(std::shared_ptr<Frame> &frame) {
    auto it2 = rots.begin();
    auto vectors = vectorSelector->calculateVectors(frame);

    for (auto it1 = vectors.begin(); it1 != vectors.end(); ++it1, ++it2) {
        it2->push_back(*it1);
    }
}

void RotAcf::print(std::ostream &os) {
    vector<double> acf;
    switch (LegendrePolynomial) {
        case 1:
            acf = calculate([](auto x) { return x; });
            break;
        case 2:
            acf = calculate([](auto x) { return 0.5 * (3 * x * x - 1); });
            break;
    }

    vector<double> integration = integrate(acf);

    os << "*********************************************************\n";

    os << " rotational autocorrelation function\n";
    vectorSelector->print(os);
    os << "Legendre Polynomial : ";
    switch (LegendrePolynomial) {
        case 1:
            os << "P1 = x\n";
            break;
        case 2:
            os << "P2 = (1/2)(3x^2 -1)\n";
            break;
    }
    os << "    Time Gap      ACF               integrate\n";
    os << "      (ps)                            (ps)\n";

    for (std::size_t t = 0; t < acf.size(); t++) {
        os << boost::format("%12.2f%18.14f%15.5f\n") % (t * time_increment_ps) % acf[t] % integration[t];
    }
    os << "*********************************************************\n";

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
class ParallelBody {
public:
    const vector<vector<tuple<double, double, double>>> &rots;
    double *acf = nullptr;
    unsigned long long *ntime = nullptr;

    size_t array_length;
    size_t max_time_grap_step;
    Function f;


    explicit ParallelBody(const vector<vector<tuple<double, double, double>>> &rots,
                          size_t max_time_grap_step, size_t array_length, Function f)
            : rots(rots), acf(new double[array_length]{}),
              ntime(new unsigned long long[array_length]{}),
              array_length(array_length),
              max_time_grap_step(max_time_grap_step), f(f) {}

    ParallelBody(const ParallelBody &rhs, tbb::split)
            : rots(rhs.rots), acf(new double[rhs.array_length]{}),
              ntime(new unsigned long long[rhs.array_length]{}),
              array_length(rhs.array_length),
              max_time_grap_step(rhs.max_time_grap_step), f(rhs.f) {}

    void join(const ParallelBody &rhs) {
        for (size_t i = 1; i < array_length; i++) {
            acf[i] += rhs.acf[i];
            ntime[i] += rhs.ntime[i];
        }
    }

    virtual ~ParallelBody() {
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

    ParallelBody parallelBody(rots, max_time_grap_step, min(rots[0].size(), max_time_grap_step + 1), f);

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
    std::cout << "1. P1 = x\n";
    std::cout << "2. P2 = (1/2)(3x^2 -1)\n";
    LegendrePolynomial = choose(1, 2, "select > ");
    this->time_increment_ps = choose(0.0, std::numeric_limits<double>::max(),
                                     "Enter the Time Increment in Picoseconds [0.1]:", true, 0.1);
    this->max_time_grap = choose(0.0, std::numeric_limits<double>::max(),
                                 "Enter the Max Time Grap in Picoseconds :");
}

void RotAcf::readAST(const RotAcfNode &ast) {
    if (!ast->vectorSelctor) {
        throw runtime_error("vector not vaild");
    } else {
        vectorSelector = VectorSelectorFactory::getVectorSelectorByAST(ast->vectorSelctor.value());
    }
    if (!ast->Legendre) {
        throw runtime_error("legendre polynomial is empty");
    }
    if (!(ast->Legendre.value() == 1 or ast->Legendre.value() == 2)) {
        throw runtime_error("legendre polynomial must be 1 or 2");
    }

    LegendrePolynomial = ast->Legendre.value();

    if (ast->time_increment_ps) {
        this->time_increment_ps = 0.1;
    } else if (ast->time_increment_ps.value() <= 0) {
        throw runtime_error("`time_increment_ps` must be postive");
    } else {
        this->time_increment_ps = ast->time_increment_ps.value();
    }

    if (!ast->max_time_grap_ps) {
        throw runtime_error("`max_time_grap_ps is empty");
    } else if (ast->max_time_grap_ps.value() <= 0) {
        throw runtime_error("`max_time_grap_ps` must be postive");
    } else if (ast->max_time_grap_ps.value() <= this->time_increment_ps) {
        throw runtime_error("`max_time_grap_ps` must be larger than `time_increment_ps`");
    } else {
        this->max_time_grap = ast->max_time_grap_ps.value();
    }
    outfilename = ast->outfilename;
    boost::trim(outfilename);
    if (outfilename.empty()) {
        throw runtime_error("outfilename cannot empty");
    }
}

