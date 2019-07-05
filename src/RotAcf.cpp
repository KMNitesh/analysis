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


using namespace std;

void RotAcf::processFirstFrame(std::shared_ptr<Frame> &frame) {
    rots.resize(vectorSelector->initialize(frame));
}

void RotAcf::process(std::shared_ptr<Frame> &frame) {
    auto it2 = rots.begin();
    auto vectors = vectorSelector->calcaulteVectors(frame);

    for (auto it1 = vectors.begin(); it1 != vectors.end(); ++it1, ++it2) {
        it2->push_back(*it1);
    }
}

void RotAcf::print(std::ostream &os) {

    vector<double> acf = calculate();

    vector<double> integration = integrate(acf);

    os << "*********************************************************\n";

    os << " rotational autocorrelation function\n";
    vectorSelector->print(os);
    os << "    Time Gap      ACF               intergrate\n";
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

vector<double> RotAcf::calculate() const {

    class ParallelBody {
    public:
        const vector<vector<tuple<double, double, double>>> &rots;
        vector<double> acf;
        vector<int> ntime;

        explicit ParallelBody(const vector<vector<tuple<double, double, double>>> &rots)
                : rots(rots), acf(rots[0].size(), 0.0), ntime(rots[0].size(), 0) {}

        ParallelBody(const ParallelBody &body, tbb::split)
                : rots(body.rots), acf(body.rots[0].size(), 0.0), ntime(body.rots[0].size(), 0) {}

        void join(const ParallelBody &body) {
            for (size_t i = 1; i < acf.size(); i++) {
                acf[i] += body.acf[i];
                ntime[i] += body.ntime[i];
            }
        }

        void operator()(const tbb::blocked_range<int> &range) {
            for (int index = range.begin(); index != range.end(); index++) {
                auto &_vector = rots[index];
                auto total_size = _vector.size();

                for (size_t i = 0; i < total_size - 1; i++) {
                    for (size_t j = i + 1; j < total_size; j++) {
                        auto m = j - i;

                        assert(i < _vector.size());
                        assert(j < _vector.size());

                        auto[xr1, yr1, zr1] = _vector[i];
                        auto[xr2, yr2, zr2] = _vector[j];

                        double cos = xr1 * xr2 + yr1 * yr2 + zr1 * zr2;

                        acf[m] += cos;
                        ntime[m]++;
                    }
                }
            }
        }
    } parallelBody(rots);

    tbb::parallel_reduce(tbb::blocked_range<int>(0, rots.size()), parallelBody, tbb::auto_partitioner());

    for (size_t i = 1; i < parallelBody.acf.size(); i++) {
        assert(parallelBody.ntime[i] > 0);
        parallelBody.acf[i] /= parallelBody.ntime[i];
    }

    parallelBody.acf[0] = 1.0;
    return parallelBody.acf;
}

void RotAcf::readInfo() {

    vectorSelector = VectorSelectorFactory::getVectorSelector();
    vectorSelector->readInfo();
    this->time_increment_ps = choose(0.0, std::numeric_limits<double>::max(),
                                     "Enter the Time Increment in Picoseconds [0.1]:", true, 0.1);
}

