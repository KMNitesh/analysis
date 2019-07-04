//
// Created by xiamr on 7/3/19.
//

#include <tbb/tbb.h>
#include "RotAcf.hpp"
#include "frame.hpp"
#include "molecule.hpp"
#include "ThrowAssert.hpp"


using namespace std;

void RotAcf::processFirstFrame(std::shared_ptr<Frame> &frame) {
    for (auto &mol : frame->molecule_list) {
        shared_ptr<Atom> atom1, atom2, atom3;
        for (auto &atom : mol->atom_list) {
            if (Atom::is_match(atom, ids1)) {
                atom1 = atom;
            } else if (Atom::is_match(atom, ids2)) {
                atom2 = atom;
            } else if (Atom::is_match(atom, ids3)) {
                atom3 = atom;
            }
        }


        throw_assert((atom1 && atom2 && atom3) or (!atom1 && !atom2 && !atom3), "Atom selection semantic error");
        if (atom1 && atom2 && atom3) {
            pairs.emplace_back(atom1, atom2, atom3);
        }
    }

    throw_assert(!pairs.empty(), "Can not empty");
    rots.resize(pairs.size());
}

void RotAcf::process(std::shared_ptr<Frame> &frame) {
    auto it2 = rots.begin();
    for (auto it1 = pairs.begin(); it1 != pairs.end(); ++it1, ++it2) {
        it2->push_back(calVector(*it1, frame));
    }
}

void RotAcf::print(std::ostream &os) {

    vector<double> acf = calculate();

    // intergrate;

    vector<double> integration = integrate(acf);

    os << "*********************************************************\n";
    os << "Group1 > " << ids1 << " Group2 > " << ids2 << " Group3 > " << ids3 << '\n';
    os << " rotational autocorrelation function\n";

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

    tbb::parallel_reduce(tbb::blocked_range<int>(0, rots.size() - 1), parallelBody, tbb::auto_partitioner());

    for (size_t i = 1; i < parallelBody.acf.size(); i++) {
        assert(parallelBody.ntime[i] > 0);
        parallelBody.acf[i] /= parallelBody.ntime[i];
    }

    parallelBody.acf[0] = 1.0;
    return parallelBody.acf;
}

void RotAcf::readInfo() {
    Atom::select1group(ids1, "Please Enter for Atom1 > ");
    Atom::select1group(ids2, "Please Enter for Atom2 > ");
    Atom::select1group(ids3, "Please Enter for Atom3 > ");

    this->time_increment_ps = choose(0.0, std::numeric_limits<double>::max(),
                                     "Enter the Time Increment in Picoseconds [0.1]:", true, 0.1);
}

std::tuple<double, double, double> RotAcf::calVector(
        std::tuple<std::shared_ptr<Atom>, std::shared_ptr<Atom>, std::shared_ptr<Atom>> &atoms,
        std::shared_ptr<Frame> &frame) {

    auto &[atom_i, atom_j, atom_k] = atoms;

    auto u1 = atom_i->x - atom_j->x;
    auto u2 = atom_i->y - atom_j->y;
    auto u3 = atom_i->z - atom_j->z;

    auto v1 = atom_k->x - atom_j->x;
    auto v2 = atom_k->y - atom_j->y;
    auto v3 = atom_k->z - atom_j->z;

    frame->image(u1, u2, u3);
    frame->image(v1, v2, v3);

    auto xv3 = u2 * v3 - u3 * v2;
    auto yv3 = u3 * v1 - u1 * v3;
    auto zv3 = u1 * v2 - u2 * v1;

    auto len = std::sqrt(xv3 * xv3 + yv3 * yv3 + zv3 * zv3);

    assert(len > 0);

    return {xv3 / len, yv3 / len, zv3 / len};
}
