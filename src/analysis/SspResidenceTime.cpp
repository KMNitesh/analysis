//
// Created by xiamr on 9/16/19.
//
#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm.hpp>
#include <tbb/tbb.h>
#include "SspResidenceTime.hpp"
#include "common.hpp"

void SspResidenceTime::print(std::ostream &os) {
    calculateSSP();

    os << "# " << title() << '\n';
    os << "# center_atom_mask: " << center_atom_mask << ",  Ow_atom_mask: " << Ow_atom_mask << '\n';
    os << "# dist_cutoff (Ang) = " << dis_cutoff << '\n';

    os << format("#%15s %15s\n", "Frame", "SSP R(t)");
    for (const auto &ele : Rt_array | boost::adaptors::indexed(1)) {
        os << format(" %15d %15.6f\n", ele.index(), 1 - ele.value());
    }

}

void SspResidenceTime::calculateSSP() {
    fill_mark();

    std::vector<std::vector<int>> hydrationed_atoms(steps);

    for (int step = 0; step < steps; ++step) {
        for (int atom = 0; atom < atom_num; atom++) {
            if (mark(step, atom)) {
                hydrationed_atoms[step].emplace_back(atom);
            }
        }
    }

    auto[Rt, ntime] = tbb::parallel_reduce(
            tbb::blocked_range<int>(0, steps - 1),
            std::make_pair(std::vector<int, tbb::tbb_allocator<int>>(steps - 1),
                           std::vector<int, tbb::tbb_allocator<int>>(steps - 1)),
            [this, &hydrationed_atoms](const tbb::blocked_range<int> &r, auto init) {
                auto &[Rt_array, ntime]  = init;
                for (auto i = r.begin(); i != r.end(); ++i) {
                    for (auto j = i + 1; j < steps; ++j) {
                        auto n = j - i - 1;
                        for (auto atom : hydrationed_atoms[i]) {
                            ++ntime[n];
                            if (!mark(j, atom)) ++Rt_array[n];
                        }
                    }
                }
                return init;
            }, [](const auto &lhs, const auto &rhs) {
                std::vector<int, tbb::tbb_allocator<int>> Rt;
                std::vector<int, tbb::tbb_allocator<int>> ntime;

                Rt.reserve(std::get<0>(lhs).size());
                ntime.reserve(std::get<0>(lhs).size());

                boost::transform(std::get<0>(lhs), std::get<0>(rhs), std::back_inserter(Rt), std::plus<>());
                boost::transform(std::get<1>(lhs), std::get<1>(rhs), std::back_inserter(ntime), std::plus<>());
                return std::make_pair(Rt, ntime);
            });

    Rt_array.reserve(steps - 1);
    boost::transform(Rt, ntime, std::back_inserter(Rt_array), std::divides<double>());
}

void SspResidenceTime::readInfo() {
    readBasicSetting();
}
