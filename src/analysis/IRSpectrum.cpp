//
// Created by xiamr on 8/11/19.
//

#include <tbb/tbb.h>
#include <boost/range/algorithm.hpp>
#include <boost/range/numeric.hpp>
#include "IRSpectrum.hpp"
#include "common.hpp"
#include "frame.hpp"
#include "molecule.hpp"

IRSpectrum::IRSpectrum() {
    enable_outfile = true;
    enable_tbb = true;
}

void IRSpectrum::process(std::shared_ptr<Frame> &frame) {
    dipole_evolution.emplace_back(getDipole(frame));
}

void IRSpectrum::print(std::ostream &os) {
    printData(os, dipole_evolution, time_increment_ps, selected_mols_mask);
}

void IRSpectrum::printData(std::ostream &os,
                           const std::deque<std::tuple<double, double, double>> &dipole_evolution,
                           double time_increment_ps, boost::optional<AmberMask> mask) {
    auto acf = calculateAcf(dipole_evolution);
    os << std::string(50, '#') << '\n';
    os << "# " << title() << '\n';
    os << "# time_increment_ps > " << time_increment_ps << '\n';
    if (mask) {
        os << "# AmberMask > " << mask.value() << '\n';
    }
    os << std::string(50, '#') << '\n';

    os << boost::format("%15s %15s\n") % "Time(ps)" % "ACF";

    for (std::size_t i = 0; i < acf.size(); ++i) {
        os << boost::format("%15.4f %15.5f\n") % (time_increment_ps * i) % acf[i];
    }

    os << std::string(50, '#') << '\n';

    std::vector<double> intense = calculateIntense(acf, time_increment_ps);

    os << boost::format("%15s %15s\n") % "Frequency (cm-1)" % "Intensity";

    for (std::size_t i = 0; i < intense.size(); ++i) {
        os << boost::format("%15.3f %15.5f\n") % (i + 1) % intense[i];
    }

    os << std::string(50, '#') << '\n';
}

std::vector<double> IRSpectrum::calculateIntense(const std::vector<double> &acf, double time_increment_ps) {
    constexpr double lightSpeed = 2.99792458E-2;

    const double factor = 2.0 * pi * lightSpeed;

    std::vector<double> intense(5000, 0);
    for (size_t f = 0; f < intense.size(); ++f) {
        auto freq = f + 1;
        for (size_t k = 0; k < acf.size(); ++k) {
            intense[f] += acf[k] * cos(freq * factor * time_increment_ps * k);
        }
    }

    auto max_intense = *std::max_element(intense.begin(), intense.end());

    std::for_each(intense.begin(), intense.end(), [max_intense](auto &i) { i /= max_intense; });
    return intense;
}


template<typename Container>
std::vector<double> IRSpectrum::calculateAcf(const Container &evolution) {

    auto max_calc_length = evolution.size() / 3;

    auto[acf, ntime] = tbb::parallel_reduce(
            tbb::blocked_range(std::size_t(0), evolution.size()),
            std::make_pair(std::vector<double>(max_calc_length), std::vector<long>(max_calc_length)),
            [&evolution, max_calc_length](const tbb::blocked_range<size_t> &range, auto init) {
                auto &[acf, ntime] = init;
                for (size_t i = range.begin(); i != range.end(); ++i) {
                    for (size_t j = i; j < std::min(evolution.size(), i + max_calc_length); ++j) {
                        auto n = j - i;
                        ++ntime[n];
                        acf[n] += dot_multiplication(evolution[i], evolution[j]);
                    }
                }
                return init;
            },
            [](const auto &lhs, const auto &rhs) {
                auto &[lhs_acf, lhs_ntime] = lhs;
                auto &[rhs_acf, rhs_ntime] = rhs;
                std::vector<double> acf(lhs_acf.size());
                std::vector<long> ntime(lhs_ntime.size());
                boost::transform(lhs_acf, rhs_acf, acf.begin(), std::plus<>());
                boost::transform(lhs_ntime, rhs_ntime, ntime.begin(), std::plus<>());
                return std::make_pair(acf, ntime);
            }
    );

    acf[0] /= ntime[0];

    for (size_t i = 1; i < acf.size(); ++i) {
        acf[i] /= (ntime[i] * acf[0]);
    }
    acf[0] = 1.0;
    return acf;
}

template std::vector<double> IRSpectrum::calculateAcf(const std::deque<double> &);

void IRSpectrum::readInfo() {
    time_increment_ps = choose(0.0, 100.0, "time_increment_ps [0.1 ps] :", Default(0.1));
    Atom::select1group(selected_mols_mask, " Enter molecule mask for dipole calculation > ");
}

void IRSpectrum::calculateSpectrum(const std::string &out) {
    auto time_increment_ps = choose(0.0, 100.0, "time_increment_ps [0.001 ps] :", Default(0.001));
    std::string file = choose_file("Enter Dipole Evolution Data File : ").isExist(true);

    std::deque<std::tuple<double, double, double>> dipole_evolution;

    std::ifstream ifstream(file);

    double dplx, dply, dplz;
    while (ifstream) {
        ifstream >> dplx >> dply >> dplz;
        dipole_evolution.emplace_back(dplx, dply, dplz);
    }

    std::ofstream ofstream(out);

    printData(ofstream, dipole_evolution, time_increment_ps);
}

void IRSpectrum::processFirstFrame(std::shared_ptr<Frame> &frame) {
    boost::for_each(frame->atom_list,
                    [this](std::shared_ptr<Atom> &atom) {
                        if (Atom::is_match(atom, selected_mols_mask)) {
                            selected_mols.insert(atom->molecule.lock());
                        }
                    });
}

std::tuple<double, double, double> IRSpectrum::getDipole(std::shared_ptr<Frame> &frame) {

    return boost::accumulate(selected_mols, std::tuple<double, double, double>{},
                             [&frame](auto &init, auto &mol) {
                                 return init + mol->calc_dipole(frame);
                             });
}
