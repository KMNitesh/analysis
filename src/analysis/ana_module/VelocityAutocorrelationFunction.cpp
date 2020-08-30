//
// Created by xiamr on 8/13/19.
//

#include <boost/range/algorithm.hpp>
#include <boost/range/combine.hpp>

#include <tbb/tbb.h>

#include "VelocityAutocorrelationFunction.hpp"

#include "data_structure/frame.hpp"
#include "utils/common.hpp"

VelocityAutocorrelationFunction::VelocityAutocorrelationFunction() {
    enable_outfile = true;
    enable_read_velocity = true;
    enable_tbb = true;
}

void VelocityAutocorrelationFunction::processFirstFrame(std::shared_ptr<Frame> &frame) {
    velocities.resize(frame->atom_list.size());
}

void VelocityAutocorrelationFunction::process(std::shared_ptr<Frame> &frame) {
    if (!frame->has_velocity) {
        std::cerr << "Trajectory does not have velocity\n";
        std::exit(EXIT_FAILURE);
    }

    auto it1 = velocities.begin();
    auto it2 = frame->atom_list.begin();
    for (; it1 != velocities.end(); ++it1, ++it2) {
        it1->emplace_back((*it2)->getVelocities());
    }
}

void VelocityAutocorrelationFunction::print(std::ostream &os) {
    acf = calculateAcf(velocities, std::ceil(max_time_gap_ps / time_increment_ps));

    os << std::string(50, '#') << '\n';
    os << "# " << title() << '\n';
    os << "# time_increment_ps (ps) > " << time_increment_ps << '\n';
    os << "# max_time_gap_ps  (ps) > " << max_time_gap_ps << '\n';
    os << std::string(50, '#') << '\n';

    os << boost::format("%15s %15s\n") % "Time(ps)" % "ACF";

    for (std::size_t i = 0; i < acf.size(); ++i) {
        os << boost::format("%15.4f %15.5f\n") % (time_increment_ps * i) % acf[i];
    }

    os << std::string(50, '#') << '\n';
}

std::vector<double> VelocityAutocorrelationFunction::calculateAcf(
    const std::vector<std::deque<std::tuple<double, double, double>>> &velocities, int max_time_gap_frame) {
    const auto max_calculate_length = std::min<std::size_t>(velocities.at(0).size(), max_time_gap_frame + 1);

    auto acf = tbb::parallel_reduce(
        tbb::blocked_range2d<size_t>(0, velocities.size(), 0, velocities.front().size()),
        std::vector<double>(max_calculate_length),
        [&velocities, max_time_gap_frame](const tbb::blocked_range2d<size_t> &range, auto acf) {
            for (auto index = range.rows().begin(); index != range.rows().end(); ++index) {
                auto &vel = velocities[index];
                for (auto i = range.cols().begin(); i != range.cols().end(); ++i) {
                    for (auto j = i; j < std::min<std::size_t>(vel.size(), max_time_gap_frame + i + 1); ++j) {
                        acf[j - i] += dot_multiplication(vel[i], vel[j]);
                    }
                }
            }
            return acf;
        },
        [](const auto &lhs, const auto &rhs) {
            std::vector<double> acf(lhs.size());
            boost::transform(lhs, rhs, acf.begin(), std::plus<>());
            return acf;
        });
    auto total_mols = velocities.size();
    auto total_size = velocities[0].size();
    acf[0] /= total_mols * total_size;
    for (size_t i = 1; i < acf.size(); ++i) {
        acf[i] /= total_mols * (total_size - i) * acf[0];
    }
    acf[0] = 1.0;
    return acf;
}

void VelocityAutocorrelationFunction::readInfo() {
    time_increment_ps = choose(100.0, "time_increment_ps [0.1 ps] :", Default(0.1));
    max_time_gap_ps = choose(10000.0, "max_time_gap_ps [100 ps] :", Default(100.0));
}

void VelocityAutocorrelationFunction::setParameters(double time_increment_ps, double max_time_gap_ps,
                                                    const std::string &outfilename) {
    this->time_increment_ps = time_increment_ps;
    this->max_time_gap_ps = max_time_gap_ps;
    setOutFilename(outfilename);
}