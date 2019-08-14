//
// Created by xiamr on 8/13/19.
//

#include <boost/range/combine.hpp>

#include "VelocityAutocorrelationFunction.hpp"
#include "common.hpp"
#include "frame.hpp"

VelocityAutocorrelationFunction::VelocityAutocorrelationFunction() {
    enable_outfile = true;
    enable_read_velocity = true;
}

void VelocityAutocorrelationFunction::processFirstFrame(std::shared_ptr<Frame> &frame) {
    velocities.resize(frame->atom_list.size());
}


void VelocityAutocorrelationFunction::process(std::shared_ptr<Frame> &frame) {
    auto it1 = velocities.begin();
    auto it2 = frame->atom_list.begin();
    for (; it1 != velocities.end(); ++it1, ++it2) {
        it1->emplace_back((*it2)->vx, (*it2)->vy, (*it2)->vz);
    }
}

void VelocityAutocorrelationFunction::print(std::ostream &os) {
    acf = calculateAcf(velocities, ceil(max_time_grap_ps / time_increment_ps));

    os << std::string(50, '#') << '\n';
    os << "# " << title() << '\n';
    os << "# time_increment_ps (ps) > " << time_increment_ps << '\n';
    os << "# max_time_grap_ps  (ps) > " << max_time_grap_ps << '\n';
    os << std::string(50, '#') << '\n';

    os << boost::format("%15s %15s\n") % "Time(ps)" % "ACF";

    for (std::size_t i = 0; i < acf.size(); ++i) {
        os << boost::format("%15.4f %15.5f\n") % (time_increment_ps * i) % acf[i];
    }

    os << std::string(50, '#') << '\n';
}

std::vector<long double>
VelocityAutocorrelationFunction::calculateAcf(
        const std::vector<std::deque<std::tuple<double, double, double>>> &velocities, int max_time_grap_frame) {
    std::vector<long double> acf(std::min<std::size_t>(velocities.at(0).size(), max_time_grap_frame + 1), 0);
    std::vector<long> ntime(std::min<std::size_t>(velocities.at(0).size(), max_time_grap_frame + 1), 0);

    for (auto &vel : velocities) {
        for (size_t i = 0; i < vel.size(); ++i) {
            for (size_t j = i; j < std::min<std::size_t>(vel.size(), max_time_grap_frame + i + 1); ++j) {
                auto n = j - i;
                ++ntime[n];
                acf[n] += dot_multiplication(vel[i], vel[j]);
            }
        }
    }
    acf[0] /= ntime[0];

    for (size_t i = 1; i < acf.size(); ++i) {
        acf[i] /= (ntime[i] * acf[0]);
    }
    acf[0] = 1.0;
    return acf;
}

void VelocityAutocorrelationFunction::readInfo() {
    time_increment_ps = choose(0.0, 100.0, "time_increment_ps [0.1 ps] :", true, 0.1);
    max_time_grap_ps = choose(0.0, 10000.0, "max_time_grap_ps [100 ps] :", true, 100.0);
}

