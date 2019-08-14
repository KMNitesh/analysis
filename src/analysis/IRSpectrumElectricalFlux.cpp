//
// Created by xiamr on 8/13/19.
//

#include "IRSpectrumElectricalFlux.hpp"
#include "common.hpp"
#include "frame.hpp"
#include "VelocityAutocorrelationFunction.hpp"
#include "IRSpectrum.hpp"

IRSpectrumElectricalFlux::IRSpectrumElectricalFlux() {
    enable_outfile = true;
    enable_read_velocity = true;
}

void IRSpectrumElectricalFlux::process(std::shared_ptr<Frame> &frame) {
    auto flux = std::accumulate(frame->atom_list.begin(), frame->atom_list.end(),
                                std::tuple<double, double, double>{}, [](auto flux, auto &atom) {
                return flux + std::make_tuple(atom->vx, atom->vy, atom->vz) * atom->charge.value();
            });
    electricalFlux.emplace_back(flux);
}

void IRSpectrumElectricalFlux::print(std::ostream &os) {
    auto acf = VelocityAutocorrelationFunction::calculateAcf({electricalFlux}, std::numeric_limits<int>::max());
    os << std::string(50, '#') << '\n';
    os << "# " << title() << '\n';
    os << "# time_increment_ps > " << time_increment_ps << '\n';
    os << std::string(50, '#') << '\n';

    os << boost::format("%15s %15s\n") % "Time(ps)" % "ACF";

    for (std::size_t i = 0; i < acf.size(); ++i) {
        os << boost::format("%15.4f %15.5f\n") % (time_increment_ps * i) % acf[i];
    }

    os << std::string(50, '#') << '\n';

    acf.resize(acf.size() / 3);

    auto intense = IRSpectrum::calculateIntense(acf, time_increment_ps);

    os << boost::format("%15s %15s\n") % "Frequency (cm-1)" % "Intensity";

    for (std::size_t i = 0; i < intense.size(); ++i) {
        os << boost::format("%15.3f %15.5f\n") % (i + 400) % intense[i];
    }

    os << std::string(50, '#') << '\n';

}

void IRSpectrumElectricalFlux::readInfo() {
    time_increment_ps = choose(0.0, 100.0, "time_increment_ps [0.1 ps] :", true, 0.1);
}
