//
// Created by xiamr on 6/30/19.
//

#include "DipoleAxisDistribution.hpp"
#include "frame.hpp"

using namespace std;

void DipoleAxisDistribution::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->ids)) this->group.insert(atom);
                  });
}

void DipoleAxisDistribution::process(std::shared_ptr<Frame> &frame) {
    for (auto &atom: group) {
        auto mol = atom->molecule.lock();

        auto[dipole_x, dipole_y, dipole_z] = mol->calc_dipole(frame);

        auto dipole_scalar = std::sqrt(dipole_x * dipole_x + dipole_y * dipole_y + dipole_z * dipole_z);

        auto _cos = (xr * dipole_x + yr * dipole_y + zr * dipole_z) / dipole_scalar;

        auto angle = acos(_cos) * radian;

        hist.update(angle);
    }
}

void DipoleAxisDistribution::print(std::ostream &os) {
    os << string(50, '#') << '\n';
    os << "# " << DipoleAxisDistribution::title() << '\n';
    os << "# Group > " << ids << '\n';
    os << "# angle_width(degree) > " << hist.getWidth() << '\n';
    os << "# " << (xr == 1.0 ? "X-Axis" : (yr == 1.0 ? "Y-Axis" : "Z-Axis")) << '\n';
    os << string(50, '#') << '\n';
    os << format("#%15s %15s\n", "Angle(degree)", "Probability Density(% degree-1)");

    printData(os);

    os << string(50, '#') << '\n';
}

void DipoleAxisDistribution::readInfo() {
    Atom::select1group(ids);
    double angle_max = choose(0.0, 180.0, "Enter Maximum Angle to Accumulate[180.0 degree]:", true, 180.0);
    auto angle_width = choose(0.0, 180.0, "Enter Width of Angle Bins [0.5 degree]:", true, 0.5);

    auto menu = [] {
        cout << "Axis selection \n";
        cout << "  1. X-Axis\n";
        cout << "  2. Y-Axis\n";
        cout << "  3. Z-Axis\n";

        return choose(1, 3, " select > ");
    };

    switch (menu()) {
        case 1:
            xr = 1.0;
            break;
        case 2:
            yr = 1.0;
            break;
        case 3:
            zr = 1.0;
            break;
    }

    hist.initialize(angle_max, angle_width);
}

void DipoleAxisDistribution::printData(std::ostream &os) const {
    for (auto[grid, value] : hist.getDistribution()) {
        os << format("%15.3f %15.3f\n", grid, 100 * value);
    }
}