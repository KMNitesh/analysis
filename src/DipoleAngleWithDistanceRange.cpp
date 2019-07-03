//
// Created by xiamr on 6/30/19.
//

#include "DipoleAngleWithDistanceRange.hpp"
#include "frame.hpp"
#include "ThrowAssert.hpp"

using namespace std;


void DipoleAngleWithDistanceRange::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->ids1)) this->group1.insert(atom);
                      if (Atom::is_match(atom, this->ids2)) this->group2.insert(atom);
                  });
}

void DipoleAngleWithDistanceRange::process(std::shared_ptr<Frame> &frame) {
    for (auto &ref : group1) {
        for (auto &atom: group2) {
            auto mol = atom->molecule.lock();

            double xr = atom->x - ref->x;
            double yr = atom->y - ref->y;
            double zr = atom->z - ref->z;

            frame->image(xr, yr, zr);

            double distance = sqrt(xr * xr + yr * yr + zr * zr);

            if (cutoff1 <= distance and distance < cutoff2) {

                auto[dipole_x, dipole_y, dipole_z] = mol->calc_dipole(frame);

                double dipole_scalar = std::sqrt(dipole_x * dipole_x + dipole_y * dipole_y + dipole_z * dipole_z);

                double _cos = (xr * dipole_x + yr * dipole_y + zr * dipole_z) / (distance * dipole_scalar);

                double angle = acos(_cos) * radian;

                int i_angle_bin = int(angle / angle_width) + 1;

                if (i_angle_bin <= angle_bins) {
                    hist[i_angle_bin] += 1;
                }
            }
        }
    }
}

void DipoleAngleWithDistanceRange::print(std::ostream &os) {

    os << string(50, '#') << '\n';
    os << "# " << DipoleAngleWithDistanceRange::title() << '\n';
    os << "# Group1 > " << ids1 << '\n';
    os << "# Group2 > " << ids2 << '\n';
    os << "# angle_width(degree) > " << angle_width << '\n';
    os << "# Cutoff1(Ang) > " << cutoff1 << '\n';
    os << "# Cutoff2(Ang) > " << cutoff2 << '\n';
    os << string(50, '#') << '\n';
    os << format("#%15s %15s\n", "Angle(degree)", "Probability Density(% degree-1)");

    printData(os);

    os << string(50, '#') << '\n';
}

void DipoleAngleWithDistanceRange::printData(ostream &os) const {
    double total = 0.0;

    for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
        total += hist.at(i_angle);
    }

    for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
        os << format("%15.3f %15.3f\n", (i_angle - 0.5) * angle_width, 100 * (hist.at(i_angle) / total) / angle_width);
    }

}

void DipoleAngleWithDistanceRange::readInfo() {
    Atom::select2group(ids1, ids2);
    double angle_max = choose(0.0, 180.0, "Enter Maximum Angle to Accumulate[180.0 degree]:", true, 180.0);
    angle_width = choose(0.0, 180.0, "Enter Width of Angle Bins [0.5 degree]:", true, 0.5);

    cutoff1 = choose(0.0, 100.0, "Cutoff1 [Angstrom]:", false);
    cutoff2 = choose(0.0, 100.0, "Cutoff2 [Angstrom]:", false);

    throw_assert(cutoff1 < cutoff2, "Cutoff1 must less than Cutoff2");

    angle_bins = int(angle_max / angle_width);

    for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
        hist[i_angle] = 0;
    }

}
