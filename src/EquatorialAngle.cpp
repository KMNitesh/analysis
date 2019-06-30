//
// Created by xiamr on 6/30/19.
//

#include "EquatorialAngle.hpp"
#include "frame.hpp"
#include "ThrowAssert.hpp"

using namespace std;

void EquatorialAngle::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->ids1)) this->group1.insert(atom);
                      if (Atom::is_match(atom, this->ids2)) this->group2.insert(atom);
                      if (Atom::is_match(atom, this->ids3)) this->group3.insert(atom);
                  });
}

void EquatorialAngle::process(std::shared_ptr<Frame> &frame) {
    for (auto &vec1_atom1 : group1) {
        for (auto &vec1_atom2: group2) {
            if (vec1_atom1->molecule.lock() == vec1_atom2->molecule.lock()) {

                double xr = vec1_atom2->x - vec1_atom1->x;
                double yr = vec1_atom2->y - vec1_atom1->y;
                double zr = vec1_atom2->z - vec1_atom1->z;

                frame->image(xr, yr, zr);

                double leng1 = sqrt(xr * xr + yr * yr + zr * zr);

                for (auto &atom3 : group3) {

                    double xr2 = atom3->x - vec1_atom1->x;
                    double yr2 = atom3->y - vec1_atom1->y;
                    double zr2 = atom3->z - vec1_atom1->z;

                    frame->image(xr2, yr2, zr2);

                    double distance = sqrt(xr2 * xr2 + yr2 * yr2 + zr2 * zr2);

                    if (cutoff1 <= distance and distance < cutoff2) {

                        double _cos = (xr * xr2 + yr * yr2 + zr * zr2) / (distance * leng1);

                        double angle = acos(_cos) * radian;

                        int i_angle_bin = int(angle / angle_width) + 1;

                        if (i_angle_bin <= angle_bins) {
                            hist[i_angle_bin] += 1;
                        }
                    }
                }
            }
        }
    }
}

void EquatorialAngle::print(std::ostream &os) {
    os << string(50, '#') << '\n';
    os << "# " << EquatorialAngle::title() << '\n';
    os << "# Group1 > " << ids1 << '\n';
    os << "# Group2 > " << ids2 << '\n';
    os << "# Group3 > " << ids3 << '\n';
    os << "# angle_width > " << angle_width << '\n';
    os << string(50, '#') << '\n';
    os << format("#%15s %15s\n", "Angle(degree)", "Probability(%)");

    printData(os);

    os << string(50, '#') << '\n';

}

void EquatorialAngle::readInfo() {
    Atom::select1group(ids1, "Enter mask for atom1 : ");
    Atom::select1group(ids2, "Enter mask for atom2(in same mol with atom1) : ");
    Atom::select1group(ids3, "Enter mask for atom3 : ");

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

void EquatorialAngle::printData(std::ostream &os) const {
    double total = 0.0;

    for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
        total += hist.at(i_angle);
    }

    for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
        os << format("%15.3f %15.3f\n", (i_angle - 0.5) * angle_width, 100 * (hist.at(i_angle) / total) / angle_width);
    }
}
