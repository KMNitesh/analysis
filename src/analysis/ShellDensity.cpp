//
// Created by xiamr on 6/14/19.
//

#include "ShellDensity.hpp"

#include "frame.hpp"

using namespace std;

void ShellDensity::process(std::shared_ptr<Frame> &frame) {
    nframe++;
    for (auto &ref :group1) {
        for (auto &atom : group2) {
            int ibin = int(atom_distance(ref, atom, frame) / distance_width) + 1;
            if (ibin <= distance_bins) {
                hist[ibin]++;
            }
        }

    }
}

void ShellDensity::print(std::ostream &os) {
    os << "************************************************" << endl;
    os << "***** Shell Density Function ****" << endl;

    os << "First Type : " << ids1 << " Second Type : " << ids2 << endl;

    os << "************************************************" << endl;
    os << "Bin    Distance    Densitry (count / Ang3 / frame)" << endl;

    for (int i = 1; i <= distance_bins; i++) {
        double dv = (4.0 / 3.0) * M_PI * (pow(i * distance_width, 3) - pow((i - 1) * distance_width, 3));
        os << boost::format("%d      %.4f      %g \n") %
              i % ((i - 0.5) * distance_width) % (hist[i] / (nframe * dv));
    }

    os << "************************************************" << endl;
}


void ShellDensity::readInfo() {
    Atom::select2group(ids1, ids2);
    double rmax = choose(0.0, std::numeric_limits<double>::max(), "Enter Maximum Distance to Accumulate[10.0 Ang]:",
                         Default(10.0));
    distance_width = choose(0.0, std::numeric_limits<double>::max(), "Enter Width of Distance Bins [0.01 Ang]:",
                            Default(0.01));
    distance_bins = int(rmax / distance_width);
    for (int i = 1; i <= distance_bins; i++) {
        hist[i] = 0;
    }
}

void ShellDensity::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->ids1)) this->group1.insert(atom);
                      if (Atom::is_match(atom, this->ids2)) this->group2.insert(atom);
                  });
}
