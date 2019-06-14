//
// Created by xiamr on 6/14/19.
//

#include "DipoleAngle.hpp"

#include "frame.hpp"
#include "atom.hpp"
#include "molecule.hpp"

using namespace std;

void DipoleAngle::process(std::shared_ptr<Frame> &frame) {

    shared_ptr<Atom> ref;

    for (auto &mol : frame->molecule_list) {
        mol->bExculde = true;
    }

    for (auto &atom : group1) {
        ref = atom;
        if (!atom->molecule.expired()) atom->molecule.lock()->calc_mass();
    }

    for (auto &atom : group2) {
        if (!atom->molecule.expired()) {
            auto mol = atom->molecule.lock();
            mol->ow = atom;
            mol->calc_mass();
            mol->bExculde = false;
        }
    }

    if (!ref) {
        std::cerr << "reference atom not found" << std::endl;
        exit(5);
    }
    double ref_x = ref->x;
    double ref_y = ref->y;
    double ref_z = ref->z;
    for (auto &mol: frame->molecule_list) {
        if (!mol->bExculde) {

            double x1 = mol->ow->x;
            double y1 = mol->ow->y;
            double z1 = mol->ow->z;

            double xr = x1 - ref_x;
            double yr = y1 - ref_y;
            double zr = z1 - ref_z;

            frame->image(xr, yr, zr);

            double distance = sqrt(xr * xr + yr * yr + zr * zr);

            auto dipole = mol->calc_dipole(frame);

            double dipole_scalar = sqrt(pow(get<0>(dipole), 2) + pow(get<1>(dipole), 2) + pow(get<2>(dipole), 2));

            double _cos =
                    (xr * get<0>(dipole) + yr * get<1>(dipole) + zr * get<2>(dipole)) / (distance * dipole_scalar);

            double angle = acos(_cos) * 180.0 / 3.1415926;

            int i_distance_bin = int(distance / distance_width) + 1;
            int i_angle_bin = int(angle / angle_width) + 1;

            if (i_distance_bin <= distance_bins and i_angle_bin <= angle_bins) {
                hist[make_pair(i_distance_bin, i_angle_bin)] += 1;
            }
        }
    }
}

void DipoleAngle::print() {

//    outfile << boost::format("%5s") % "Dist";
//    for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
//        outfile << boost::format("%7.1f") % (i_angle  * angle_width);
//    }
//    outfile << endl;
//    for (int i_distance = 1; i_distance < distance_bins; i_distance++) {
//        outfile << boost::format("%5.3f") % (i_distance * distance_width);
//        size_t total = 0;
//        for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
//            total += hist[make_pair(i_distance,i_angle)];
//        }
//        for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
//            outfile << boost::format("%7.4f") % ( total == 0 ? 0.0 : double(hist[make_pair(i_distance,i_angle)]) / total);
//        }
//        outfile << endl;
//    }
    for (int i_distance = 1; i_distance < distance_bins; i_distance++) {
        size_t total = 0;
        for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
            total += hist[make_pair(i_distance, i_angle)];
        }
        for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
            outfile << boost::format("%5.3f%10.3f%10.4f\n")
                       % ((i_distance - 0.5) * distance_width)
                       % ((i_angle - 0.5) * angle_width)
                       % (total == 0 ? 0.0 : double(hist[make_pair(i_distance, i_angle)]) / total);
        }
    }
}

void DipoleAngle::readInfo() {
    Atom::select2group(ids1, ids2);
    double rmax = choose(0.0, std::numeric_limits<double>::max(), "Enter Maximum Distance to Accumulate[10.0 Ang]:",
                         true, 10.0);
    distance_width = choose(0.0, std::numeric_limits<double>::max(), "Enter Width of Distance Bins [0.01 Ang]:", true,
                            0.01);
    double angle_max = choose(0.0, 180.0, "Enter Maximum Angle to Accumulate[180.0 degree]:", true, 180.0);
    angle_width = choose(0.0, 180.0, "Enter Width of Angle Bins [0.5 degree]:", true, 0.5);

    distance_bins = int(rmax / distance_width);
    angle_bins = int(angle_max / angle_width);

    for (int i_distance = 1; i_distance <= distance_bins; i_distance++) {
        for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
            hist[make_pair(i_distance, i_angle)] = 0;
        }
    }
}

void DipoleAngle::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->ids1)) this->group1.insert(atom);
                      if (Atom::is_match(atom, this->ids2)) this->group2.insert(atom);
                  });
}
