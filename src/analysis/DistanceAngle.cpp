//
// Created by xiamr on 7/4/19.
//

#include "DistanceAngle.hpp"
#include "frame.hpp"
#include "molecule.hpp"
#include "ThrowAssert.hpp"


using namespace std;

void DistanceAngle::processFirstFrame(std::shared_ptr<Frame> &frame) {
    for (auto &mol : frame->molecule_list) {
        shared_ptr<Atom> atom1, atom2;
        for (auto &atom : mol->atom_list) {
            if (Atom::is_match(atom, id1)) {
                atom1 = atom;
            } else if (Atom::is_match(atom, id2)) {
                atom2 = atom;
            } else if (Atom::is_match(atom, id3)) {
                group3.insert(atom);
            }
        }

        throw_assert((atom1 && atom2) || (!atom1 && !atom2), "Atom selection semtatic error");

        if (atom1 && atom2) {
            pairs.emplace_back(atom1, atom2);
        }
    }
    throw_assert(!pairs.empty(), "Atom selection semtatic error ! not atom1 & atom2 selected!");
}

void DistanceAngle::process(std::shared_ptr<Frame> &frame) {
    for (auto &[ref, atom2] : pairs) {

        auto r = atom2->getCoordinate() - ref->getCoordinate();
        frame->image(r);

        r /= vector_norm(r);

        for (auto &atom3 : group3) {

            auto v = atom3->getCoordinate() - ref->getCoordinate();
            frame->image(v);

            auto dist = vector_norm(v);

            v /= dist;

            auto _cos = dot_multiplication(r, v);

            auto angle = std::abs(acos(_cos) * radian - 90.0);

            int i_distance_bin = int(dist / distance_width) + 1;
            int i_angle_bin = int(angle / angle_width) + 1;

            if (i_distance_bin <= distance_bins and i_angle_bin <= angle_bins) {
                hist[{i_distance_bin, i_angle_bin}] += 1;
            }
        }

    }
}

void DistanceAngle::print(std::ostream &os) {
    os << string(50, '#') << '\n';
    os << "# " << DistanceAngle::title() << '\n';
    os << "# Group1 > " << id1 << '\n';
    os << "# Group2 > " << id2 << '\n';
    os << "# Group3 > " << id3 << '\n';
    os << "# distance_width(Ang) > " << distance_width << '\n';
    os << "# angle_width(degree) > " << angle_width << '\n';
    os << "# Temperature(K) > " << temperature << '\n';
    os << string(50, '#') << '\n';
    os << format("#%15s %15s %15s\n", "Distance(Ang)", "Angle(degree)", "Energy(kcal/mol)");

    printData(os);

    os << string(50, '#') << '\n';
}

void DistanceAngle::printData(ostream &os) const {
    const double factor = -kb * temperature * avogadro_constant / 4184.0;
    double max_value = 0.0;

    for (int i_distance = 1; i_distance < distance_bins; i_distance++) {
        double dv = pow(i_distance * distance_width, 3) - pow((i_distance - 1) * distance_width, 3);
        for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
            max_value = max(max_value, hist.at(make_pair(i_distance, i_angle)) / (dv));
        }
    }
    for (int i_distance = 1; i_distance < distance_bins; i_distance++) {
        double dv = pow(i_distance * distance_width, 3) - pow((i_distance - 1) * distance_width, 3);
        for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
            double pop = double(hist.at(make_pair(i_distance, i_angle))) / (max_value * dv);
            os << boost::format("%15.3f %15.3f %15.6f\n")
                  % ((i_distance - 0.5) * distance_width)
                  % ((i_angle - 0.5) * angle_width)
                  % (pop == 0.0 ? 100.0 : factor * log(pop));
        }
    }
}

void DistanceAngle::readInfo() {
    Atom::select1group(id1, "Please Enter mask for Atom1(An) > ");
    Atom::select1group(id2, "Please Enter mask for Atom2(OAn) > ");
    Atom::select1group(id3, "Please Enter mask for Atom3(Ow) > ");

    double rmax = choose(0.0, std::numeric_limits<double>::max(), "Enter Maximum Distance to Accumulate[10.0 Ang]:",
                         Default(10.0));
    distance_width = choose(0.0, std::numeric_limits<double>::max(), "Enter Width of Distance Bins [0.01 Ang]:",
                            Default(0.01));
    double angle_max = choose(0.0, 180.0, "Enter Maximum Angle to Accumulate[180.0 degree]:", Default(180.0));
    angle_width = choose(0.0, 180.0, "Enter Width of Angle Bins [0.5 degree]:", Default(0.5));
    temperature = choose(0.0, 10000.0, "Temperature [298] (K):", Default(298.0));

    distance_bins = int(rmax / distance_width);
    angle_bins = int(angle_max / angle_width);

    for (int i_distance = 1; i_distance <= distance_bins; i_distance++) {
        for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
            hist[make_pair(i_distance, i_angle)] = 0;
        }
    }
}
