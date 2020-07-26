
#include "DistanceAngle.hpp"

#include "data_structure/frame.hpp"
#include "data_structure/molecule.hpp"
#include "utils/ThrowAssert.hpp"
#include "utils/common.hpp"

DistanceAngle::DistanceAngle() { enable_outfile = true; }

void DistanceAngle::processFirstFrame(std::shared_ptr<Frame> &frame) {
    for (auto &mol : frame->molecule_list) {
        std::shared_ptr<Atom> atom1, atom2;
        for (auto &atom : mol->atom_list) {
            if (is_match(atom, mask1)) {
                atom1 = atom;
            } else if (is_match(atom, mask2)) {
                atom2 = atom;
            } else if (is_match(atom, mask3)) {
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
    os << std::string(50, '#') << '\n';
    os << "# " << title() << '\n';
    os << "# Group1 > " << mask1 << '\n';
    os << "# Group2 > " << mask2 << '\n';
    os << "# Group3 > " << mask3 << '\n';
    os << "# distance_width(Ang) > " << distance_width << '\n';
    os << "# angle_width(degree) > " << angle_width << '\n';
    os << "# Temperature(K) > " << temperature << '\n';
    os << std::string(50, '#') << '\n';
    os << format("#%15s %15s %15s\n", "Distance(Ang)", "Angle(degree)", "Energy(kcal/mol)");

    printData(os);

    os << std::string(50, '#') << '\n';
}

void DistanceAngle::printData(std::ostream &os) const {
    const double factor = -kb * temperature * avogadro_constant / 4184.0;
    double max_value = 0.0;

    for (int i_distance = 1; i_distance < distance_bins; i_distance++) {
        double dv = pow(i_distance * distance_width, 3) - pow((i_distance - 1) * distance_width, 3);
        for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
            max_value = std::max(max_value, hist.at(std::make_pair(i_distance, i_angle)) / (dv));
        }
    }
    const boost::format fmt("%15.3f %15.3f %15.6f\n");
    for (int i_distance = 1; i_distance < distance_bins; i_distance++) {
        double dv = pow(i_distance * distance_width, 3) - pow((i_distance - 1) * distance_width, 3);
        for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
            double pop = double(hist.at(std::make_pair(i_distance, i_angle))) / (max_value * dv);
            os << boost::format(fmt) % ((i_distance - 0.5) * distance_width) % ((i_angle - 0.5) * angle_width) %
                      (pop == 0.0 ? 100.0 : factor * log(pop));
        }
    }
}

void DistanceAngle::readInfo() {
    select1group(mask1, "Please Enter mask for Atom1(An) > ");
    select1group(mask2, "Please Enter mask for Atom2(OAn) > ");
    select1group(mask3, "Please Enter mask for Atom3(Ow) > ");

    double rmax = choose(0.0, std::numeric_limits<double>::max(),
                         "Enter Maximum Distance to Accumulate[10.0 Ang]:", Default(10.0));
    distance_width =
        choose(0.0, std::numeric_limits<double>::max(), "Enter Width of Distance Bins [0.01 Ang]:", Default(0.01));
    double angle_max = choose(0.0, 180.0, "Enter Maximum Angle to Accumulate[180.0 degree]:", Default(180.0));
    angle_width = choose(0.0, 180.0, "Enter Width of Angle Bins [0.5 degree]:", Default(0.5));
    temperature = choose(0.0, 10000.0, "Temperature [298] (K):", Default(298.0));

    distance_bins = int(rmax / distance_width);
    angle_bins = int(angle_max / angle_width);

    for (int i_distance = 1; i_distance <= distance_bins; i_distance++) {
        for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
            hist[std::make_pair(i_distance, i_angle)] = 0;
        }
    }
}
