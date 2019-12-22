//
// Created by xiamr on 6/30/19.
//

#include "PlaneAngle.hpp"
#include "data_structure/frame.hpp"
#include "utils/ThrowAssert.hpp"
#include "utils/common.hpp"
#include "json.hpp"

using namespace std;

PlaneAngle::PlaneAngle() {
    enable_outfile = true;
}

void PlaneAngle::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->mask1)) this->group1.insert(atom);
                      if (Atom::is_match(atom, this->mask2)) this->group2.insert(atom);
                      if (Atom::is_match(atom, this->mask3)) this->group3.insert(atom);
                      if (Atom::is_match(atom, this->mask4)) this->group4.insert(atom);
                  });

    for (auto &vec1_atom1 : group1) {
        for (auto &vec1_atom2: group2) {
            if (vec1_atom1->molecule.lock() == vec1_atom2->molecule.lock()) {
                for (auto &vec1_atom3 : group3) {
                    if (vec1_atom1->molecule.lock() == vec1_atom3->molecule.lock()) {
                        pairs.emplace_back(vec1_atom1, vec1_atom2, vec1_atom3);
                    }
                }
            }
        }
    }
}

void PlaneAngle::process(std::shared_ptr<Frame> &frame) {

    for (auto&[vec1_atom1, vec1_atom2, vec1_atom3] : pairs) {

        double u1 = vec1_atom2->x - vec1_atom1->x;
        double u2 = vec1_atom2->y - vec1_atom1->y;
        double u3 = vec1_atom2->z - vec1_atom1->z;

        double v1 = vec1_atom3->x - vec1_atom1->x;
        double v2 = vec1_atom3->y - vec1_atom1->y;
        double v3 = vec1_atom3->z - vec1_atom1->z;

        frame->image(u1, u2, u3);
        frame->image(v1, v2, v3);


        auto xr = u2 * v3 - u3 * v2;
        auto yr = u1 * v3 - u3 * v1;
        auto zr = u1 * v2 - u2 * v1;

        double leng1 = sqrt(xr * xr + yr * yr + zr * zr);

        for (auto &atom4 : group4) {

            double xr2 = atom4->x - vec1_atom1->x;
            double yr2 = atom4->y - vec1_atom1->y;
            double zr2 = atom4->z - vec1_atom1->z;

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

void PlaneAngle::print(std::ostream &os) {
    os << string(50, '#') << '\n';
    os << "# " << PlaneAngle::title() << '\n';
    os << "# Group1 > " << mask1 << '\n';
    os << "# Group2 > " << mask2 << '\n';
    os << "# Group3 > " << mask3 << '\n';
    os << "# Group4 > " << mask4 << '\n';
    os << "# angle_width(degree) > " << angle_width << '\n';
    os << "# Cutoff1(Ang) > " << cutoff1 << '\n';
    os << "# Cutoff2(Ang) > " << cutoff2 << '\n';
    os << string(50, '#') << '\n';
    os << format("#%15s %15s\n", "Angle(degree)", "Probability Density(% degree-1)");

    printData(os);

    os << string(50, '#') << '\n';

    os << ">>>JSON<<<\n";
    saveJson(os);
    os << "<<<JSON>>>\n";
}

void PlaneAngle::readInfo() {
    Atom::select1group(mask1, "Enter mask for atom1 (Ow) : ");
    Atom::select1group(mask2, "Enter mask for atom2 (Hw) : ");
    Atom::select1group(mask3, "Enter mask for atom3 (Hw) : ");
    Atom::select1group(mask4, "Enter mask for atom4 (Metal) : ");

    double angle_max = choose(0.0, 180.0, "Enter Maximum Angle to Accumulate[180.0 degree]:", Default(180.0));
    angle_width = choose(0.0, 180.0, "Enter Width of Angle Bins [0.5 degree]:", Default(0.5));

    cutoff1 = choose(0.0, 100.0, "Cutoff1 [Angstrom]:");
    cutoff2 = choose(0.0, 100.0, "Cutoff2 [Angstrom]:");

    throw_assert(cutoff1 < cutoff2, "Cutoff1 must less than Cutoff2");

    angle_bins = int(angle_max / angle_width);

    for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
        hist[i_angle] = 0;
    }
}

void PlaneAngle::printData(std::ostream &os) const {
    double total = 0.0;

    for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
        total += hist.at(i_angle);
    }

    for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
        os << format("%15.3f %15.3f\n", (i_angle - 0.5) * angle_width, 100 * (hist.at(i_angle) / total) / angle_width);
    }
}

void PlaneAngle::saveJson(std::ostream &os) const {

    nlohmann::json json;

    json["title"] = title();
    json["group1"] = to_string(mask1);
    json["group2"] = to_string(mask2);
    json["group3"] = to_string(mask3);
    json["group4"] = to_string(mask4);
    json["angle_width"] = {{"value", angle_width},
                           {"unit",  "degree"}};

    json["cutoff1"] = {{"value", cutoff1},
                       {"unit",  "Ang"}};

    json["cutoff2"] = {{"value", cutoff2},
                       {"unit",  "Ang"}};

    json["Probability Density"] = {
            {"meta", {
                             {"X", "Angle(degree)"},
                             {"Y", "Probability Density(degree-1)"}
                     }}
    };

    double total = 0.0;
    for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
        total += hist.at(i_angle);
    }

    for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
        json["Probability Density"]["values"].push_back(
                {(i_angle - 0.5) * angle_width, 100 * (hist.at(i_angle) / total) / angle_width});
    }

    os << json;
}
