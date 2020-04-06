//
// Created by xiamr on 6/14/19.
//

#include "ShellDensity.hpp"

#include "data_structure/frame.hpp"
#include "nlohmann/json.hpp"
#include "utils/common.hpp"

using namespace std;

ShellDensity::ShellDensity() { enable_outfile = true; }

void ShellDensity::process(std::shared_ptr<Frame> &frame) {
    nframe++;
    for (auto &ref : group1) {
        for (auto &atom : group2) {
            int ibin = int(atom_distance(ref, atom, frame) / distance_width) + 1;
            if (ibin <= distance_bins) {
                hist[ibin]++;
            }
        }
    }
}

void ShellDensity::print(std::ostream &os) {
    os << std::string(50, '#') << '\n';
    os << "# " << title() << " # \n";
    os << "# First Mask : " << mask1 << " Second Mask : " << mask2 << endl;
    os << std::string(50, '#') << '\n';
    os << "# Bin    Distance(Ang)    Density (count / Ang3 / frame)" << endl;
    for (int i = 1; i <= distance_bins; i++) {
        double dv = (4.0 / 3.0) * M_PI * (pow(i * distance_width, 3) - pow((i - 1) * distance_width, 3));
        os << boost::format("%d      %.4f      %g \n") % i % ((i - 0.5) * distance_width) % (hist[i] / (nframe * dv));
    }
    os << std::string(50, '#') << '\n';

    os << ">>>JSON<<<\n";
    saveJson(os);
    os << "<<<JSON>>>\n";
}

void ShellDensity::saveJson(std::ostream &os) const {
    nlohmann::json json;

    json["title"] = title();
    json["mask1"] = to_string(mask1);
    json["mask2"] = to_string(mask2);

    json["Density"] = {{"X", {{"name", "Bin"}}},
                       {"Y1", {{"name", "Distance"}, {"unit", "Ang"}}},
                       {"Y2", {{"name", "Density"}, {"unit", "count / Ang3 / frame"}}}};
    for (int i = 1; i <= distance_bins; i++) {
        double dv = (4.0 / 3.0) * M_PI * (pow(i * distance_width, 3) - pow((i - 1) * distance_width, 3));

        json["Density"]["X"]["values"].push_back(i);
        json["Density"]["Y1"]["values"].push_back((i - 0.5) * distance_width);
        json["Density"]["Y2"]["values"].push_back(hist.at(i) / (nframe * dv));
    }

    os << json;
}

void ShellDensity::readInfo() {
    select2group(mask1, mask2);
    double rmax = choose(0.0, std::numeric_limits<double>::max(),
                         "Enter Maximum Distance to Accumulate[10.0 Ang]:", Default(10.0));
    distance_width =
        choose(0.0, std::numeric_limits<double>::max(), "Enter Width of Distance Bins [0.01 Ang]:", Default(0.01));
    distance_bins = int(rmax / distance_width);
    for (int i = 1; i <= distance_bins; i++) {
        hist[i] = 0;
    }
}

void ShellDensity::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(), [this](shared_ptr<Atom> &atom) {
        if (is_match(atom, mask1)) group1.insert(atom);
        if (is_match(atom, mask2)) group2.insert(atom);
    });
}
