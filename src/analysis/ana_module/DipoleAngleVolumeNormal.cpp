//
// Created by xiamr on 6/14/19.
//

#include "DipoleAngleVolumeNormal.hpp"
#include "data_structure/frame.hpp"
#include "utils/common.hpp"

using namespace std;

void DipoleAngleVolumeNormal::print(std::ostream &os) {
    os << "DipoleAngleVolumeNormal : " << std::endl;
    size_t total = 0;
    for (int i_distance = 1; i_distance < distance_bins; i_distance++) {
        for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
            total += hist[make_pair(i_distance, i_angle)];
        }

    }
    double factor = (4.0 / 3.0) * M_PI;
    for (int i_distance = 1; i_distance < distance_bins; i_distance++) {
        double dv = factor * (pow(i_distance * distance_width, 3) - pow((i_distance - 1) * distance_width, 3));
        for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
            os << boost::format("%5.3f   %10.3f   %g\n")
                  % ((i_distance - 0.5) * distance_width)
                  % ((i_angle - 0.5) * angle_width)
                  % (total == 0 ? 0.0 : double(hist[make_pair(i_distance, i_angle)]) / (total * dv * angle_width));
        }
    }
}

