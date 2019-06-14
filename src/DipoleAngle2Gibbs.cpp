//
// Created by xiamr on 6/14/19.
//

#include "DipoleAngle2Gibbs.hpp"

using namespace std;

void DipoleAngle2Gibbs::print() {

    double factor = -kb * temperature * avogadro_constant / 4184.0;
    double max_value = 0.0;

    for (int i_distance = 1; i_distance < distance_bins; i_distance++) {
        double dv = pow(i_distance * distance_width, 3) - pow((i_distance - 1) * distance_width, 3);
        for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
            max_value = max(max_value, hist[make_pair(i_distance, i_angle)] / (dv));
        }
    }
    for (int i_distance = 1; i_distance < distance_bins; i_distance++) {
        double dv = pow(i_distance * distance_width, 3) - pow((i_distance - 1) * distance_width, 3);
        for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
            double pop = double(hist[make_pair(i_distance, i_angle)]) / (max_value * dv);
            outfile << boost::format("%5.3f   %10.3f   %15.6f\n")
                       % ((i_distance - 0.5) * distance_width)
                       % ((i_angle - 0.5) * angle_width)
                       % (pop == 0.0 ? 100.0 : factor * log(pop));
            //               %  pop;
        }
    }

    /*for (int i_distance = 1; i_distance < distance_bins; i_distance++) {
        for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
            max_value = max(max_value,hist[make_pair(i_distance, i_angle)] / (angle_width));
        }
    }
    for (int i_distance = 1; i_distance < distance_bins; i_distance++) {
        double dv = pow(i_distance * distance_width, 3) - pow((i_distance - 1) * distance_width, 3);
        for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
            double pop = double(hist[make_pair(i_distance, i_angle)]) / (max_value * angle_width);
            outfile << boost::format("%5.3f%7.1f%10.6f\n")
                       % (i_distance * distance_width)
                       % (i_angle * angle_width)
                     //  % (pop == 0.0 ? 10000.0 : factor * log(pop));
                       %  pop;
        }
    }
    */
}

void DipoleAngle2Gibbs::readInfo() {
    DipoleAngle::readInfo();
    temperature = choose(0.0, 10000.0, "Temperature [298] (K):", true, 298.0);
}

