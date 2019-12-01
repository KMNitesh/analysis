//
// Created by xiamr on 6/14/19.
//

#include "DipoleAngle2Gibbs.hpp"
#include "utils/common.hpp"

using namespace std;

void DipoleAngle2Gibbs::print(std::ostream &os) {


    os << string(50, '#') << '\n';
    os << "# " << DipoleAngle2Gibbs::title() << '\n';
    os << "# Group1 > " << ids1 << '\n';
    os << "# Group2 > " << ids2 << '\n';
    os << "# distance_width(Ang) > " << distance_width << '\n';
    os << "# angle_width > " << angle_width << '\n';
    os << "# Temperature(K) > " << temperature << '\n';
    os << string(50, '#') << '\n';
    os << format("#%15s %15s %15s\n", "Distance(Ang)", "Angle(degree)", "Energy(kcal/mol)");

    printData(os);

    os << string(50, '#') << '\n';
}

void DipoleAngle2Gibbs::printData(ostream &os) const {
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

void DipoleAngle2Gibbs::readInfo() {
    DipoleAngle::readInfo();
    temperature = choose(0.0, 10000.0, "Temperature [298] (K):", Default(298.0));
}

