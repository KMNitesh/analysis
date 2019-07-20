//
// Created by xiamr on 7/8/19.
//

#include "DipoleAngleAxis3D.hpp"
#include "frame.hpp"

using namespace std;

void DipoleAngleAxis3D::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->ids)) this->group.insert(atom->molecule.lock());
                  });
}

void DipoleAngleAxis3D::process(std::shared_ptr<Frame> &frame) {
    constexpr tuple<double, double, double> xaxis = {1, 0, 0}, yaxis = {0, 1, 0}, zaxis = {0, 0, 1};
    for (auto &mol: group) {
        auto dipole = mol->calc_dipole(frame);
        dipole /= vector_norm(dipole);

        auto angle_x = acos(dot_multiplication(dipole, xaxis)) * radian;
        auto angle_y = acos(dot_multiplication(dipole, yaxis)) * radian;
        auto angle_z = acos(dot_multiplication(dipole, zaxis)) * radian;

        distributions.emplace_back(angle_x, angle_y, angle_z);
    }
}

void DipoleAngleAxis3D::print(std::ostream &os) {
    os << string(50, '#') << '\n';
    os << "# " << title() << '\n';
    os << "# Group > " << ids << '\n';
    os << string(50, '#') << '\n';
    os << format("#%15s %15s %15s\n", "Angle_X(degree)", "Angle_Y(degree)", "Angle_Z(degree)");

    printData(os);

    os << string(50, '#') << '\n';
}

void DipoleAngleAxis3D::readInfo() {
    Atom::select1group(ids);
}

void DipoleAngleAxis3D::printData(std::ostream &os) const {
    for (auto[angle_x, angle_y, angle_z] : distributions) {
        os << format("%15.3f %15.3f %15.3f\n", angle_x, angle_y, angle_z);
    }
}
