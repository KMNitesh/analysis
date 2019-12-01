//
// Created by xiamr on 8/5/19.
//

#include "SpatialOrientationDistribution.hpp"
#include "data_structure/frame.hpp"
#include "utils/ThrowAssert.hpp"
#include "data_structure/molecule.hpp"
#include "utils/common.hpp"

using namespace std;

SpatialOrientationDistribution::SpatialOrientationDistribution() {
    enable_forcefield = true;
    enable_outfile = true;
}

void SpatialOrientationDistribution::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->ids)) this->group.insert(atom->molecule.lock());
                  });
}

void SpatialOrientationDistribution::process(std::shared_ptr<Frame> &frame) {
    for (auto &mol : group) {
        auto dipole = mol->calc_dipole(frame);
        dipole /= vector_norm(dipole);
        normal_vectors.push_back(dipole);
    }
}

void SpatialOrientationDistribution::print(std::ostream &os) {
    os << string(50, '#') << '\n';
    os << "# " << title() << '\n';
    os << "# Group > " << ids << '\n';
    os << string(50, '#') << '\n';
    os << format("#%15s %15s %15s\n", "X", "Y", "Z");

    for (auto[x, y, z] : normal_vectors) {
        os << format("%15.3f %15.3f %15.3f\n", x, y, z);
    }

    os << string(50, '#') << '\n';
}

void SpatialOrientationDistribution::readInfo() {
    Atom::select1group(ids);
}

double SpatialOrientationDistribution::calculatePhiAngle(const std::tuple<double, double, double> &vector) {

    constexpr std::tuple<double, double, double> z_axis = {0.0, 0.0, 1.0};

    auto norm = vector_norm(vector);

    throw_assert(norm > 0, "Zero Vector");

    return radian * acos(dot_multiplication(vector, z_axis) / norm);
}

double SpatialOrientationDistribution::calculateThetaAngle(const std::tuple<double, double, double> &vector) {
    constexpr std::tuple<double, double, double> x_axis = {1.0, 0.0, 0.0};
    constexpr std::tuple<double, double, double> y_axis = {0.0, 1.0, 0.0};

    std::tuple<double, double, double> vector_project_to_xyplane = {std::get<0>(vector), std::get<1>(vector), 0.0};
    auto norm = vector_norm(vector_project_to_xyplane);

    throw_assert(norm > 0, "Zero Vector");
    vector_project_to_xyplane /= norm;

    auto angle_with_xaxis = radian * acos(dot_multiplication(vector_project_to_xyplane, x_axis));

    auto angle_with_yaxis = radian * acos(dot_multiplication(vector_project_to_xyplane, y_axis));

    if (angle_with_yaxis > 90) {
        angle_with_xaxis = 360 - angle_with_xaxis;
    }
    return angle_with_xaxis;
}




