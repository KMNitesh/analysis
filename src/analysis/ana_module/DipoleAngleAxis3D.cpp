
#include "DipoleAngleAxis3D.hpp"

#include <boost/range/algorithm.hpp>

#include "data_structure/frame.hpp"
#include "data_structure/molecule.hpp"
#include "utils/common.hpp"

DipoleAngleAxis3D::DipoleAngleAxis3D() {
    enable_forcefield = true;
    enable_outfile = true;
}

void DipoleAngleAxis3D::processFirstFrame(std::shared_ptr<Frame> &frame) {
    boost::for_each(frame->atom_list, [this](std::shared_ptr<Atom> &atom) {
        if (is_match(atom, this->amberMask)) this->group.insert(atom->molecule.lock());
    });
}

void DipoleAngleAxis3D::process(std::shared_ptr<Frame> &frame) {
    constexpr std::tuple<double, double, double> xaxis = {1, 0, 0}, yaxis = {0, 1, 0}, zaxis = {0, 0, 1};
    for (auto &mol : group) {
        auto dipole = mol->calc_dipole(frame);
        dipole /= vector_norm(dipole);

        auto angle_x = acos(dot_multiplication(dipole, xaxis)) * radian;
        auto angle_y = acos(dot_multiplication(dipole, yaxis)) * radian;
        auto angle_z = acos(dot_multiplication(dipole, zaxis)) * radian;

        distributions.emplace_back(angle_x, angle_y, angle_z);
    }
}

void DipoleAngleAxis3D::print(std::ostream &os) {
    os << std::string(50, '#') << '\n';
    os << "# " << title() << '\n';
    os << "# Group > " << amberMask << '\n';
    os << std::string(50, '#') << '\n';
    os << format("#%15s %15s %15s\n", "Angle_X(degree)", "Angle_Y(degree)", "Angle_Z(degree)");

    printData(os);

    os << std::string(50, '#') << '\n';
}

void DipoleAngleAxis3D::readInfo() { select1group(amberMask); }

void DipoleAngleAxis3D::printData(std::ostream &os) const {
    const boost::format fmt("%15.3f %15.3f %15.3f\n");
    for (auto [angle_x, angle_y, angle_z] : distributions) {
        os << boost::format(fmt) % angle_x % angle_y % angle_z;
    }
}
