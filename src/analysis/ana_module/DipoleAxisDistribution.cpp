//
// Created by xiamr on 6/30/19.
//

#include "DipoleAxisDistribution.hpp"
#include "data_structure/frame.hpp"
#include "data_structure/molecule.hpp"
#include "utils/common.hpp"

DipoleAxisDistribution::DipoleAxisDistribution() {
    enable_forcefield = true;
    enable_outfile = true;
}

void DipoleAxisDistribution::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](std::shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->ids)) this->group.insert(atom->molecule.lock());
                  });
}

void DipoleAxisDistribution::process(std::shared_ptr<Frame> &frame) {
    for (auto &mol: group) {
        auto dipole = mol->calc_dipole(frame);
        dipole /= vector_norm(dipole);
        auto angle = acos(dot_multiplication(dipole, axis_vector)) * radian;
        hist.update(angle);
    }
}

void DipoleAxisDistribution::print(std::ostream &os) {
    os << std::string(50, '#') << '\n';
    os << "# " << title() << '\n';
    os << "# Group > " << ids << '\n';
    os << "# angle_width(degree) > " << hist.getWidth() << '\n';
    const std::unordered_map<int, std::string> mapping{
            {1, "X-Axis"},
            {2, "Y-Axis"},
            {3, "Z-Axis"}
    };
    os << "# " << mapping.at(axis_num) << '\n';
    os << std::string(50, '#') << '\n';
    os << format("#%15s %15s\n", "Angle(degree)", "Probability Density(% degree-1)");

    printData(os);

    os << std::string(50, '#') << '\n';
}

void DipoleAxisDistribution::readInfo() {
    Atom::select1group(ids);
    double angle_max = choose(0.0, 180.0, "Enter Maximum Angle to Accumulate[180.0 degree]:", Default(180.0));
    auto angle_width = choose(0.0, 180.0, "Enter Width of Angle Bins [0.5 degree]:", Default(0.5));

    auto menu = [] {
        std::cout << "Axis selection \n"
                  << "  1. X-Axis\n"
                  << "  2. Y-Axis\n"
                  << "  3. Z-Axis\n";

        return choose(1, 3, " select > ");
    };

    const std::unordered_map<int, std::tuple<double, double, double>> mapping{
            {1, xaxis_vector},
            {2, yaxis_vector},
            {3, zaxis_vector}
    };
    axis_num = menu();
    axis_vector = mapping.at(axis_num);
    hist.initialize(angle_max, angle_width);
}

void DipoleAxisDistribution::printData(std::ostream &os) const {
    const boost::format fmt("%15.3f %15.3f\n");
    for (auto[grid, value] : hist.getDistribution()) {
        os << boost::format(fmt) % grid % (100 * value);
    }
}