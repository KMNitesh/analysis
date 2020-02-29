
#include "ana_module/Angle.hpp"
#include "data_structure/frame.hpp"
#include "utils/common.hpp"
#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm.hpp>

Angle::Angle() {
    enable_outfile = true;
    enable_forcefield = true;
}

void Angle::processFirstFrame(std::shared_ptr<Frame> &frame) {

    boost::for_each(frame->atom_list, [this](std::shared_ptr<Atom> &atom) {
        if (Atom::is_match(atom, mask1))
            atom_group1.push_back(atom);
        if (Atom::is_match(atom, mask2))
            atom_group2.push_back(atom);
    });
}

Eigen::Matrix3d Angle::calculate_inertia(
        std::shared_ptr<Frame> &frame, const std::vector<std::shared_ptr<Atom>> &atom_group) {

    std::vector<std::tuple<double, double, double>> atom_group_positions;
    atom_group_positions.reserve(atom_group.size());

    std::tuple<double, double, double> mass_center{};
    auto pre = atom_group.front()->getCoordinate();
    auto total_mass = 0.0;

    for (auto &atom : atom_group) {
        auto r = atom->getCoordinate() - pre;
        frame->image(r);
        pre += r;
        mass_center += atom->mass.value() * pre;
        total_mass += atom->mass.get();
        atom_group_positions.push_back(pre);
    }
    mass_center /= total_mass;

    Eigen::Matrix3d inertia = Eigen::Matrix3d::Zero();
    for (const auto &element : atom_group | boost::adaptors::indexed()) {
        const auto &position = atom_group_positions[element.index()];
        auto mass = element.value()->mass.get();
        auto[x, y, z] = position - mass_center;
        inertia(0, 0) += mass * ( y * y + z * z);
        inertia(0, 1) -= mass * x * y;
        inertia(0, 2) -= mass * x * z;
        inertia(1, 1) += mass * ( x * x + z * z);
        inertia(1, 2) -= mass * y * z;
        inertia(2, 2) += mass * ( x * x + y * y);
    }
    inertia(1,0) = inertia(0,1);
    inertia(2,0) = inertia(0,2);
    inertia(2,1) = inertia(1,2);
    return inertia;
}

std::tuple<double, double, double> Angle::calculate_axis(AxisType type, const Eigen::Matrix3d &inertia) const {
    Eigen::EigenSolver<Eigen::Matrix3d> es(inertia);
    const auto &eigenvalues = es.eigenvalues();
    int index = 0;
    double value = eigenvalues[0].real();
    for (int i = 0; i < 3; ++i) {
        auto e = eigenvalues[i].real();
        assert(e >= 0.0);
        if (type == AxisType::MAX ? e > value : e < value) {
            index = i;
            value = e;
        }
    }

    auto ev = es.eigenvectors().col(index);
    return {ev[0].real(), ev[1].real(), ev[2].real()};

}

void Angle::process(std::shared_ptr<Frame> &frame) {
    auto vector1 = calculate_axis(type1, calculate_inertia(frame, atom_group1));
    auto vector2 = calculate_axis(type2, calculate_inertia(frame, atom_group2));

    auto angle = radian * std::acos(std::abs(dot_multiplication(vector1 ,vector2)) / std::sqrt(vector_norm2(vector1) * vector_norm2(vector2)));

    deque.push_back(angle);
}

void Angle::readInfo() {

    auto get_type = []{
        std::cout << "(0) MIX\n(1) MAX\n";
        return static_cast<AxisType>(choose(0,1,"choose : "));
    };
    Atom::select1group(mask1, "Enter mask for group1 > ");
    type1 = get_type();

    Atom::select1group(mask2, "Enter mask for group2 > ");
    type2 = get_type();

}

void Angle::print(std::ostream &os) {

    os << std::string(50, '#') << '\n';
    os << "# " << title() << '\n';
    os << "# mask1 > " << mask1 << '\n';
    os << "# type1 > " << (type1 == AxisType::MIN ? "MIN" : "MAX") << '\n';
    os << "# mask2 > " << mask2 << '\n';
    os << "# type2 > " << (type2 == AxisType::MIN ? "MIN" : "MAX") << '\n';
    os << std::string(50, '#') << '\n';
    os << std::setw(10) << "#Frame" << std::setw(15) << "Angle(degree)" << '\n';
    os << std::setprecision(3);
    for (const auto &ele : deque | boost::adaptors::indexed(0)) {
        os << std::setw(10) << ele.index() << std::setw(15) << ele.value() << '\n';
    }
}
