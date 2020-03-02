
#include "ana_module/Angle.hpp"

#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm.hpp>

#include "ana_module/RMSDCal.hpp"
#include "data_structure/frame.hpp"
#include "data_structure/molecule.hpp"
#include "utils/PBCUtils.hpp"
#include "utils/common.hpp"

Angle::Angle() {
    enable_outfile = true;
    enable_forcefield = true;
}

void Angle::processFirstFrame(std::shared_ptr<Frame> &frame) {
    boost::for_each(frame->atom_list, [this](std::shared_ptr<Atom> &atom) {
        if (Atom::is_match(atom, mask1)) atom_group1.push_back(atom);
        if (Atom::is_match(atom, mask2)) atom_group2.push_back(atom);
    });
    mol1 = PBCUtils::calculate_intermol(atom_group1, frame);
    mol2 = PBCUtils::calculate_intermol(atom_group2, frame);
}

Eigen::Matrix3d Angle::calculate_inertia(std::shared_ptr<Frame> &frame,
                                         const std::vector<std::shared_ptr<Atom>> &atom_group,
                                         PBCUtils::MolPair &mols) {
    PBCUtils::move(mols, frame);

    std::tuple<double, double, double> mass_center{};
    auto total_mass = 0.0;

    for (auto &atom : atom_group) {
        mass_center += atom->mass.value() * atom->getCoordinate();
        total_mass += atom->mass.get();
    }
    mass_center /= total_mass;

    Eigen::Matrix3d inertia = Eigen::Matrix3d::Zero();
    for (const auto &atom : atom_group) {
        const auto mass = atom->mass.get();
        const auto [x, y, z] = atom->getCoordinate() - mass_center;
        inertia(0, 0) += mass * (y * y + z * z);
        inertia(0, 1) -= mass * x * y;
        inertia(0, 2) -= mass * x * z;
        inertia(1, 1) += mass * (x * x + z * z);
        inertia(1, 2) -= mass * y * z;
        inertia(2, 2) += mass * (x * x + y * y);
    }
    inertia(1, 0) = inertia(0, 1);
    inertia(2, 0) = inertia(0, 2);
    inertia(2, 1) = inertia(1, 2);
    return inertia;
}

std::tuple<double, double, double> Angle::calculate_axis(AxisType type, const Eigen::Matrix3d &inertia) const {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(inertia);
    Eigen::Vector3d eigenvalues = es.eigenvalues();
    Eigen::Matrix3d::Index row;
    type == AxisType::MAX ? eigenvalues.maxCoeff(&row) : eigenvalues.minCoeff(&row);

    Eigen::Vector3d ev = es.eigenvectors().col(row);
    return {ev[0], ev[1], ev[2]};
}

void Angle::process(std::shared_ptr<Frame> &frame) {
    auto vector1 = calculate_axis(type1, calculate_inertia(frame, atom_group1, mol1));
    auto vector2 = calculate_axis(type2, calculate_inertia(frame, atom_group2, mol2));

    auto angle = radian * std::acos(std::abs(dot_multiplication(vector1, vector2)));

    deque.push_back(angle);
}

void Angle::readInfo() {
    auto get_type = [] {
        std::cout << "(0) MIX\n(1) MAX\n";
        return static_cast<AxisType>(choose(0, 1, "choose : "));
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
