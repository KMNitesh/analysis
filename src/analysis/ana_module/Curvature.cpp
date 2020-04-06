

#include "Curvature.hpp"
#include "data_structure/frame.hpp"
#include "utils/common.hpp"
#include <Eigen/Eigen>
#include <boost/range/adaptors.hpp>

Curvature::Curvature() { enable_outfile = true; }

void Curvature::processFirstFrame(std::shared_ptr<Frame> &frame) {
    atoms = PBCUtils::find_atoms(mask, frame);

    mol = PBCUtils::calculate_intermol(atoms, frame);
}

double Curvature::calculate_curvature() {
    Eigen::Vector4d x0 = {0.0, 0.0, 0.0, 40.0};

    // residuals
    Eigen::VectorXd r = Eigen::VectorXd::Zero(atoms.size());

    // Jacobi Matrix
    Eigen::MatrixXd Dr = Eigen::MatrixXd::Zero(atoms.size(), 4);

    int iteration = 0;
    do {
        for (auto i : boost::irange(atoms.size())) {
            const auto &atom = atoms[i];

            auto dx = x0[0] - atom->x;
            auto dy = x0[1] - atom->y;
            auto dz = x0[2] - atom->z;

            auto dr2 = dx * dx + dy * dy + dz * dz;

            // auto dr = std::sqrt(dr2);
            r[i] = dr2 - x0[3];

            Dr(i, 0) = 2 * dx;
            Dr(i, 1) = 2 * dy;
            Dr(i, 2) = 2 * dz;
            Dr(i, 3) = -1;
        }
        Eigen::MatrixXd DrT = Dr.transpose();

        Eigen::MatrixXd LeftA = DrT * Dr;
        Eigen::MatrixXd b = -DrT * r;

        Eigen::Vector4d v = LeftA.colPivHouseholderQr().solve(b);

        x0 += v;

        ++iteration;

        // std::cout << boost::format("it = %|5$-5| (x,y,z,R) = (%1$.10f,%2$.10f,%3$.10f,%4$.10f) %|50t| |v| = %6%\n") %
        //                  x0[0] % x0[1] % x0[2] % x0[3] % iteration % v.norm();
        if (std::abs(v.norm()) < 1e-5)
            break;
    } while (iteration < 100);
    if (iteration == 100) {
        std::cerr << "optimization error\n";
        std::exit(0);
    }
    return std::sqrt(x0[3]);
}

void Curvature::process(std::shared_ptr<Frame> &frame) {
    PBCUtils::move(mol, frame);
    curvature_lists.push_back(calculate_curvature());
}

void Curvature::print(std::ostream &os) {

    os << std::string(20, '#') << '\n';
    os << "# mask" << mask << '\n';
    os << '#' << title() << '\n';
    os << std::string(20, '#') << '\n';
    os << boost::format("#%14s %14s %14s\n") % "Frame" % "R" % "Curvature(Ang-1)";
    for (const auto &element : curvature_lists | boost::adaptors::indexed(1)) {
        os << boost::format("%15d %14g %14g\n") % element.index() % element.value() % (1 / element.value());
    }
    os << std::string(20, '#') << '\n';
}

void Curvature::readInfo() { select1group(mask, "select atoms for curvature > "); }