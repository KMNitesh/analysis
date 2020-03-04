
#include "others/HOOH_Calculator.hpp"

#include <Eigen/Eigen>

#include "utils/common.hpp"

void HOOH_Calculator::process() {
    Eigen::Vector3d A, B, C;

    std::string line;
    std::stringstream ss;

    std::cout << "input OH2: ";
    std::getline(std::cin, line);
    ss.str(line);
    ss >> A[0] >> A[1] >> A[2];

    std::cout << "input OH1: ";
    std::getline(std::cin, line);
    ss.str(line);
    ss >> B[0] >> B[1] >> B[2];

    std::cout << "input HH1: ";
    std::getline(std::cin, line);
    ss.str(line);
    ss >> C[0] >> C[1] >> C[2];

    double l_BA = choose(0.0, "input length(O-O)[0.1474 nm]: ", Default(0.1474));
    double l_BC = choose(0.0, "input length(O-H)[0.0950 nm]: ", Default(0.0950));

    double theta = choose(0.0, "input angle(O-O-H)[94.8]: ", Default(94.8));
    double phi = choose(0.0, "input diheral(H-O-O-H)[115.5]: ", Default(115.5));

    auto [D1, D2] = calculate(A, B, C, l_BA, l_BC, theta, phi);

    std::cout << std::fixed << std::setprecision(3);
    using std::setw;
    std::cout << " Results\n";
    std::cout << "OH2 > " << setw(8) << A[0] << setw(8) << A[1] << setw(8) << A[2] << '\n';
    std::cout << "HH1 > " << setw(8) << C[0] << setw(8) << C[1] << setw(8) << C[2] << '\n';

    std::cout << "HH2  > " << setw(8) << D1[0] << setw(8) << D1[1] << setw(8) << D1[2] << '\n';
    std::cout << "HH2' > " << setw(8) << D2[0] << setw(8) << D2[1] << setw(8) << D2[2] << '\n';
}

std::array<Eigen::Vector3d, 2> HOOH_Calculator::calculate(Eigen::Vector3d &A, const Eigen::Vector3d &B,
                                                          Eigen::Vector3d &C, double l_BA, double l_BC, double theta,
                                                          double phi) {
    using namespace Eigen;

    Vector3d BA = A - B;
    Vector3d x = BA / BA.norm();

    Vector3d c = x.cross(C - B);
    Vector3d y = c / c.norm();

    Vector3d z = x.cross(y);

    A = l_BA * x + B;

    C = l_BC * (std::cos(theta / radian) * x - std::sin(theta / radian) * z) + B;

    Vector3d D1 =
        l_BC * (std::cos((180 - theta) / radian) * x - std::sin((180 - theta) / radian) * std::sin(phi / radian) * y -
                std::sin((180 - theta) / radian) * std::cos(phi / radian) * z) +
        A;

    Vector3d D2 =
        l_BC * (std::cos((180 - theta) / radian) * x + std::sin((180 - theta) / radian) * std::sin(phi / radian) * y -
                std::sin((180 - theta) / radian) * std::cos(phi / radian) * z) +
        A;
    return {D1, D2};
}
