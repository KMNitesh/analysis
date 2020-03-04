
#include <gmock/gmock.h>

#include "others/HOOH_Calculator.hpp"
#include "utils/common.hpp"

using namespace testing;
using namespace Eigen;

TEST(Test_HOOH_Calculator, defaultValue) {
    Vector3d A = {4.119, 5.621, 5.206};  // OH2
    Vector3d B = {4.173, 5.756, 5.196};  // OH1
    Vector3d C = {4.146, 5.808, 5.277};  // HH1

    double l_BA = 0.146;
    double l_BC = 0.1;

    double theta = 98;
    double phi = 120;

    auto [D1, D2] = HOOH_Calculator::calculate(A, B, C, l_BA, l_BC, theta, phi);

    ASSERT_THAT((A - B).norm(), DoubleNear(0.146, 0.001));
    ASSERT_THAT((C - B).norm(), DoubleNear(0.100, 0.001));
    ASSERT_THAT((A - D1).norm(), DoubleNear(0.100, 0.001));
    ASSERT_THAT((A - D2).norm(), DoubleNear(0.100, 0.001));

    ASSERT_THAT(std::acos(((C - B).dot(A - B)) / ((C - B).norm() * (A - B).norm())) * radian, DoubleNear(theta, 0.1));

    ASSERT_THAT(std::acos(((D1 - A).dot(B - A)) / ((D1 - A).norm() * (B - A).norm())) * radian, DoubleNear(theta, 0.1));

    ASSERT_THAT(std::acos(((D2 - A).dot(B - A)) / ((D2 - A).norm() * (B - A).norm())) * radian, DoubleNear(theta, 0.1));

    Vector3d v1 = (A - B).cross(C - B);
    Vector3d v2 = (D1 - A).cross(B - A);
    Vector3d v3 = (D2 - A).cross(B - A);

    ASSERT_THAT(std::acos(v1.dot(v2) / (v1.norm() * v2.norm())) * radian, DoubleNear(phi, 0.1));
    ASSERT_THAT(std::acos(v1.dot(v3) / (v1.norm() * v3.norm())) * radian, DoubleNear(phi, 0.1));
}
