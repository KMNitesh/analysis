//
// Created by xiamr on 10/11/19.
//

#include "std.hpp"
#include <gmock/gmock.h>
#include <Eigen/Eigen>
#include <boost/format.hpp>
#include <boost/range/algorithm.hpp>
#include "common.hpp"
#include "CoordinationStructureMatch.hpp"

using namespace testing;

TEST(TestNonlinearLestSquaresLearn, ThreeCircles) {

    Eigen::Vector2d x0 = {0, 0};

    double x[] = {-1, 1, 1};
    double y[] = {0, 0.5, -0.5};
    double R[] = {1, 0.5, 0.5};

    double S[3];

    //residuals
    Eigen::Vector3d r;

    //Jacobi Matrix
    Eigen::Matrix<double, 3, 2> Dr;

    int iteration = 0;
    do {
        for (int i = 0; i < 3; ++i) {
            auto dx = x0[0] - x[i];
            auto dy = x0[1] - y[i];
            S[i] = std::sqrt(dx * dx + dy * dy);
            r[i] = S[i] - R[i];

            Dr(i, 0) = dx / S[i];
            Dr(i, 1) = dy / S[i];
        }

        auto DrT = Dr.transpose();

        auto LeftA = DrT * Dr;
        auto b = -DrT * r;

        auto v = LeftA.colPivHouseholderQr().solve(b);

        x0 += v;

        ++iteration;

        std::cout << boost::format("it = %|3$-5| (x,y) = (%1$.10f,%2$.10f) %|50t| |v| = %4%\n")
                     % x0[0] % x0[1] % iteration % v.norm();

    } while (iteration < 10);

}


class TestNonlinearLestSquares : public Test {
protected:

    std::vector<std::tuple<double, double, double>> casp_coord{
            // CASP
            {-1.32663504, -0.50082192, -0.84008400},
            {1.74902742,  0.05760216,  -0.28925004},
            {0.86807377,  -0.49465605, 1.35448181},
            {0.42678461,  -0.49484127, -1.49144034},
            {-0.91779692, -0.81890269, 1.03464236},
            {0.39656900,  -1.72124744, -0.07942184},
            {-1.13306883, 0.92261401,  0.47113115},
            {0.65280186,  1.24686064,  0.79097061},
            {0.02582338,  1.27802768,  -1.04991219}
    };

    std::vector<std::tuple<double, double, double>> tctp_coord{
            // TCTP
            {-1.29095254, 0.81022141,  -0.18248319},
            {0.05763523,  -1.21493328, 1.20098413},
            {0.67796084,  -0.69851645, -1.47902186},
            {-0.07730040, 0.42368972,  0.89991192},
            {0.29448331,  0.73319712,  -0.70631320},
            {1.10274104,  -0.48055178, 0.12284911},
            {-0.95570364, -0.45317847, -1.22429238},
            {-0.14744591, -1.66692738, -0.39513006},
            {-1.32748735, -0.76268587, 0.38193275}

    };
};

TEST_F(TestNonlinearLestSquares, TCTP) {
    ASSERT_THAT(CoordinationStructureMatch::testTCTP(tctp_coord), Lt(CoordinationStructureMatch::testCASP(tctp_coord)));
}


TEST_F(TestNonlinearLestSquares, CASP) {
    ASSERT_THAT(CoordinationStructureMatch::testCASP(casp_coord), Lt(CoordinationStructureMatch::testTCTP(casp_coord)));
}

