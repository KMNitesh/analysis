//
// Created by xiamr on 8/5/19.
//

#include <gmock/gmock.h>
#include "SpatialOrientationDistribution.hpp"
#include "ThrowAssert.hpp"

using namespace std;
using namespace testing;

TEST(SpatialOrientationDistributionTest, calculatePhiAngleWithNomalOrientation) {
    ASSERT_THAT(SpatialOrientationDistribution::calculatePhiAngle({1, 1, 1}),
                DoubleEq(radian * acos(1 / sqrt(3))));
}

TEST(SpatialOrientationDistributionTest, calculatePhiAngleWithZeroDegree) {
    ASSERT_THAT(SpatialOrientationDistribution::calculatePhiAngle({0, 0, 1}),
                DoubleEq(0));
}

TEST(SpatialOrientationDistributionTest, calculatePhiAngleWith180Degree) {
    ASSERT_THAT(SpatialOrientationDistribution::calculatePhiAngle({0, 0, -1}),
                DoubleEq(180));
}

TEST(SpatialOrientationDistributionTest, calculatePhiAngleWithZeroVector) {
    ASSERT_THROW(SpatialOrientationDistribution::calculatePhiAngle({0, 0, 0}), AssertionFailureException);
}

TEST(SpatialOrientationDistributionTest, calculateThetaAngleWithZeroDegree) {
    ASSERT_THAT(SpatialOrientationDistribution::calculateThetaAngle({1, 0, 0}),
                DoubleEq(0));
}

TEST(SpatialOrientationDistributionTest, calculateThetaAngleWith45DegreeinPlane) {
    ASSERT_THAT(SpatialOrientationDistribution::calculateThetaAngle({1, 1, 0}),
                DoubleEq(45));
}

TEST(SpatialOrientationDistributionTest, calculateThetaAngleWith90DegreeinPlane) {
    ASSERT_THAT(SpatialOrientationDistribution::calculateThetaAngle({0, 1, 0}),
                DoubleEq(90));
}

TEST(SpatialOrientationDistributionTest, calculateThetaAngleWith135DegreeinPlane) {
    ASSERT_THAT(SpatialOrientationDistribution::calculateThetaAngle({-1, 1, 0}),
                DoubleEq(135));
}

TEST(SpatialOrientationDistributionTest, calculateThetaAngleWith180DegreeinPlane) {
    ASSERT_THAT(SpatialOrientationDistribution::calculateThetaAngle({-1, 0, 0}),
                DoubleEq(180));
}

TEST(SpatialOrientationDistributionTest, calculateThetaAngleWith225DegreeinPlane) {
    ASSERT_THAT(SpatialOrientationDistribution::calculateThetaAngle({-1, -1, 0}),
                DoubleEq(225));
}

TEST(SpatialOrientationDistributionTest, calculateThetaAngleWith270DegreeinPlane) {
    ASSERT_THAT(SpatialOrientationDistribution::calculateThetaAngle({0, -1, 0}),
                DoubleEq(270));
}

TEST(SpatialOrientationDistributionTest, calculateThetaAngleWith315DegreeinPlane) {
    ASSERT_THAT(SpatialOrientationDistribution::calculateThetaAngle({1, -1, 0}),
                DoubleEq(315));
}

TEST(SpatialOrientationDistributionTest, calculateThetaAngleWith315DegreeSpatial) {
    ASSERT_THAT(SpatialOrientationDistribution::calculateThetaAngle({1, -1, 1}),
                DoubleEq(315));
}

TEST(SpatialOrientationDistributionTest, calculateThetaAngleWithZeroVector) {
    ASSERT_THROW(SpatialOrientationDistribution::calculateThetaAngle({0, 0, 1}), AssertionFailureException);
}










