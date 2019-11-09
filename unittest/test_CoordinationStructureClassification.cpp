//
// Created by xiamr on 9/22/19.
//

#include "std.hpp"
#include <gmock/gmock.h>

#include <random>
#include "CoordinationStructureClassification.hpp"
#include <boost/range/algorithm.hpp>

using namespace testing;

class TestCoordinationStructureClassification : public testing::Test {
public:
    std::vector<std::tuple<double, double, double>> c1{
            {-1.91, 0.21,  -1.44},
            {2.02,  -1.35, 0.09},
            {0.73,  -0.39, 2.25},
            {0.03,  -1.66, -1.84},
            {-1.98, 0.13,  1.41},
            {-0.62, -2.24, 0.7},
            {-0.59, 2.33,  0.17},
            {1.88,  1.55,  0.54},
            {0.96,  1.08,  -2.06}
    };
};

TEST_F(TestCoordinationStructureClassification, calculateRmsdOfTwoStructsWithSameCoorindation) {
    auto c2 = c1;
    ASSERT_THAT(CoordinationStructureClassification::calculateRmsdOfTwoStructs(c1, c2), DoubleEq(0));
}

TEST(CoordinationStructureClassification, calculateRmsdOfTwoStructsWithRandomCoorindation) {
    std::vector<std::tuple<double, double, double>> c2(8);
    std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> distribution(-4, 4);
    auto dice = std::bind(distribution, generator);

    for (int i = 0; i < 10000; ++i) {
        for (auto &element : c2) {
            element = {dice(), dice(), dice()};
        }
        auto c3 = c2;
        ASSERT_THAT(CoordinationStructureClassification::calculateRmsdOfTwoStructs(c3, c2), DoubleEq(0));
    }
}

TEST_F(TestCoordinationStructureClassification, calculateRmsdOfTwoStructsWithRandomShuffle) {
    auto c2 = c1;
    for (int i = 0; i < 100000; ++i) {
        std::shuffle(c2.begin(), c2.end(), std::mt19937(std::random_device()()));
        ASSERT_THAT(CoordinationStructureClassification::calculateRmsdOfTwoStructs(c1, c2), DoubleEq(0));
    }
}


TEST(CoordinationStructureClassification, calculateRmsdOfTwoStructsWithRandomCoorindationWithRandomSuffle) {
    std::vector<std::tuple<double, double, double>> c2(8);
    std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> distribution(-4, 4);
    auto dice = std::bind(distribution, generator);

    for (int i = 0; i < 100; ++i) {
        for (auto &element : c2) {
            element = {dice(), dice(), dice()};
        }
        auto c3 = c2;
        for (int j = 0; j < 100; ++j) {
            std::shuffle(c3.begin(), c3.end(), std::mt19937(std::random_device()()));
            std::shuffle(c2.begin(), c2.end(), std::mt19937(std::random_device()()));
            ASSERT_THAT(CoordinationStructureClassification::calculateRmsdOfTwoStructs(c3, c2), DoubleEq(0));
        }
    }
}

TEST_F(TestCoordinationStructureClassification, calculateRmsdOfTwoStructsWithAllPermutation) {
    auto c2 = c1;
    do {
        ASSERT_THAT(CoordinationStructureClassification::calculateRmsdOfTwoStructs(c1, c2), DoubleEq(0));
    } while (std::next_permutation(c2.begin(), c2.end()));
}

TEST_F(TestCoordinationStructureClassification, calculateRmsdOfTwoStructsWithAllPermutationOfSmallShift) {
    auto c2 = c1;
    for (auto &coord : c2) {
        std::get<0>(coord) += 1;
    }
    boost::sort(c2);
    do {
        ASSERT_THAT(CoordinationStructureClassification::calculateRmsdOfTwoStructs(c1, c2),
                    DoubleNear(1.40435, 0.00001));
    } while (boost::next_permutation(c2));
}

TEST_F(TestCoordinationStructureClassification, calculateRmsdOfTwoStructsWithAllPermutationOfSmallShift2) {
    std::vector<std::tuple<double, double, double>> c2{
            {-1.91, 0.21,  -1.44},
            {2.92,  -1.35, 0.09},
            {0.73,  -0.39, 2.05},
            {0.03,  -1.66, -1.84},
            {-1.20, 0.13,  1.41},
            {-0.62, -2.24, 0.7},
            {-0.59, 2.5,   0.17},
            {1.09,  1.65,  0.54},
            {0.96,  1.08,  -2.16}
    };
    boost::sort(c2);
    do {
        ASSERT_THAT(CoordinationStructureClassification::calculateRmsdOfTwoStructs(c1, c2),
                    DoubleNear(0.758613, 0.000001));
    } while (boost::next_permutation(c2));
}

