//
// Created by xiamr on 8/30/19.
//

#include "utils/std.hpp"
#include <gmock/gmock.h>
#include "others/CrossCorrelation.hpp"

using namespace testing;


TEST(CrossCorrelation, PostiveRelationForCalculateCrossCorrelationDWang) {
    std::deque<std::tuple<double, double, double>> series = {
            {1,  1,  11},
            {2,  2,  12},
            {3,  3,  13},
            {4,  4,  14},
            {5,  5,  15},
            {6,  6,  16},
            {7,  7,  17},
            {8,  8,  18},
            {9,  9,  19},
            {10, 10, 20}
    };

    ASSERT_THAT(CrossCorrelation::calculateCrossCorrelationDWang(series), Gt(0));
    ASSERT_THAT(CrossCorrelation::calculateCrossCorrelationDWang2(series), Gt(0));
}


TEST(CrossCorrelation, NegativeRelationForCalculateCrossCorrelationDWang) {
    std::deque<std::tuple<double, double, double>> series = {
            {1,  1,  20},
            {2,  2,  19},
            {3,  3,  18},
            {4,  4,  17},
            {5,  5,  16},
            {6,  6,  15},
            {7,  7,  14},
            {8,  8,  13},
            {9,  9,  12},
            {10, 10, 11}
    };

    ASSERT_THAT(CrossCorrelation::calculateCrossCorrelationDWang(series), Lt(0));
    ASSERT_THAT(CrossCorrelation::calculateCrossCorrelationDWang2(series), Lt(0));
}


