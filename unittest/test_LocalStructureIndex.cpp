//
// Created by xiamr on 9/23/19.
//

#include "std.hpp"
#include <gmock/gmock.h>
#include "LocalStructureIndex.hpp"

using namespace testing;

TEST(TestLocalStructureIndex, calculateLSIwithZeroLSI) {
    std::deque<double> rs{2.1, 2.9, 2.6, 2.4, 2.5, 2.7, 2.3, 2.8, 2.2};
    ASSERT_THAT(LocalStructureIndex::calculateLSI(rs), DoubleNear(0.0, 1e-10));
    ASSERT_THAT(rs.back(), DoubleEq(2.9));
}

TEST(TestLocalStructureIndex, calculateLSIwithNormalLSI) {
    std::deque<double> rs{2.1, 2.1, 2.3, 2.45, 2.65, 2.7, 2.3, 2.8, 2.2};
    ASSERT_THAT(LocalStructureIndex::calculateLSI(rs), DoubleNear(4.21875e-3, 1e-8));
    ASSERT_THAT(rs.back(), DoubleEq(2.8));
}





