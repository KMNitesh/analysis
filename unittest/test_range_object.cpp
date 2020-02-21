//
// Created by xiamr on 6/14/19.
//

#include <gmock/gmock.h>
#include <type_traits>
#include <string>
#include <list>
#include "utils/common.hpp"

using namespace testing;

TEST(range, WithThreeParameters) {
    std::list<int> l;
    for (auto i : range(-1, 6, 2)) {
        l.push_back(i);
    }
    ASSERT_THAT(l, ElementsAre(-1, 1, 3, 5));
}

TEST(range, WithTwoParameters) {
    std::list<int> l;
    for (auto i : range(1, 4)) {
        l.push_back(i);
    }
    ASSERT_THAT(l, ElementsAre(1, 2, 3));
}


TEST(range, WithSingleParameter) {
    std::list<int> l;
    for (auto i : range(4)) {
        l.push_back(i);
    }
    ASSERT_THAT(l, ElementsAre(0, 1, 2, 3));
}


TEST(range, WithNegativeIncrement) {
    std::list<int> l;
    for (auto i : range(4, -1, -2)) {
        l.push_back(i);
    }
    ASSERT_THAT(l, ElementsAre(4, 2, 0));
}

TEST(range, WithNoneResult) {
    for ([[maybe_unused]] auto i : range(10, 10)) {
        FAIL();
    }
}

TEST(range, ZeroIncrement) {
    ASSERT_THROW(range(10, 12, 0), std::runtime_error);
}

TEST(range, WrongDirectionNegative) {
    ASSERT_THROW(range(10, 12, -2), std::runtime_error);
}

TEST(range, WrongDirectionPositive) {
    ASSERT_THROW(range(12, 5, 2), std::runtime_error);
}


TEST(range, WithNegativeEnd) {
    ASSERT_THROW(range(-1), std::runtime_error);
}
