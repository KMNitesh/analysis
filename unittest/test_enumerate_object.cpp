//
// Created by xiamr on 6/14/19.
//

#include <gmock/gmock.h>
#include <type_traits>
#include <string>
#include <list>
#include "../src/common.hpp"

using namespace testing;


TEST(enumerate, DefaultPostive) {
    std::list<std::pair<int, int>> l;
    auto iter = {1, 3, 4, 5};
    for (auto value : enumerate(iter).start(1)) {
        l.push_back(value);
    }
    ASSERT_THAT(l, ElementsAre(Pair(1, 1), Pair(2, 3), Pair(3, 4), Pair(4, 5)));
}

TEST(enumerate, NegativeIncrement) {
    std::list<std::pair<int, int>> l;
    auto iter = {1, 3, 4, 5};
    for (auto value : enumerate(iter).step(-2)) {
        l.push_back(value);
    }
    ASSERT_THAT(l, ElementsAre(Pair(0, 1), Pair(-2, 3), Pair(-4, 4), Pair(-6, 5)));
}

TEST(enumerate, NegativeIncrementWithNoZeroStart) {
    std::list<std::pair<int, int>> l;
    auto iter = {1, 3, 4, 5};
    for (auto value : enumerate(iter).start(3).step(-2)) {
        l.push_back(value);
    }
    ASSERT_THAT(l, ElementsAre(Pair(3, 1), Pair(1, 3), Pair(-1, 4), Pair(-3, 5)));
}

TEST(enumerate, ZeroIncrement) {
    auto iter = {1, 3, 4, 5};
    ASSERT_THROW(enumerate(iter).step(0), std::runtime_error);
}