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
    for (auto value : enumerate(iter)) {
        l.push_back(value);
    }
    ASSERT_THAT(l.size(), Eq(4u));
    auto[index, value] = l.front();
    ASSERT_TRUE(index == 0 && value == 1);
    l.pop_front();

    std::tie(index, value) = l.front();
    ASSERT_TRUE(index == 1 && value == 3);
    l.pop_front();

    std::tie(index, value) = l.front();
    ASSERT_TRUE(index == 2 && value == 4);
    l.pop_front();

    std::tie(index, value) = l.front();
    ASSERT_TRUE(index == 3 && value == 5);
    l.pop_front();

    ASSERT_TRUE(l.empty());
}

TEST(enumerate, NegativeIncrement) {
    std::list<std::pair<int, int>> l;
    auto iter = {1, 3, 4, 5};
    for (auto value : enumerate(iter).step(-2)) {
        l.push_back(value);
    }
    ASSERT_THAT(l.size(), Eq(4u));
    auto[index, value] = l.front();
    ASSERT_TRUE(index == 0 && value == 1);
    l.pop_front();

    std::tie(index, value) = l.front();
    ASSERT_TRUE(index == -2 && value == 3);
    l.pop_front();

    std::tie(index, value) = l.front();
    ASSERT_TRUE(index == -4 && value == 4);
    l.pop_front();

    std::tie(index, value) = l.front();
    ASSERT_TRUE(index == -6 && value == 5);
    l.pop_front();

    ASSERT_TRUE(l.empty());
}

TEST(enumerate, NegativeIncrementWithNoZeroStart) {
    std::list<std::pair<int, int>> l;
    auto iter = {1, 3, 4, 5};
    for (auto value : enumerate(iter).start(3).step(-2)) {
        l.push_back(value);
    }
    ASSERT_THAT(l.size(), Eq(4u));
    auto[index, value] = l.front();
    ASSERT_TRUE(index == 3 && value == 1);
    l.pop_front();

    std::tie(index, value) = l.front();
    ASSERT_TRUE(index == 1 && value == 3);
    l.pop_front();

    std::tie(index, value) = l.front();
    ASSERT_TRUE(index == -1 && value == 4);
    l.pop_front();

    std::tie(index, value) = l.front();
    ASSERT_TRUE(index == -3 && value == 5);
    l.pop_front();

    ASSERT_TRUE(l.empty());
}

TEST(enumerate, ZeroIncrement) {
    auto iter = {1, 3, 4, 5};
    ASSERT_THROW(enumerate(iter).step(0), std::runtime_error);
}