//
// Created by xiamr on 6/14/19.
//

#include <gmock/gmock.h>
#include <type_traits>
#include <string>
#include <list>
#include "../src/common.hpp"

using namespace testing;

TEST(range, WithThreeParameters) {
    std::list<int> l;
    for (auto i : range(-1, 6, 2)) {
        l.push_back(i);
    }
    ASSERT_TRUE(l.size() == 4);
    
    int i;
    i = l.front();
    ASSERT_TRUE(i == -1);
    l.pop_front();
    i = l.front();
    ASSERT_TRUE(i == 1);
    l.pop_front();
    i = l.front();
    ASSERT_TRUE(i == 3);
    l.pop_front();
    i = l.front();
    ASSERT_TRUE(i == 5);
    l.pop_front();

    ASSERT_TRUE(l.empty());
}

TEST(range, WithTwoParameters) {
    std::list<int> l;
    for (auto i : range(1, 4)) {
        l.push_back(i);
    }
    ASSERT_TRUE(l.size() == 3);
    
    int i;
    i = l.front();
    ASSERT_TRUE(i == 1);
    l.pop_front();
    i = l.front();
    ASSERT_TRUE(i == 2);
    l.pop_front();
    i = l.front();
    ASSERT_TRUE(i == 3);
    l.pop_front();
    ASSERT_TRUE(l.empty());
}


TEST(range, WithSingleParameter) {
    std::list<int> l;
    for (auto i : range(4)) {
        l.push_back(i);
    }
    ASSERT_TRUE(l.size() == 4);

    int i;
    i = l.front();
    ASSERT_TRUE(i == 0);
    l.pop_front();
    i = l.front();
    ASSERT_TRUE(i == 1);
    l.pop_front();
    i = l.front();
    ASSERT_TRUE(i == 2);
    l.pop_front();
    i = l.front();
    ASSERT_TRUE(i == 3);
    l.pop_front();

    ASSERT_TRUE(l.empty());
}


TEST(range, WithNegativeIncrement) {
    std::list<int> l;
    for (auto i : range(4, -1, -2)) {
        l.push_back(i);
    }
    ASSERT_TRUE(l.size() == 3);

    int i;
    i = l.front();
    ASSERT_TRUE(i == 4);
    l.pop_front();
    i = l.front();
    ASSERT_TRUE(i == 2);
    l.pop_front();
    i = l.front();
    ASSERT_TRUE(i == 0);
    l.pop_front();

    ASSERT_TRUE(l.empty());
}

TEST(range, WithNoneResult) {
    for (auto i : range(10, 10)) {
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