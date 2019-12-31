//
// Created by xiamr on 6/24/19.
//

#include <string>
#include <gmock/gmock.h>
#include <boost/range/adaptors.hpp>
#include "utils/common.hpp"

using namespace testing;
using namespace std;


TEST(PushIterable, SimpleValue) {
    std::vector<int> vec;
    auto push = {1, 2, 3, 4};
    for (auto v : PushIterable(push)) {
        vec.push_back(v);
    }
    ASSERT_THAT(vec, ElementsAre(1, 2, 3, 4));
}

TEST(PushIterable, SimpleValueWithPushBack) {
    std::vector<int> vec;
    auto push = {1, 2, 3, 4};
    auto it = PushIterable(push);
    it.push_back(5);
    for (auto v : it) {
        vec.push_back(v);
    }
    ASSERT_THAT(vec, ElementsAre(5, 1, 2, 3, 4));
}

TEST(combine_seq, Simple) {
    std::vector<std::string> vec;
    auto push = {1, 5, 9, 11, 12, 13, 15, 16};
    for (auto v : combine_seq(push)) {
        vec.push_back(v);
    }
    ASSERT_THAT(vec, ContainerEq(std::vector<std::string>{"1", "5", "9", "11-13", "15", "16"}));
}

TEST(combine_seq, SingleValue) {
    std::vector<std::string> vec;
    auto push = {5};
    for (auto v : combine_seq(push)) {
        vec.push_back(v);
    }
    ASSERT_THAT(vec, ContainerEq(std::vector<std::string>{"5"}));
}

TEST(combine_seq, Nagitve) {
    std::vector<std::string> vec;
    auto push = {-1, -3};
    for (auto v : combine_seq(push)) {
        vec.push_back(v);
    }
    ASSERT_THAT(vec, ContainerEq(std::vector<std::string>{"-1", "-3"}));
}

TEST(combine_seq, Empty) {
    std::vector<std::string> vec;
    std::vector<int> seq;
    for (auto v : combine_seq(seq)) {
        vec.push_back(v);
    }
    ASSERT_THAT(vec, ContainerEq(std::vector<std::string>{}));
}

TEST(combine_seq, STL_Vector) {
    std::vector<std::string> vec;
    std::vector<int> seq{1, 2, 3, 4, 5, 6, 7};
    for (auto v : combine_seq(seq)) {
        vec.push_back(v);
    }
    ASSERT_THAT(vec, ContainerEq(std::vector<std::string>{"1-7"}));
}


TEST(combine_seq, RealProblem) {
    std::vector<std::string> vec;
    auto push = {45, 46, 57, 58, 59, 60};
    for (auto v : combine_seq(push)) {
        vec.push_back(v);
    }
    ASSERT_THAT(vec, ContainerEq(std::vector<std::string>{"45", "46", "57-60"}));
}

TEST(combine_seq, Tradition) {
    std::vector<std::string> vec;
    auto push = {45, 46, 57, 58, 59, 60};
    auto seq = combine_seq(push);
    auto it = seq.begin();
    for (; it != seq.end(); ++it) {
        vec.push_back(*it);
    }
    ASSERT_THAT(vec, ContainerEq(std::vector<std::string>{"45", "46", "57-60"}));
}

TEST(combine_seq, TraditionReference) {
    std::vector<std::string> vec;
    auto push = {45, 46, 57, 58, 59, 60};
    auto seq = combine_seq(push);
    auto &it = seq.begin();
    for (; it != seq.end(); ++it) {
        vec.push_back(*it);
    }
    ASSERT_THAT(vec, ContainerEq(std::vector<std::string>{"45", "46", "57-60"}));
}

