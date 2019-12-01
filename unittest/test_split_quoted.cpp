//
// Created by xiamr on 7/18/19.
//

#include <gmock/gmock.h>

#include "utils/common.hpp"

using namespace std;
using namespace testing;

TEST(split_quoted, NoQuoted) {
    ASSERT_THAT(split_quoted(" aa  bb 11 "), ContainerEq(vector<string>{"aa", "bb", "11"}));
}

TEST(split_quoted, SingleQuoted) {
    ASSERT_THAT(split_quoted(R"( aa  bb 11 " Hello " 22.1 )"),
                ContainerEq(vector<string>{"aa", "bb", "11", " Hello ", "22.1"}));
}

TEST(split_quoted, DoubleQuoted) {
    ASSERT_THAT(split_quoted(R"( aa  bb 11 " Hello " 22.1 "OK  OK" )"),
                ContainerEq(vector<string>{"aa", "bb", "11", " Hello ", "22.1", "OK  OK"}));
}