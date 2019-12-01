//
// Created by xiamr on 6/14/19.
//

#include <gmock/gmock.h>
#include <type_traits>
#include <string>
#include <list>
#include "utils/common.hpp"


using namespace testing;

TEST(format, SingleIntger) {
    ASSERT_EQ(format("%d", 2).str(), (boost::format("%d") % 2).str());
}

TEST(format, SingleStr) {
    ASSERT_EQ(format("%s", "Hello").str(), (boost::format("%s") % "Hello").str());
}

TEST(format, TwoIntgers) {
    ASSERT_EQ(format("%d %d", 2, 3).str(), (boost::format("%d %d") % 2 % 3).str());
}

TEST(format, IntergerWithAdd) {
    ASSERT_EQ(format("%d ", 2 + 3).str(), (boost::format("%d ") % (2 + 3)).str());
    ASSERT_EQ(format("%d %d ", 2 + 3, 5 + 6).str(), (boost::format("%d %d ") % (2 + 3) % (5 + 6)).str());
}

