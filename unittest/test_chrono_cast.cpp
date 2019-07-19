//
// Created by xiamr on 7/19/19.
//

#include <gmock/gmock.h>
#include <chrono>
#include "common.hpp"

using namespace std;
using namespace testing;

TEST(chrono_cast, ZeroSecond) {
    ASSERT_THAT(chrono_cast(chrono::seconds(0)), Eq("0 sec"));
}

TEST(chrono_cast, TenSecs) {
    ASSERT_THAT(chrono_cast(chrono::seconds(10)), StrEq("10 secs"));
}

TEST(chrono_cast, ExactlyOneMinite) {
    ASSERT_THAT(chrono_cast(chrono::seconds(60)), StrEq("1 min"));
}

TEST(chrono_cast, ExactlyTwoMinites) {
    ASSERT_THAT(chrono_cast(chrono::seconds(120)), StrEq("2 mins"));
}

TEST(chrono_cast, SeventSeconds) {
    ASSERT_THAT(chrono_cast(chrono::seconds(70)), StrEq("1 min 10 secs"));
}

TEST(chrono_cast, OneHour) {
    ASSERT_THAT(chrono_cast(chrono::hours(1)), StrEq("1 hour"));
}

TEST(chrono_cast, TwoHours) {
    ASSERT_THAT(chrono_cast(chrono::hours(2)), StrEq("2 hours"));
}

TEST(chrono_cast, OneHourOneMinteOneSecond) {
    ASSERT_THAT(chrono_cast(chrono::hours(1) + chrono::minutes(1) + chrono::seconds(1)), StrEq("1 hour 1 min 1 sec"));
}

TEST(chrono_cast, TwoHourThreeMinteFourSecond) {
    ASSERT_THAT(chrono_cast(chrono::hours(2) + chrono::minutes(3) + chrono::seconds(4)),
                StrEq("2 hours 3 mins 4 secs"));
}





