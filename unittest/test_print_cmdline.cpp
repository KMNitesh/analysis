//
// Created by xiamr on 6/20/19.
//

#include "utils/common.hpp"
#include <gmock/gmock.h>
#include <type_traits>

using namespace testing;

TEST(print_comdline, NegativeArgc) {
    const char *const argv[] = {};
    int argc = -1;
    ASSERT_THAT(print_cmdline(argc, argv), Eq(""));
}

TEST(print_comdline, Empty) {
    const char *const argv[] = {};
    int argc = std::extent_v<decltype(argv)>;
    ASSERT_THAT(argc, Eq(0));
    ASSERT_THAT(print_cmdline(argc, argv), Eq(""));
}

TEST(print_comdline, NoArgument) {
    const char *const argv[] = {"analysis"};
    int argc = std::extent_v<decltype(argv)>;
    ASSERT_THAT(print_cmdline(argc, argv), Eq("analysis"));
}

TEST(print_comdline, WithOneOption) {
    const char *const argv[] = {"analysis", "-p", "Hello", "World"};
    int argc = std::extent_v<decltype(argv)>;
    ASSERT_THAT(print_cmdline(argc, argv), Eq("analysis -p Hello World"));
}

TEST(print_comdline, WithTwoOption) {
    const char *const argv[] = {"analysis", "-p", "Hello", "World", "-f", "OK"};
    int argc = std::extent_v<decltype(argv)>;
    ASSERT_THAT(print_cmdline(argc, argv), Eq("analysis -p Hello World -f OK"));
}
