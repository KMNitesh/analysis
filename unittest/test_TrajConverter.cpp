
#include <gmock/gmock.h>
#include <boost/fusion/sequence.hpp>
#include "others/TrajConverter.hpp"

using namespace testing;

TEST(TestTrajConverter, parseNormalHeaderLine) {
    ASSERT_THAT(TrajConverter::parse_header("block=coordinates records=14"), Eq(14));
}

TEST(TestTrajConverter, parseNormalAtomLine) {
    auto attr = TrajConverter::parse_atom("B 1.024232 2.539793 0.000000");

    ASSERT_THAT(boost::fusion::at_c<0>(attr), Eq("B"));
    ASSERT_THAT(boost::fusion::at_c<1>(attr), DoubleEq(1.024232));
    ASSERT_THAT(boost::fusion::at_c<2>(attr), DoubleEq(2.539793));
    ASSERT_THAT(boost::fusion::at_c<3>(attr), DoubleEq(0.000000));
}