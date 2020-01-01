
#include <gmock/gmock.h>
#include "others/NBOOrbitalComposition.hpp"
#include <boost/fusion/sequence.hpp>

using namespace testing;


TEST(TestNBOOrbitalComposition, RegularExpressionMatch) {
    ASSERT_THAT(NBOOrbitalComposition::match("   1. (1.00000) CR ( 1)Am  1            s(100.00%)"), Eq(true));
}

TEST(TestNBOOrbitalComposition, SpiritParseLine) {
    namespace bf = boost::fusion;
    std::optional<bf::vector<int, std::string, std::string, std::string, double>> attribute
            = NBOOrbitalComposition::parseLine("       1    1(Am)    s        Cor(5s)      0.000000%");

    ASSERT_THAT(attribute.has_value(), Eq(true));

    ASSERT_THAT(bf::at_c<0>(*attribute), Eq(1));
    ASSERT_THAT(bf::at_c<1>(*attribute), Eq("Am"));
    ASSERT_THAT(bf::at_c<2>(*attribute), Eq("Cor"));
    ASSERT_THAT(bf::at_c<3>(*attribute), Eq("5s"));
    ASSERT_THAT(bf::at_c<4>(*attribute), DoubleEq(0.0));
}

TEST(TestNBOOrbitalComposition, SpiritParseLineFail) {

    auto attribute = NBOOrbitalComposition::parseLine("    NAO#   Center   Label      Type      Composition");

    ASSERT_THAT(attribute.has_value(), Eq(false));
}