//
// Created by xiamr on 11/1/19.
//

#include "utils/std.hpp"
#include <gmock/gmock.h>
#include <boost/core/ignore_unused.hpp>
#include "others/NBOSpin.hpp"

using namespace testing;

TEST(NBOSpinTest, OneLine) {

    std::string line = "      B  1      [core]2s( 0.44)2p( 1.32)3p( 0.01)3d( 0.01)";

    ASSERT_THAT(NBOSpin::total_spin(line), DoubleEq(1.78));
}

TEST(NBOSpinTest, ThrowExceptionWhenWithoutRightParentheses) {

    std::string line = "      B  1      [core]2s( 0.44)2p( 1.32)3p( 0.01)3d( 0.01";

    ASSERT_THROW(boost::ignore_unused(NBOSpin::total_spin(line)), std::runtime_error);
}

TEST(NBOSpinTest, ThrowExceptionWhenWithoutMiddleParentheses) {

    std::string line = "      B  1      [core]2s( 0.442p( 1.32)3p( 0.01)3d( 0.01)";

    ASSERT_ANY_THROW(boost::ignore_unused(NBOSpin::total_spin(line)));
}

TEST(NBOSpinTest, WholeTable) {
    std::stringstream ss("                                 Natural Population\n"
                         " ---------------------------------------------------------\n"
                         "   Core                      33.99063 ( 99.9724% of   34)\n"
                         "   Valence                   92.72161 ( 99.7007% of   93)\n"
                         "   Natural Minimal Basis    126.71224 ( 99.7734% of  127)\n"
                         "   Natural Rydberg Basis      0.28776 (  0.2266% of  127)\n"
                         " ---------------------------------------------------------\n"
                         "\n"
                         "    Atom No         Natural Electron Configuration\n"
                         " ----------------------------------------------------------------------------\n"
                         "      B  1      [core]2s( 0.44)2p( 1.32)3p( 0.01)3d( 0.01)\n"
                         "      O  2      [core]2s( 1.67)2p( 5.22)3d( 0.01)\n"
                         "\n"
                         "    Atom No         Natural Electron Configuration\n"
                         " ----------------------------------------------------------------------------\n"
                         "      B  1      [core]2s( 0.22)2p( 0.66)3p( 0.01)\n"
                         "      O  2      [core]2s( 0.83)2p( 2.61)3d( 0.01)\n"
                         "\n"
                         "    Atom No         Natural Electron Configuration\n"
                         " ----------------------------------------------------------------------------\n"
                         "      B  1      [core]2s( 0.22)2p( 0.66)3p( 0.01)\n"
                         "      O  2      [core]2s( 0.83)2p( 2.61)3d( 0.01)\n"
                         "\n");

    auto table = NBOSpin::getElectronSpin(ss);
    ASSERT_THAT(table.size(), Eq(2));
    ASSERT_THAT(table.at(1).first, Eq("B"));
    ASSERT_THAT(table.at(1).second, Pointwise(DoubleEq(), std::array<double, 3>{1.78, 0.89, 0.89}));

    ASSERT_THAT(table.at(2).first, Eq("O"));
    ASSERT_THAT(table.at(2).second, Pointwise(DoubleEq(), std::array<double, 3>{6.90, 3.45, 3.45}));
}