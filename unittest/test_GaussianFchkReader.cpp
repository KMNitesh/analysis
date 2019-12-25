
#include <gmock/gmock.h>
#include <sstream>
#include <string>
#include "others/GaussianFchkReader.hpp"

using namespace testing;

TEST(TEST_GaussianFchkReader, ReaderAtomCartesian) {

    std::stringstream ss("  1.00000000E+00  1.00000000E+00  8.00000000E+00  1.00000000E+00  1.00000000E+00\n"
                         "Current cartesian coordinates              R   N=          75\n"
                         "  9.49651332E-03  2.47881007E-02  2.35123949E-01 -3.85786281E+00  1.37651239E-01\n"
                         " -2.51052097E+00 -3.86024897E+00 -2.65805097E-02 -4.34131488E+00 -4.91199818E+00\n"
                         "  1.58782973E+00 -2.08954118E+00  3.74426028E+00  2.94034545E+00  1.65438207E+00\n"
                         "  4.13499513E+00  3.34046593E+00  3.40501762E+00  5.31640482E+00  2.39586646E+00\n"
                         "  8.69214325E-01 -3.56127551E+00  3.07904527E+00  1.54746993E+00 -4.28154082E+00\n"
                         "  3.10767207E+00  3.23880460E+00 -3.04833199E+00  4.80723340E+00  1.17226267E+00\n"
                         " -3.48827134E+00 -2.89330299E+00  1.77203295E+00 -5.15882210E+00 -2.55477556E+00\n"
                         "  1.08140291E+00 -3.73097405E+00 -3.33971973E+00  3.53897706E+00  2.58000082E-01\n"
                         "  4.29827368E+00 -1.80839382E+00  1.81431981E+00  5.09767734E+00 -1.22387875E+00\n"
                         "  1.17344295E-01  4.63518207E+00 -3.61032078E+00  3.99082455E+00 -3.03316248E-01\n"
                         " -2.35186859E+00  4.08035314E+00 -1.43811693E-02 -4.16567474E+00  4.92528738E+00\n"
                         " -1.85185907E+00 -2.00988278E+00 -2.34982882E-01 -4.05389494E+00 -2.17747216E+00\n"
                         " -1.79129140E+00 -4.97295278E+00 -1.82969151E+00  9.05540943E-02 -4.22838383E+00\n"
                         " -3.97855360E+00  3.05434661E+00 -3.41155621E+00  1.74402106E+00  3.59551950E+00\n"
                         " -3.59833923E+00  3.49104885E+00  2.56594867E+00 -5.08576664E+00  1.15815131E+00\n"
                         "Force Field                                I                0\n"
                         "Int Atom Types                             I   N=          25\n");


    std::vector<double> expected{
            9.49651332E-03, 2.47881007E-02, 2.35123949E-01, -3.85786281E+00, 1.37651239E-01,
            -2.51052097E+00, -3.86024897E+00, -2.65805097E-02, -4.34131488E+00, -4.91199818E+00,
            1.58782973E+00, -2.08954118E+00, 3.74426028E+00, 2.94034545E+00, 1.65438207E+00,
            4.13499513E+00, 3.34046593E+00, 3.40501762E+00, 5.31640482E+00, 2.39586646E+00,
            8.69214325E-01, -3.56127551E+00, 3.07904527E+00, 1.54746993E+00, -4.28154082E+00,
            3.10767207E+00, 3.23880460E+00, -3.04833199E+00, 4.80723340E+00, 1.17226267E+00,
            -3.48827134E+00, -2.89330299E+00, 1.77203295E+00, -5.15882210E+00, -2.55477556E+00,
            1.08140291E+00, -3.73097405E+00, -3.33971973E+00, 3.53897706E+00, 2.58000082E-01,
            4.29827368E+00, -1.80839382E+00, 1.81431981E+00, 5.09767734E+00, -1.22387875E+00,
            1.17344295E-01, 4.63518207E+00, -3.61032078E+00, 3.99082455E+00, -3.03316248E-01,
            -2.35186859E+00, 4.08035314E+00, -1.43811693E-02, -4.16567474E+00, 4.92528738E+00,
            -1.85185907E+00, -2.00988278E+00, -2.34982882E-01, -4.05389494E+00, -2.17747216E+00,
            -1.79129140E+00, -4.97295278E+00, -1.82969151E+00, 9.05540943E-02, -4.22838383E+00,
            -3.97855360E+00, 3.05434661E+00, -3.41155621E+00, 1.74402106E+00, 3.59551950E+00,
            -3.59833923E+00, 3.49104885E+00, 2.56594867E+00, -5.08576664E+00, 1.15815131E+00
    };

    auto cartesians = GaussianFchkReader::readCartesian(ss);

    ASSERT_THAT(cartesians, Pointwise(DoubleEq(), expected));


    auto coord = GaussianFchkReader::getCoordinateOfAtom(1, cartesians);
    ASSERT_THAT(std::get<0>(coord), DoubleEq(9.49651332E-03));
    ASSERT_THAT(std::get<1>(coord), DoubleEq(2.47881007E-02));
    ASSERT_THAT(std::get<2>(coord), DoubleEq(2.35123949E-01));

    ASSERT_NO_THROW(GaussianFchkReader::getCoordinateOfAtom(25, cartesians));
    ASSERT_ANY_THROW(GaussianFchkReader::getCoordinateOfAtom(26, cartesians));
}