//
// Created by xiamr on 7/24/19.
//

#include <gmock/gmock.h>
#include "TinkerDynReader.hpp"

using namespace std;
using namespace testing;


TEST(TinkerDynReaderTest, readContent) {
    auto inputstream = make_shared<std::istringstream>(
            " Number of Atoms and Title :\n"
            "  1539  BUILT WITH PACKMOL\n"
            " Periodic Box Dimensions :\n"
            "    0.2484379303840836D+02    0.2484379303840836D+02    0.2484379303840836D+02\n"
            "    0.9000000000000000D+02    0.9000000000000000D+02    0.9000000000000000D+02\n"
            " Current Atomic Positions :\n"
            "   -0.4961401991829424D+01   -0.7444276104052326D+01    0.1154951382088674D+02\n"
            "   -0.4539746555161510D+01   -0.7868545106902119D+01   -0.7868545106902119D+01\n"
            "   -0.5179750544034674D+01   -0.6483757363157239D+01    0.1127877357752338D+02\n"
            " Current Atomic Velocities :\n"
            "    0.9277195775170808D+01    0.2457196471429321D+01   -0.2009712267283242D+01\n"
    );

    TinkerDynReader reader(inputstream);
    reader.readContent();
    ASSERT_THAT(reader.getPBCBox().xbox, DoubleEq(0.2484379303840836E+02));
    ASSERT_THAT(reader.getPBCBox().ybox, DoubleEq(0.2484379303840836E+02));
    ASSERT_THAT(reader.getPBCBox().zbox, DoubleEq(0.2484379303840836E+02));

    ASSERT_THAT(reader.getPBCBox().alpha, DoubleEq(90.0));
    ASSERT_THAT(reader.getPBCBox().beta, DoubleEq(90.0));
    ASSERT_THAT(reader.getPBCBox().gamma, DoubleEq(90.0));

    ASSERT_THAT(reader.getAtomicPosition().coordinates.size(), Eq(3));

    //ATOM1
    ASSERT_THAT(get<0>(reader.getAtomicPosition().coordinates.at(0)), DoubleEq(-0.4961401991829424E+01));
    ASSERT_THAT(get<1>(reader.getAtomicPosition().coordinates.at(0)), DoubleEq(-0.7444276104052326E+01));
    ASSERT_THAT(get<2>(reader.getAtomicPosition().coordinates.at(0)), DoubleEq(0.1154951382088674E+02));
    //ATOM2
    ASSERT_THAT(get<0>(reader.getAtomicPosition().coordinates.at(1)), DoubleEq(-0.4539746555161510E+01));
    ASSERT_THAT(get<1>(reader.getAtomicPosition().coordinates.at(1)), DoubleEq(-0.7868545106902119E+01));
    ASSERT_THAT(get<2>(reader.getAtomicPosition().coordinates.at(1)), DoubleEq(-0.7868545106902119E+01));
    //ATOM3
    ASSERT_THAT(get<0>(reader.getAtomicPosition().coordinates.at(2)), DoubleEq(-0.5179750544034674E+01));
    ASSERT_THAT(get<1>(reader.getAtomicPosition().coordinates.at(2)), DoubleEq(-0.6483757363157239E+01));
    ASSERT_THAT(get<2>(reader.getAtomicPosition().coordinates.at(2)), DoubleEq(0.1127877357752338E+02));
}

