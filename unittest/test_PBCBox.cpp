
#include "utils/std.hpp"
#include <gmock/gmock.h>
#include "utils/PBCBox.hpp"

using namespace testing;

TEST(PBCBoxTest, orthogonal) {
    PBCBox box(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
    ASSERT_THAT(box.get_box_type(), Eq(PBCBox::Type::orthogonal));
}

TEST(PBCBoxTest, get_orthogonal) {
    PBCBox box(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
    ASSERT_THAT(box.getBoxParameter(), Pointwise(DoubleEq(), {10.0, 10.0, 10.0, 90.0, 90.0, 90.0}));
}

TEST(PBCBoxTest, orthogonal_to_box) {
    PBCBox pbc_box(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
    gmx::matrix box;
    pbc_box.getBoxParameter(box);

    ASSERT_THAT(box[0][0], FloatNear(1.0, 1E-10));
    ASSERT_THAT(box[0][1], FloatNear(0.0, 1E-10));
    ASSERT_THAT(box[0][2], FloatNear(0.0, 1E-10));

    ASSERT_THAT(box[1][0], FloatNear(0.0, 1E-10));
    ASSERT_THAT(box[1][1], FloatNear(1.0, 1E-10));
    ASSERT_THAT(box[1][2], FloatNear(0.0, 1E-10));

    ASSERT_THAT(box[2][0], FloatNear(0.0, 1E-10));
    ASSERT_THAT(box[2][1], FloatNear(0.0, 1E-10));
    ASSERT_THAT(box[2][2], FloatNear(1.0, 1E-10));
}

TEST(PBCBoxTest, box_constructor_orthogonal) {

    gmx::matrix box{
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0}
    };
    PBCBox pbc_box(box);
    ASSERT_THAT(pbc_box.get_box_type(), Eq(PBCBox::Type::orthogonal));

    ASSERT_THAT(pbc_box.getBoxParameter(), Pointwise(DoubleNear(1E-5), {10.0, 10.0, 10.0, 90.0, 90.0, 90.0}));
}

TEST(PBCBoxTest, box_constructor_octahedron) {
    gmx::matrix box{
            {1.61941e+01,  0.00000e+00, 0.00000e+00},
            {5.39803e+00,  1.52679e+01, 0.00000e+00},
            {-5.39803e+00, 7.63396e+00, 1.32224e+01}
    };

    PBCBox pbc_box(box);
    ASSERT_THAT(pbc_box.get_box_type(), Eq(PBCBox::Type::octahedron));

    ASSERT_THAT(pbc_box.getBoxParameter(), Pointwise(DoubleNear(1E-2),
                                                     {161.940994, 161.940582, 161.940704, 70.528800, 109.471246,
                                                      70.528738}));
}

TEST(PBCBoxTest, octahedron) {

    PBCBox pbc_box(161.940994, 161.940582, 161.940704, 70.528800, 109.471246, 70.528738);
    ASSERT_THAT(pbc_box.get_box_type(), Eq(PBCBox::Type::octahedron));

    gmx::matrix box;
    pbc_box.getBoxParameter(box);

    ASSERT_THAT(box[0][0], FloatNear(1.61941e+01, 1E-3));
    ASSERT_THAT(box[0][1], FloatNear(0.00000e+00, 1E-3));
    ASSERT_THAT(box[0][2], FloatNear(0.00000e+00, 1E-3));

    ASSERT_THAT(box[1][0], FloatNear(5.39803e+00, 1E-3));
    ASSERT_THAT(box[1][1], FloatNear(1.52679e+01, 1E-3));
    ASSERT_THAT(box[1][2], FloatNear(0.00000e+00, 1E-3));

    ASSERT_THAT(box[2][0], FloatNear(-5.39803e+00, 1E-3));
    ASSERT_THAT(box[2][1], FloatNear(7.63396e+00, 1E-3));
    ASSERT_THAT(box[2][2], FloatNear(1.32224e+01, 1E-3));
}
