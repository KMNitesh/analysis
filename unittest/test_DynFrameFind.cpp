//
// Created by xiamr on 7/24/19.
//

#include <gmock/gmock.h>
#include "ana_module/DynFrameFind.hpp"
#include "data_structure/frame.hpp"
#include "data_structure/atom.hpp"

using namespace std;
using namespace testing;

class DynFrameFindTest : public DynFrameFind, public Test {
protected:
    virtual void SetUp() {
        reader = make_shared<TinkerDynReader>(make_shared<std::istringstream>(
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
        ));
        eps = 1E-7;
        reader->readContent();

        // Matched Frame
        frame = make_shared<Frame>();
        frame->box = PBCBox(0.2484379303840836E+02, 0.2484379303840836E+02, 0.2484379303840836E+02, 90.0, 90.0, 90.0);

        atom1 = make_shared<Atom>();
        atom1->x = -0.4961401991829424E+01;
        atom1->y = -0.7444276104052326E+01;
        atom1->z = 0.1154951382088674E+02;
        frame->atom_list.push_back(atom1);

        atom2 = make_shared<Atom>();
        atom2->x = -0.4539746555161510E+01;
        atom2->y = -0.7868545106902119E+01;
        atom2->z = -0.7868545106902119E+01;
        frame->atom_list.push_back(atom2);

        atom3 = make_shared<Atom>();
        atom3->x = -0.5179750544034674E+01;
        atom3->y = -0.6483757363157239E+01;
        atom3->z = 0.1127877357752338E+02;
        frame->atom_list.push_back(atom3);
    }

    shared_ptr<Atom> atom1, atom2, atom3;
    shared_ptr<Frame> frame;
    static constexpr auto AnyFrameNumber = 10;
};

TEST_F(DynFrameFindTest, processFrameMatch) {
    nframe = AnyFrameNumber;
    process(frame);
    ASSERT_THAT(matched_frames.size(), Eq(1));
    ASSERT_THAT(matched_frames.front(), Eq(AnyFrameNumber + 1));
}

TEST_F(DynFrameFindTest, processFrameNotMatchWithAtomPosition) {
    nframe = AnyFrameNumber;
    atom1->x = 10.00;
    process(frame);
    ASSERT_THAT(matched_frames.size(), Eq(0));
}

TEST_F(DynFrameFindTest, processFrameNotMatchWithBoxLength) {
    nframe = AnyFrameNumber;
    frame->box = PBCBox(0.2484379303840836E+02, 9.99, 0.2484379303840836E+02, 90.0, 90.0, 90.0);
    process(frame);
    ASSERT_THAT(matched_frames.size(), Eq(0));
}

TEST_F(DynFrameFindTest, processFrameNotMatchWithBoxAngle) {
    nframe = AnyFrameNumber;
    frame->box = PBCBox(0.2484379303840836E+02, 0.2484379303840836E+02, 0.2484379303840836E+02, 90.0, 89.9999, 90.0);
    process(frame);
    ASSERT_THAT(matched_frames.size(), Eq(0));
}





