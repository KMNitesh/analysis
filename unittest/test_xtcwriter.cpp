//
// Created by xiamr on 6/28/19.
//

#include <gmock/gmock.h>

#include "utils/GromacsInterface.hpp"
#include "utils/ThrowAssert.hpp"
#include "data_structure/atom.hpp"
#include "data_structure/frame.hpp"
#include "gtest_utility.hpp"
#include "utils/xtc_writer.hpp"
#include "GromacsInterfaceMock.hpp"

using namespace testing;
using namespace std;


class XTCWriterTest : public XTCWriter {
protected:
    GromacsInterface *getGromacsImpl() override { return gromacsMock; }

public:
    GromacsMock *gromacsMock;

    XTCWriterTest(GromacsMock *gromacsMock) : gromacsMock(gromacsMock) {}
};


TEST(XTCWriter, OpenAndCloseFile) {

    GromacsMock mock;
    XTCWriterTest writer(&mock);

    std::string filename = "output.trr";
    auto fio = reinterpret_cast<gmx::t_fileio *>(1000);
    EXPECT_CALL(mock, open_xtc(StrEq(filename), StrEq("w"))).WillOnce(Return(fio));
    EXPECT_CALL(mock, close_xtc(fio)).Times(1);
    writer.open(filename);
    writer.close();
}

TEST(XTCWriter, WriteFrameToFile) {

    GromacsMock mock;
    XTCWriterTest writer(&mock);

    std::string filename = "output.trr";
    auto fio = reinterpret_cast<gmx::t_fileio *>(1000);
    auto frame = make_shared<Frame>();
    frame->box = PBCBox(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);

    auto atom1 = make_shared<Atom>();

    atom1->x = -9.90;
    atom1->y = 8.89;
    atom1->z = 2.1;
    atom1->seq = 1;

    auto atom2 = make_shared<Atom>();
    atom2->x = 0.90;
    atom2->y = 3.19;
    atom2->z = 4.1;
    atom2->seq = 2;

    frame->atom_list = {atom1, atom2};
    frame->atom_map = {{atom1->seq, atom1},
                       {atom2->seq, atom2}};

    gmx::rvec box[3] = {{1.0, 0.0, 0.0},
                        {0.0, 1.0, 0.0},
                        {0.0, 0.0, 1.0}};
    gmx::rvec x[2] = {{-0.99, 0.889, 0.21},
                      {0.09,  0.319, 0.41}};

    EXPECT_CALL(mock, open_xtc(StrEq(filename), StrEq("w"))).WillOnce(Return(fio));
    EXPECT_CALL(mock, write_xtc(fio, 2, 0, 0.0, FLOAT_ARRAY_EQ(box, 9), FLOAT_ARRAY_EQ(x, 6), 1000)).Times(1);
    EXPECT_CALL(mock, write_xtc(fio, 2, 1, 1.0, FLOAT_ARRAY_EQ(box, 9), FLOAT_ARRAY_EQ(x, 6), 1000)).Times(1);
    EXPECT_CALL(mock, close_xtc(fio)).Times(1);

    writer.open(filename);
    writer.write(frame);
    writer.write(frame);
    writer.close();
}

TEST(XTCWriter, CloseBeforeOpen) {
    GromacsMock mock;
    XTCWriterTest writer(&mock);
    ASSERT_THROW(writer.close(), AssertionFailureException);
}

TEST(XTCWriter, WriteBeforeOpen) {
    GromacsMock mock;
    XTCWriterTest writer(&mock);
    auto frame = make_shared<Frame>();
    ASSERT_THROW(writer.write(frame), AssertionFailureException);
}

TEST(XTCWriter, WriteAfterClose) {
    GromacsMock mock;
    XTCWriterTest writer(&mock);
    auto frame = make_shared<Frame>();
    EXPECT_CALL(mock, open_xtc(StrEq("Output.xtc"), StrEq("w"))).WillOnce(
            Return(reinterpret_cast<gmx::t_fileio *>(10100)));
    EXPECT_CALL(mock, close_xtc(reinterpret_cast<gmx::t_fileio *>(10100))).Times(1);
    ASSERT_NO_THROW(writer.open("Output.xtc"));
    ASSERT_NO_THROW(writer.close());
    ASSERT_THROW(writer.write(frame), AssertionFailureException);
}

TEST(XTCWriter, OpenError) {
    GromacsMock mock;
    XTCWriterTest writer(&mock);
    EXPECT_CALL(mock, open_xtc(StrEq("Output.xtc"), StrEq("w"))).WillOnce(Return(nullptr));
    ASSERT_THROW(writer.open("Output.xtc"), AssertionFailureException);
}

TEST(XTCWriter, OpenTwice) {
    GromacsMock mock;
    XTCWriterTest writer(&mock);
    EXPECT_CALL(mock, open_xtc(StrEq("Output.xtc"), StrEq("w"))).WillOnce(
            Return(reinterpret_cast<gmx::t_fileio *>(10100)));
    ASSERT_NO_THROW(writer.open("Output.xtc"));
    ASSERT_THROW(writer.open("Output.xtc"), AssertionFailureException);
}

TEST(XTCWriter, CloseTwice) {
    GromacsMock mock;
    XTCWriterTest writer(&mock);
    EXPECT_CALL(mock, open_xtc(StrEq("Output.xtc"), StrEq("w"))).WillOnce(
            Return(reinterpret_cast<gmx::t_fileio *>(10100)));
    EXPECT_CALL(mock, close_xtc(reinterpret_cast<gmx::t_fileio *>(10100))).Times(1);
    ASSERT_NO_THROW(writer.open("Output.xtc"));
    ASSERT_NO_THROW(writer.close());
    ASSERT_THROW(writer.close(), AssertionFailureException);
}


