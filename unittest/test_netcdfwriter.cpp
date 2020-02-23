//
// Created by xiamr on 6/28/19.
//

#include <gmock/gmock.h>
#include "utils/AmberNetcdf.h"
#include "utils/NetcdfInterface.hpp"
#include "utils/netcdf_writer.hpp"
#include "data_structure/atom.hpp"
#include "data_structure/frame.hpp"
#include "gtest_utility.hpp"

using namespace std;
using namespace testing;

class NetcdfMock : public NetcdfInterface {
public:
    MOCK_METHOD1(netcdfClose, int(
            struct AmberNetcdf *A));

    MOCK_METHOD4(netcdfCreate, int(
            struct AmberNetcdf *A,
            const char *filename,
            int natom,
            int isBox));

    MOCK_METHOD3(netcdfWriteNextFrame, int(
            struct AmberNetcdf *A,
            double *X,
            double *box));
};


class NetCDFWriterTest : public NetCDFWriter {
protected:
    NetcdfInterface *getNetcdfImpl() override {
        return mock;
    }

public:
    NetcdfMock *mock;

    NetCDFWriterTest(NetcdfMock *mock) : mock(mock) {}
};

TEST(NetCDFWriter, OpenAndCloseFile) {

    NetcdfMock mock;
    NetCDFWriterTest writer(&mock);

    std::string filename = "output.nc";
    writer.open(filename);
    writer.close();
}

MATCHER_P(AmberNetcdfEq, ref, "struct AmberNetcdf compare") {
    return arg->ncid == ref.ncid;
}

TEST(NetCDFWriter, WriteFrameToFile) {

    NetcdfMock mock;
    NetCDFWriterTest writer(&mock);

    std::string filename = "output.nc";
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
    atom1->seq = 2;

    frame->atom_list = {atom1, atom2};
    frame->atom_map = {{atom1->seq, atom1},
                       {atom2->seq, atom2}};

    double box[] = {10, 10, 10, 90, 90, 90};
    double x[] = {-9.9, 8.89, 2.1, 0.90, 3.19, 4.10};

    // int netcdfCreate(struct AmberNetcdf *A, const char *filename, int natom, int isBox)
    struct AmberNetcdf NC{};
    NC.ncatom = 2;
    NC.ncatom3 = 6;
    EXPECT_CALL(mock, netcdfCreate(_, StrEq("output.nc"), 2, 1)).WillOnce(DoAll(SetArgPointee<0>(NC), Return(0)));
    // int netcdfWriteNextFrame(struct AmberNetcdf *A, double *X, double *box) = 0;
    EXPECT_CALL(mock, netcdfWriteNextFrame(AmberNetcdfEq(NC), DOUBLE_ARRAY_EQ(x, 6), DOUBLE_ARRAY_EQ(box, 6)))
            .WillOnce(Return(0))
            .WillOnce(Return(0));
    EXPECT_CALL(mock, netcdfClose(AmberNetcdfEq(NC))).WillOnce(Return(0));

    writer.open(filename);
    writer.write(frame);
    writer.write(frame);
    writer.close();
}

TEST(NetCDFWriter, ErrOpenException) {

    NetcdfMock mock;
    NetCDFWriterTest writer(&mock);
    EXPECT_CALL(mock, netcdfCreate(_, StrEq("output.nc"), 0, 1)).WillOnce(Return(1));
    writer.open("output.nc");
    auto frame = make_shared<Frame>();
    ASSERT_THROW(writer.write(frame), std::runtime_error);

}
