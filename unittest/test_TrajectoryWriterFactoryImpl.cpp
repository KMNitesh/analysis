//
// Created by xiamr on 6/28/19.
//

#include <gmock/gmock.h>
#include "utils/TrajectoryWriterFactoryImpl.hpp"
#include "utils/netcdf_writer.hpp"
#include "utils/trr_writer.hpp"
#include "utils/xtc_writer.hpp"
#include "utils/gro_writer.hpp"
#include "utils/ArcWriter.hpp"

using namespace std;
using namespace testing;


TEST(TrajectoryWriterFactoryImpl, make_instace) {
    shared_ptr<TrajectoryWriterFactoryInterface> factory = make_shared<TrajectoryWriterFactoryImpl>();
    ASSERT_TRUE(typeid(*factory->make_instance(FileType::NC)) == typeid(NetCDFWriter));
    ASSERT_TRUE(typeid(*factory->make_instance(FileType::TRR)) == typeid(TRRWriter));
    ASSERT_TRUE(typeid(*factory->make_instance(FileType::XTC)) == typeid(XTCWriter));
    ASSERT_TRUE(typeid(*factory->make_instance(FileType::GRO)) == typeid(GROWriter));
    ASSERT_TRUE(typeid(*factory->make_instance(FileType::ARC)) == typeid(ArcWriter));
}

TEST(TrajectoryWriterFactoryImpl, make_instace_Error) {
    shared_ptr<TrajectoryWriterFactoryInterface> factory = make_shared<TrajectoryWriterFactoryImpl>();
    ASSERT_THROW(factory->make_instance(FileType::MOL2), runtime_error);
    ASSERT_THROW(factory->make_instance(FileType::PRM), runtime_error);
    ASSERT_THROW(factory->make_instance(FileType::TPR), runtime_error);
    ASSERT_THROW(factory->make_instance(FileType::UnKnown), runtime_error);
}