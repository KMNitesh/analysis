//
// Created by xiamr on 6/28/19.
//

#include <gmock/gmock.h>
#include "TrajectoryWriterFactoryImpl.hpp"
#include "netcdf_writer.hpp"
#include "trr_writer.hpp"
#include "xtc_writer.hpp"
#include "gro_writer.hpp"

using namespace std;
using namespace testing;


TEST(TrajectoryWriterFactoryImpl, make_instace) {
    shared_ptr<TrajectoryWriterFactoryInterface> factory = make_shared<TrajectoryWriterFactoryImpl>();
    ASSERT_TRUE(typeid(*factory->make_instance(FileType::NC)) == typeid(NetCDFWriter));
    ASSERT_TRUE(typeid(*factory->make_instance(FileType::TRR)) == typeid(TRRWriter));
    ASSERT_TRUE(typeid(*factory->make_instance(FileType::XTC)) == typeid(XTCWriter));
    ASSERT_TRUE(typeid(*factory->make_instance(FileType::GRO)) == typeid(GROWriter));
}

TEST(TrajectoryWriterFactoryImpl, make_instace_Error) {
    shared_ptr<TrajectoryWriterFactoryInterface> factory = make_shared<TrajectoryWriterFactoryImpl>();
    ASSERT_THROW(factory->make_instance(FileType::ARC), runtime_error);
    ASSERT_THROW(factory->make_instance(FileType::MOL2), runtime_error);
    ASSERT_THROW(factory->make_instance(FileType::PRM), runtime_error);
    ASSERT_THROW(factory->make_instance(FileType::TPR), runtime_error);
    ASSERT_THROW(factory->make_instance(FileType::UnKnown), runtime_error);
}