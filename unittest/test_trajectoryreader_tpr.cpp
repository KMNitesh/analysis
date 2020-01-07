
#include <gmock/gmock.h>
#include "trajectory_reader/trajectoryreader.hpp"
#include "data_structure/frame.hpp"
#include "data_structure/atom.hpp"

using namespace testing;

TEST(test_trajectoryreader_tpr, test_md_10ns_tpr) {

    TrajectoryReader reader;
    reader.add_topology("tpr_test_system1.tpr");

    auto frame = reader.readTopology();

    ASSERT_THAT(frame->atom_list.size(), Eq(7590));
    ASSERT_THAT(frame->molecule_list.size(), Eq(2030));

    auto atom1 = frame->atom_list.front();
    ASSERT_THAT(atom1->charge.get(), DoubleNear(-8.10476e-01, 0.01));

    ASSERT_THAT(atom1->atom_name, Eq("O1"));
    ASSERT_THAT(atom1->type_name, Eq("o"));

    ASSERT_THAT(atom1->residue_name.get(), Eq("C27"));

    ASSERT_THAT(atom1->residue_num.get(), Eq(1));
}

