//
// Created by xiamr on 5/22/19.
//


#define BOOST_TEST_MODULE "TrajectoryReader Gromacs TPR Topology Test"

#include <string>
#include <memory>

#include <boost/test/included/unit_test.hpp>

#include "trajectoryreader.hpp"
#include "frame.hpp"
#include "atom.hpp"

BOOST_AUTO_TEST_SUITE(test_trajectoryreader_tpr)

BOOST_AUTO_TEST_CASE(test_md_10ns_tpr) {

    class TestTrajectoryReader : public TrajectoryReader {

    public:
        std::shared_ptr<Frame> TestreadOneFrameTpr(){
            return readOneFrameTpr();
        }
    } reader;
    reader.add_topology("md-10ns.tpr");

    auto frame = reader.TestreadOneFrameTpr();

    BOOST_CHECK(frame->atom_list.size() == 7590);
    BOOST_CHECK(frame->molecule_list.size() == 2030 );

    auto atom1 = frame->atom_list.front();
    BOOST_TEST(atom1->charge == -8.10476e-01, boost::test_tools::tolerance(0.01));

    BOOST_TEST(atom1->atom_name == "O1");
    BOOST_TEST(atom1->type_name == "o");

    BOOST_TEST(atom1->residue_name.get() == "C27");

    BOOST_TEST(atom1->residue_num.get() == 1 );


}

BOOST_AUTO_TEST_SUITE_END()