//
// Created by xiamr on 10/14/19.
//

#include "utils/std.hpp"
#include <gmock/gmock.h>
#include "ana_module/RMSFCal.hpp"
#include "data_structure/frame.hpp"

using namespace testing;

class TestRMSF : public RMSFCal, public Test {
protected:
    void SetUp() override {

        atom1 = std::make_shared<Atom>();
        atom1->seq = 1;
        atom1->residue_name = "MOL";
        atom1->residue_num = 1;

        atom2 = std::make_shared<Atom>();
        atom2->seq = 2;
        atom2->residue_name = "MOL";
        atom2->residue_num = 1;

        atom3 = std::make_shared<Atom>();
        atom3->seq = 3;
        atom3->residue_name = "MOL";
        atom3->residue_num = 1;

        atom4 = std::make_shared<Atom>();
        atom4->seq = 4;
        atom4->residue_name = "TOL";
        atom4->residue_num = 2;

        atoms_for_superpose.push_back(atom1);
        atoms_for_superpose.push_back(atom2);
        atoms_for_superpose.push_back(atom3);

        atoms_for_rmsfcalc.push_back(atom4);

        frame = std::make_shared<Frame>();

        frame->atom_list.push_back(atom1);
        frame->atom_list.push_back(atom2);
        frame->atom_list.push_back(atom3);
        frame->atom_list.push_back(atom4);

        allocate_array_memory();

        std::tie(atom1->x, atom1->y, atom1->z) = std::make_tuple(0.0, 0.0, 0.0);
        std::tie(atom2->x, atom2->y, atom2->z) = std::make_tuple(1.0, 0.0, 0.0);
        std::tie(atom3->x, atom3->y, atom3->z) = std::make_tuple(0.0, 1.0, 0.0);
    }

    std::shared_ptr<Atom> atom1, atom2, atom3, atom4;
    std::shared_ptr<Frame> frame;
};

TEST_F(TestRMSF, IntegrationTest) {

    std::vector<std::tuple<double, double, double>> coordinates{
            {10, 10, 10},
            {11, 10, 10},
            {9,  10, 10},
            {8,  10, 10},
            {7,  10, 10}
    };

    for (const auto &element : coordinates) {
        std::tie(atom4->x, atom4->y, atom4->z) = element;
        process(frame);
    }

    calculate_average_structure();

    ASSERT_THAT(rmsvalue(0), DoubleEq(std::sqrt(2)));
}