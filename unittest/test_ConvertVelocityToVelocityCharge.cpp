//
// Created by xiamr on 8/17/19.
//

#include <gmock/gmock.h>
#include "ConvertVelocityToVelocityCharge.hpp"
#include "common.hpp"
#include "atom.hpp"
#include "frame.hpp"

using namespace testing;

class TRRWriterMock : public TRRWriter {
public:
    MOCK_METHOD1(open, void(
            const std::string &filename));

    MOCK_METHOD1(write, void(
            const std::shared_ptr<Frame> &frame));

    MOCK_METHOD0(close, void());

    MOCK_METHOD1(setWriteVelocities, void (
            bool writeVelocities));

    MOCK_METHOD1(setCurrentTime, void (
            gmx::real t));
};

TEST(ConvertVelocityToVelocityCharge, Default) {

    auto mock = std::make_unique<TRRWriterMock>();
    auto mock_ptr = mock.get();

    class ConvertVelocityToVelocityChargeTest : public ConvertVelocityToVelocityCharge {
    public:
        explicit ConvertVelocityToVelocityChargeTest(std::unique_ptr<TRRWriter> writer)
                : ConvertVelocityToVelocityCharge(std::move(writer)) {
            trr_vq_outfilename = "test_out.trr";
        }

        void setSelectedMols(const std::unordered_set<std::shared_ptr<Molecule>> &mols) {
            selected_mols = mols;
        }

    protected:
        void do_select_mol(std::shared_ptr<Frame> &frame) override {

        }
    } convertor(std::move(mock));

    ASSERT_THAT(enable_read_velocity, Eq(true));

    auto frame = std::make_shared<Frame>();
    frame->setCurrentTime(100.0);

    auto atom1 = std::make_shared<Atom>();
    atom1->x = 1.0;
    atom1->y = 2.0;
    atom1->z = 2.0;
    atom1->vx = 1.0;
    atom1->vy = 2.0;
    atom1->vz = -1.0;
    atom1->charge = 1.0;

    frame->atom_list.push_back(atom1);

    auto atom2 = std::make_shared<Atom>();
    atom2->x = 1.0;
    atom2->y = 3.0;
    atom2->z = 2.0;
    atom2->vx = -1.0;
    atom2->vy = -2.0;
    atom2->vz = -1.0;
    atom2->charge = 1.5;

    frame->atom_list.push_back(atom2);

    auto mol = std::make_shared<Molecule>();
    atom1->molecule = mol;
    atom2->molecule = mol;

    mol->atom_list.push_back(atom1);
    mol->atom_list.push_back(atom2);

    frame->molecule_list.push_back(mol);

    convertor.setSelectedMols({mol});

    EXPECT_CALL(*mock_ptr, open(StrEq("test_out.trr"))).Times(1);
    EXPECT_CALL(*mock_ptr, setWriteVelocities(true)).Times(1);
    EXPECT_CALL(*mock_ptr, setCurrentTime(FloatEq(100.0))).Times(1);
    EXPECT_CALL(*mock_ptr, write(_)).Times(1);
    EXPECT_CALL(*mock_ptr, close()).Times(1);


    convertor.processFirstFrame(frame);
    convertor.process(frame);
    convertor.print(std::cout);


    ASSERT_THAT(frame->atom_list.front()->x, DoubleEq(1.0));
    ASSERT_THAT(frame->atom_list.front()->y, DoubleEq(2.0));
    ASSERT_THAT(frame->atom_list.front()->z, DoubleEq(2.0));

    ASSERT_THAT(frame->atom_list.front()->vx, DoubleEq(0.0));
    ASSERT_THAT(frame->atom_list.front()->vy, DoubleEq(0.0));
    ASSERT_THAT(frame->atom_list.front()->vz, DoubleEq(0.0));

    ASSERT_THAT(frame->atom_list.back()->x, DoubleEq(1.0));
    ASSERT_THAT(frame->atom_list.back()->y, DoubleEq(3.0));
    ASSERT_THAT(frame->atom_list.back()->z, DoubleEq(2.0));

    ASSERT_THAT(frame->atom_list.back()->vx, DoubleEq(-0.5));
    ASSERT_THAT(frame->atom_list.back()->vy, DoubleEq(-1.0));
    ASSERT_THAT(frame->atom_list.back()->vz, DoubleEq(-2.5));

}



