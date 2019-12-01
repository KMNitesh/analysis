//
// Created by xiamr on 6/15/19.
//

#include <gmock/gmock.h>
#include "utils/common.hpp"

using namespace testing;

/*
 *      XTC,TRR,NC,ARC,TPR,MOL2,PRM,UnKnown
 */

TEST(getFileType, XTCFile) {
    ASSERT_EQ(getFileType("abx.xtc"), FileType::XTC);
}

TEST(getFileType, TRRFile) {
    ASSERT_EQ(getFileType("abx.trr"), FileType::TRR);
}

TEST(getFileType, NCFile) {
    ASSERT_EQ(getFileType("abx.nc"), FileType::NC);
}

TEST(getFileType, MDCRDFile) {
    ASSERT_EQ(getFileType("abx.mdcrd"), FileType::NC);
}

TEST(getFileType, TPRFile) {
    ASSERT_EQ(getFileType("abx.tpr"), FileType::TPR);
}

TEST(getFileType, MOL2File) {
    ASSERT_EQ(getFileType("abx.mol2"), FileType::MOL2);
}

TEST(getFileType, PRMFile) {
    ASSERT_EQ(getFileType("abx.prm"), FileType::PRM);
}

TEST(getFileType, GROFile) {
    ASSERT_EQ(getFileType("abx.gro"), FileType::GRO);
}

TEST(getFileType, TrajFile) {
    ASSERT_EQ(getFileType("abx.traj"), FileType::TRAJ);
}

TEST(getFileType, CapitalNameOfTrajFile) {
    ASSERT_EQ(getFileType("abx.TRAJ"), FileType::TRAJ);
}

TEST(getFileType, NoExtension) {
    ASSERT_EQ(getFileType("abx"), FileType::UnKnown);
}

TEST(getFileType, NoExtensionWithDot) {
    ASSERT_EQ(getFileType("abx."), FileType::UnKnown);
}

TEST(getFileType, WithOneDot) {
    ASSERT_EQ(getFileType("."), FileType::UnKnown);
}

TEST(getFileType, WithDoubleDot) {
    ASSERT_EQ(getFileType(".."), FileType::UnKnown);
}

TEST(getFileType, NoBaseName) {
    ASSERT_EQ(getFileType(".xtc"), FileType::XTC);
}

TEST(getFileType, WrongExtension) {
    ASSERT_EQ(getFileType(".xt"), FileType::UnKnown);
}

