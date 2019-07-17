//
// Created by xiamr on 7/17/19.
//

#include <gmock/gmock.h>
#include <boost/any.hpp>
#include "common.hpp"
#include "TypeUtility.hpp"

using namespace std;
using namespace testing;


TEST(TypeCheckTest, IsInt) {
    boost::any val = 1;
    ASSERT_THAT(TypeIs<int>()(val), Eq(true));

    val = 1.1;
    ASSERT_THAT(TypeIs<int>()(val), Eq(false));
}

TEST(TypeCheckTest, IsIntOrDouble) {
    boost::any val = 1;

    TypeIs<int, double> check{};

    ASSERT_THAT(check(val), Eq(true));

    val = 1.1;
    ASSERT_THAT(check(val), Eq(true));
}

TEST(TypePrettyNameTest, IntDoubleStringBoolAmberMask) {
    TypePrettyNames<int, double, string, bool, Atom::Node> check{};
    ASSERT_THAT(check(), ContainerEq(vector<string>{"int", "double", "string", "bool", "AmberMask"}));
}

TEST(getPrettyNameTest, BoostAny) {
    boost::any v = Atom::Node();
    ASSERT_THAT(getPrettyName(v), Eq("AmberMask"));
    v = shared_ptr<BasicAnalysis>();
    ASSERT_THAT(getPrettyName(v), Eq("BasicAnalysis"));
    v = shared_ptr<VectorSelector>();
    ASSERT_THAT(getPrettyName(v), Eq("VectorSelector"));
    v = Grid();
    ASSERT_THAT(getPrettyName(v), Eq("Grid"));
    v = 1;
    ASSERT_THAT(getPrettyName(v), Eq("int"));
    v = 1.1;
    ASSERT_THAT(getPrettyName(v), Eq("double"));
    v = true;
    ASSERT_THAT(getPrettyName(v), Eq("bool"));
}

TEST(AutoConvertTest, Int) {
    boost::any v = 1;
    int i = AutoConvert(v);
    ASSERT_THAT(i, Eq(1));
}


