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

TEST(TypePrettyNameTest, ALLAccept) {
    TypePrettyNames<int, double, string, bool, Atom::Node, Grid,
            shared_ptr<AbstractAnalysis>, shared_ptr<VectorSelector>> check{};
    ASSERT_THAT(check(), ContainerEq(vector<string>{
            "int", "double", "string", "bool", "AmberMask", "Grid", "BasicAnalysis", "VectorSelector"}));
}

TEST(getPrettyNameTest, BoostAny) {
    boost::any v = Atom::Node();
    ASSERT_THAT(getPrettyName(v), Eq("AmberMask"));
    v = shared_ptr<AbstractAnalysis>();
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

TEST(AutoConvertTest, DoubleFromInt) {
    boost::any v = 1;
    double i = AutoConvert(v);
    ASSERT_THAT(i, DoubleEq(1.0));
}

TEST(AutoConvertTest, DoubleFromDouble) {
    boost::any v = 1.0;
    double i = AutoConvert(v);
    ASSERT_THAT(i, DoubleEq(1.0));
}

TEST(AutoConvertTest, DoubleFromBool) {
    boost::any v = true;
    double i;
    ASSERT_NO_THROW((i = AutoConvert(v)));
    ASSERT_THAT(i, DoubleEq(1.0));
}

TEST(AutoConvertTest, Bool) {
    boost::any v = true;
    bool i = AutoConvert(v);
    ASSERT_THAT(i, Eq(true));
}

TEST(AutoConvertTest, String) {
    boost::any v = string("Hello");
    string i = AutoConvert(v);
    ASSERT_THAT(i, Eq("Hello"));
}

TEST(AutoConvertTest, AmberMask) {
    auto node = Atom::Node();
    boost::any v = node;
    Atom::Node i;
    ASSERT_NO_THROW((i = AutoConvert(v)));
}




