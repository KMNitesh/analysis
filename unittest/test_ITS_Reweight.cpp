
#include <gmock/gmock.h>
#include <boost/spirit/include/qi.hpp>
#include "others/ITS_Reweight.hpp"

using namespace testing;

TEST(Test_ITS_Reweight, read_fb) {

    std::stringstream ss(
            "0 : 0.00247875,0.00247875,0.00247875,0.00247875,\n"
            "100 : -0,0.00587086,0.0113129,0.0163715,\n");

    auto[ok, fb] = ITS_Reweight::read_fb(ss);

    ASSERT_TRUE(ok);
    ASSERT_THAT(fb.size(), Eq(2));

    ASSERT_THAT(boost::fusion::at_c<0>(fb[0]), Eq(0));
    ASSERT_THAT(boost::fusion::at_c<1>(fb[0]), Pointwise(DoubleEq(), {0.00247875, 0.00247875, 0.00247875, 0.00247875}));

    ASSERT_THAT(boost::fusion::at_c<0>(fb[1]), Eq(100));
    ASSERT_THAT(boost::fusion::at_c<1>(fb[1]), Pointwise(DoubleEq(), {-0.0, 0.00587086, 0.0113129, 0.0163715}));
}

TEST(Test_ITS_Reweight, read_pot) {

    std::stringstream ss(
            "#          time          totpot\n"
            "              0       -40311.64\n"
            "            0.5     -38957.7524\n");

    auto[ok, pot] = ITS_Reweight::read_pot(ss);

    ASSERT_TRUE(ok);
    ASSERT_THAT(pot[0], Pair(DoubleEq(0.0), DoubleEq(-40311.64)));
    ASSERT_THAT(pot[1], Pair(DoubleEq(0.5), DoubleEq(-38957.7524)));
}

TEST(Test_ITS_Reweight, eol) {
    std::stringstream ss(
            "0 : 0.00247875,0.00247875,0.00247875,0.00247875,\n100 : -0,0.00587086,0.0113129,0.0163715,\n");
    ss.unsetf(std::ios::skipws);

    boost::spirit::istream_iterator begin(ss);
    boost::spirit::istream_iterator end;

    using namespace boost::spirit::qi;
    using namespace boost::phoenix;

    auto parser = copy(+(omit[int_] >> ':' >> +(double_ >> ',') >> eol));

    ASSERT_TRUE(phrase_parse(begin, end, parser, ascii::space - eol));
}





