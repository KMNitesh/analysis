
#include <gmock/gmock.h>
#include <boost/optional.hpp>
#include <boost/fusion/include/at_c.hpp>
#include "others/ADCHCharge.hpp"

using namespace testing;

TEST(TestADCHCharge, read_charge) {
    auto charge = ADCHCharge::read_charge(" Atom:    1Am  Corrected charge:    1.164452  Before:    1.124097");
    ASSERT_THAT(charge.has_value(), Eq(true));
    ASSERT_THAT(boost::fusion::at_c<0>(*charge), Eq("Am"));
    ASSERT_THAT(boost::fusion::at_c<1>(*charge), DoubleEq(1.164452));
}