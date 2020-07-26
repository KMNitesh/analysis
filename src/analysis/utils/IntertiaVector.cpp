
#include "IntertiaVector.hpp"
#include "ana_module/Angle.hpp"
#include <boost/algorithm/string.hpp>

IntertiaVector::IntertiaVector(const AmberMask &mask, const std::string &axis) : mask(mask) {

    if (auto lowcase_axis = boost::to_lower_copy(axis); lowcase_axis == "max") {
        this->axis = Angle::AxisType::MAX;
    } else if (lowcase_axis == "min") {
        this->axis = Angle::AxisType::MIN;
    } else {
        throw std::runtime_error("unknown axis type : " + axis);
    }
}