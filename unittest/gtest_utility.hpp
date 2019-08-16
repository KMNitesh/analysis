//
// Created by xiamr on 6/28/19.
//

#ifndef TINKER_GTEST_UTILITY_HPP
#define TINKER_GTEST_UTILITY_HPP

#include <gmock/gmock.h>

using namespace testing;


template<typename T, typename T1, typename T2>
bool array_eq_impl(const T1 &a1, const T2 a2, int n) {
    auto p1 = (T *) a1;
    auto p2 = (T *) a2;
    while (n > 0) {
        if (abs(*(p1++) - *(p2++)) > std::numeric_limits<T>::epsilon()) {
            return false;
        }
        n--;
    }
    return true;
}

MATCHER_P2(FLOAT_ARRAY_EQ, s2, n, "float array compare") {
    assert(n > 0);
    static_assert(std::is_same_v<std::decay_t<arg_type>, std::decay_t<s2_type> >);
    return array_eq_impl<float>(arg, s2, n);
}

MATCHER_P2(DOUBLE_ARRAY_EQ, s2, n, "double array compare") {
    assert(n > 0);
    static_assert(std::is_same_v<std::decay_t<arg_type>, std::decay_t<s2_type> >);
    return array_eq_impl<double>(arg, s2, n);
}


#endif //TINKER_GTEST_UTILITY_HPP
