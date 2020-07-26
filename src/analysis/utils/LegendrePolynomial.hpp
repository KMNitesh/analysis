//
// Created by xiamr on 9/7/19.
//

#ifndef TINKER_LEGENDREPOLYNOMIAL_HPP
#define TINKER_LEGENDREPOLYNOMIAL_HPP

#include "std.hpp"

template <int I> struct LegendrePolynomial;

template <> struct LegendrePolynomial<1> {
    template <typename T> auto operator()(T x) { return x; }
};

template <> struct LegendrePolynomial<2> {
    template <typename T> auto operator()(T x) { return 0.5 * (3 * x * x - 1); }
};

template <> struct LegendrePolynomial<3> {
    template <typename T> auto operator()(T x) { return 0.5 * (5 * x * x * x - 3 * x); }
};

template <> struct LegendrePolynomial<4> {
    template <typename T> auto operator()(T x) { return 1.0 / 8.0 * (35 * x * x * x * x - 30 * x * x + 3); }
};

using LegendrePolynomialLevel1 = LegendrePolynomial<1>;
using LegendrePolynomialLevel2 = LegendrePolynomial<2>;
using LegendrePolynomialLevel3 = LegendrePolynomial<3>;
using LegendrePolynomialLevel4 = LegendrePolynomial<4>;

inline const std::unordered_map<int, std::string_view> LegendreStr{
    {1, "P1 = x"},
    {2, "P2 = (1/2)(3x^2 - 1)"},
    {3, "P3 = (1/2)(5x^3 - 3x)"},
    {4, "P4 = (1/8)(35x^4 - 30x^2 + 3)"},
};

#endif // TINKER_LEGENDREPOLYNOMIAL_HPP
