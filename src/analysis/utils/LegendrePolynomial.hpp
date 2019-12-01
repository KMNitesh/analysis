//
// Created by xiamr on 9/7/19.
//

#ifndef TINKER_LEGENDREPOLYNOMIAL_HPP
#define TINKER_LEGENDREPOLYNOMIAL_HPP

#include "std.hpp"

class LegendrePolynomialLevel1 {
public:
    template <typename T>
    T operator()(T x){
        return x;
    }
};

class LegendrePolynomialLevel2 {
public:
    template <typename T>
    T operator()(T x){
        return 0.5 * (3 * x * x - 1); 
    }
};

class LegendrePolynomialLevel3 {
public:
    template <typename T>
    T operator()(T x){
       return 0.5 * (5 * x * x * x - 3 * x); 
    }
};

class LegendrePolynomialLevel4 {
public:
    template <typename T>
    T operator()(T x){
       return 1.0 / 8.0 * (35 * x * x * x * x - 30 * x * x + 3);
    }
};

extern const std::unordered_map<int, std::string> LegendreStr;

#endif //TINKER_LEGENDREPOLYNOMIAL_HPP
