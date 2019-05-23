//
// Created by xiamr on 3/17/19.
//

#ifndef TINKER_COMMON_HPP
#define TINKER_COMMON_HPP

#include "config.h"
#include <string>
#include <vector>
#include <type_traits>
#include <iostream>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

class Atom;
class Frame;
class Forcefield;



// global variables

extern bool enable_read_velocity;
extern bool enable_tbb;
extern bool enable_outfile;

extern Forcefield forcefield;
extern bool enable_forcefield;

bool file_exist(const std::string &name);

std::vector<std::string> split(const std::string &str, const std::string &sep);
std::vector<std::string> split(const std::string &str);

std::string input(const std::string &prompt = "");

template<typename T>
struct type_name_string;

template<>
struct type_name_string<int> {
    constexpr static auto value = "int";
};

template<>
struct type_name_string<double> {
    constexpr static auto value = "double";
};


template<typename T, typename = std::enable_if_t<std::is_same<T, int>::value or std::is_same<T, double>::value>>
T choose(T min, T max, const std::string &prompt, bool hasdefault = false, T value = T()) {
    while (true) {
        std::string input_line = input(prompt);
        boost::trim(input_line);
        if (input_line.empty()) {
            if (!hasdefault) continue;
            return value;
        }
        try {
            int option = boost::lexical_cast<T>(input_line);
            if (option >= min and option <= max) return option;

            std::cerr << "must be a " << type_name_string<T>::value << " range " << min << " and " << max
                 << "! please retype!\n";
        } catch (boost::bad_lexical_cast &e) {
            std::cerr << "must be a " << type_name_string<T>::value << " ! please retype!" << e.what() << std::endl;
        }
    }
}

std::string ext_filename(const std::string &filename);

std::string choose_file(const std::string &prompt, bool exist, std::string ext = "", bool can_empty = false);

template<typename T>
T sign(T &x, T &y) { return y > 0 ? std::abs(x) : -std::abs(x); }

double
atom_distance(const std::shared_ptr<Atom> &atom1, const std::shared_ptr<Atom> &atom2, std::shared_ptr<Frame> &frame);

double
atom_distance2(const std::shared_ptr<Atom> &atom1, const std::shared_ptr<Atom> &atom2, std::shared_ptr<Frame> &frame);
#endif //TINKER_COMMON_HPP
