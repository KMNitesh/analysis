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
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/phoenix.hpp>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

#include <boost/format.hpp>

class Atom;

class Frame;

class Forcefield;

constexpr int ATOM_MAX = 10000;

constexpr double radian = 57.29577951308232088;

// global variables

extern bool enable_read_velocity;
extern bool enable_tbb;
extern bool enable_outfile;

extern Forcefield forcefield;
extern bool enable_forcefield;
extern std::fstream outfile;

enum class FileType {
    XTC,
    TRR,
    NC,
    ARC,
    TPR,
    MOL2,
    PRM,
    UnKnown
};

FileType getFileType(const std::string &filename);

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


template<typename T, typename = std::enable_if_t<std::is_same_v<T, int> or std::is_same_v<T, double>>>
T choose(T min, T max, const std::string &prompt, bool hasdefault = false, T value = T()) {
    while (true) {
        std::string input_line = input(prompt);
        boost::trim(input_line);
        if (input_line.empty()) {
            if (!hasdefault) continue;
            return value;
        }
        try {
            auto option = boost::lexical_cast<T>(input_line);
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

template<typename T>
struct make_shared_f {
    template<typename... A>
    struct result {
        typedef std::shared_ptr<T> type;
    };

    template<typename... A>
    typename result<A...>::type operator()(A &&... a) const {
        return std::make_shared<T>(std::forward<A>(a)...);
    }
};


template<typename T, typename... _Args>
inline auto make_shared_(_Args &&... __args) {
    return boost::phoenix::function<make_shared_f<T>>()(std::forward<_Args>(__args)...);
}


po::options_description make_program_options();


class range_object {
private:
    int _end;
    int _step;
    int _curr;

public:

    range_object(int start, int end, int step = 1) : _end(end), _step(step), _curr(start) {};

    const range_object &begin() const { return *this; }

    const range_object &end() const { return *this; }

    bool operator!=(const range_object &) const {
        return _step > 0 ? _curr < _end : _curr > _end;
    }

    void operator++() {
        _curr += _step;
    }

    int operator*() const {
        return _curr;
    }
};


/*
 *  Python-like range function in C++
 */
inline range_object range(int start, int end, int step = 1) {
    if (step == 0) {
        throw std::runtime_error("zero increment step size!!");
    }
    if ((start > end && step > 0) || (start < end && step < 0)) {
        throw std::runtime_error("wrong direction increment step size!!");
    }
    return range_object(start, end, step);
}

inline range_object range(int end) {
    if (end < 0) {
        throw std::runtime_error("wrong end value !!");
    }
    return range_object(0, end, 1);
}


template<typename Iterable>
class enumerate_object {
private:
    Iterable _iter;
    int _size;
    int _step;
    decltype(std::begin(_iter)) _begin;
    const decltype(std::end(_iter)) _end;

public:
    enumerate_object(Iterable iter, int start = 0, int step = 1) :
            _iter(iter),
            _size(start),
            _step(step),
            _begin(std::begin(iter)),
            _end(std::end(iter)) {}

    const enumerate_object &begin() const { return *this; }

    const enumerate_object &end() const { return *this; }

    bool operator!=(const enumerate_object &) const {
        return _begin != _end;
    }

    void operator++() {
        ++_begin;
        _size += _step;
    }

    auto operator*() const
    -> std::pair<std::size_t, decltype(*_begin)> {
        return {_size, *_begin};
    }

    enumerate_object &start(const int &start) {
        _size = start;
        return *this;
    }

    enumerate_object &step(const int &step) {
        if (step == 0) {
            throw std::runtime_error("Zero Increment Not Allowed");
        }
        _step = step;
        return *this;
    }
};

/*
 *  Python-like enumerate function in C++
 */

template<typename Iterable>
auto enumerate(Iterable &&iter)
-> enumerate_object<Iterable> {
    return {std::forward<Iterable>(iter)};
}

template<typename T, typename... Args>
auto format(T &&s, Args &&... args) {
    return (boost::format(std::forward<T>(s)) %  ... % (std::forward<Args>(args)));
}

std::string print_cmdline(int argc, const char *const argv[]);

#endif //TINKER_COMMON_HPP
