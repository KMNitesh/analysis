//
// Created by xiamr on 3/17/19.
//

#ifndef TINKER_COMMON_HPP
#define TINKER_COMMON_HPP

#include "std.hpp"
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/phoenix.hpp>
#include <boost/optional.hpp>
#include <boost/program_options.hpp>
#include <boost/type_index.hpp>
#include <boost/math/constants/constants.hpp>

namespace po = boost::program_options;

#include <boost/format.hpp>

class Atom;

class Frame;

class Forcefield;

constexpr int ATOM_MAX = 10000;

constexpr double radian = 57.29577951308232088;
constexpr double pi = boost::math::constants::pi<double>();

constexpr double avogadro_constant = 6.022140857e23;
constexpr double kb = 1.380649e-23; // unit: J/K

// global variables

extern bool enable_read_velocity;
extern bool enable_tbb;
extern bool enable_outfile;

extern Forcefield forcefield;
extern bool enable_forcefield;


enum class FileType {
    XTC,
    TRR,
    NC,
    ARC,
    TPR,
    MOL2,
    PRM,
    GRO,
    TRAJ,
    UnKnown
};

FileType getFileType(const std::string &filename);

std::vector<std::string> split(const std::string &str, const std::string &sep);

std::vector<std::string> split(const std::string &str);

std::vector<std::string> split_quoted(const std::string &str);

std::string input(const std::string &prompt = "", std::istream &in = std::cin, std::ostream &out = std::cout);


template<typename T, typename = std::enable_if_t<std::is_same_v<T, int> or std::is_same_v<T, double>>>
T choose(T min, T max, const std::string &prompt, bool hasdefault = false, T value = T(),
         std::istream &in = std::cin, std::ostream &out = std::cout) {
    while (true) {
        std::string input_line = input(prompt, in, out);
        boost::trim(input_line);
        if (input_line.empty()) {
            if (!hasdefault) continue;
            return value;
        }
        try {
            auto option = boost::lexical_cast<T>(input_line);
            if (option >= min and option <= max) return option;

            out << "must be a " << boost::typeindex::type_id<T>().pretty_name() << " range " << min << " and " << max
                << "! please retype!\n";
        } catch (boost::bad_lexical_cast &e) {
            out << "must be a " << boost::typeindex::type_id<T>().pretty_name() << " ! please retype!" << e.what()
                << '\n';
        }
    }
}

template<typename T>
class Default {
public:
    explicit Default(T x) : x(x) {}

    T getValue() { return x; }

private:
    T x;
};

bool choose_bool(const std::string &prompt, boost::optional<Default<bool>> defaultValue = {},
                 std::istream &in = std::cin, std::ostream &out = std::cout);

std::string ext_filename(const std::string &filename);

std::string choose_file(const std::string &prompt, bool exist, std::string ext = "", bool can_empty = false,
                        std::istream &in = std::cin, std::ostream &out = std::cout);

template<typename T>
T sign(const T &x, const T &y) { return y > 0 ? std::abs(x) : -std::abs(x); }

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
class enumerate_object;

template<typename Iterable>
class enumerate_iterator {
    decltype(std::begin(Iterable())) it;
    int _start;
    int _step;
public:
    enumerate_iterator(decltype(std::begin(Iterable())) it, int start, int step) : it(it), _start(start), _step(step) {}

    void operator++() {
        _start += _step;
        ++it;
    }


    bool operator!=(const enumerate_iterator<Iterable> &other) {
        return this->it != other.it;
    }

    auto operator*() const {
        return std::make_pair(_start, *it);
    }

};

template<typename Iterable>
class enumerate_object {
    friend class enumerate_iterator<Iterable>;

    Iterable &_iter;
    int _start;
    int _step;

public:
    explicit enumerate_object(Iterable &iter, int start = 0, int step = 1) :
            _iter(iter),
            _start(start),
            _step(step) {
//    std::cout << "Constructor1\n";
    }

    explicit enumerate_object(Iterable &&iter, int start = 0, int step = 1) :
            _iter(iter),
            _start(start),
            _step(step) {
//    std::cout << "Constructor2\n";
    }

    virtual ~enumerate_object() {
//    std::cout << "Destructor\n";
    }

    auto begin() { return enumerate_iterator<Iterable>(std::begin(_iter), _start, _step); }

    auto end() { return enumerate_iterator<Iterable>(std::end(_iter), _start, _step); }

    enumerate_object &start(int start) {
        _start = start;
        return *this;
    }

    enumerate_object &step(int step) {
        if (step == 0) {
            throw std::runtime_error("ERROR !! step cannot zero!");
        }
        _step = step;
        return *this;
    }
};

/*
 *    Python-like enumerate function in C++
 *
 */

template<typename Iterable>
auto enumerate(Iterable &&iter) {
    return enumerate_object<std::remove_reference_t<Iterable>>(std::forward<Iterable>(iter));
}

template<typename T>
auto enumerate(std::initializer_list<T> &&iter) {
    return enumerate_object<std::initializer_list<T>>(std::forward<std::initializer_list<T>>(iter));
}

template<typename T, typename... Args>
auto format(T &&s, Args &&... args) {
    return (boost::format(std::forward<T>(s)) %  ... % (std::forward<Args>(args)));
}

std::string print_cmdline(int argc, const char *const argv[]);


template<typename Iterable>
class PushIterable_object {
public:
    using value_type = typename Iterable::value_type;
private:
    Iterable &_iter;

    decltype(std::begin(_iter)) _begin;
    const decltype(std::end(_iter)) _end;

    std::stack<value_type> queue;

public:
    explicit PushIterable_object(Iterable &iter) :
            _iter(iter),
            _begin(std::begin(iter)),
            _end(std::end(iter)) {}

    const PushIterable_object &begin() const { return *this; }

    const PushIterable_object &end() const { return *this; }


    bool operator!=(const PushIterable_object &) const {
        return _begin != _end || !queue.empty();
    }

    bool operator==(const PushIterable_object &) const {
        return _begin == _end && queue.empty();
    }

    void operator++() {
        if (!queue.empty()) {
            queue.pop();
        } else {
            ++_begin;
        }

    }

    auto operator*() const {
        if (!queue.empty()) {
            return queue.top();
        } else {
            assert(_begin != _end);
            return *_begin;
        }
    }

    void push_back(value_type value) {
        queue.push(value);
    }

    bool empty() const {
        return _begin == _end && queue.empty();
    }

    auto next() {
        if (!queue.empty()) {
            auto value = queue.top();
            queue.pop();
            return value;
        } else {
            assert(_begin != _end);
            auto value = *_begin;
            _begin++;
            return value;
        }
    }
};

template<typename Iterable, typename = std::enable_if_t<
        std::is_integral_v<typename Iterable::value_type> && !std::is_same_v<typename Iterable::value_type, bool> >>
class CombineSeq {
public:
    using value_type = typename Iterable::value_type;
private:
    PushIterable_object<Iterable> _iter;

    boost::optional<std::string> _curr;

public:
    explicit CombineSeq(Iterable &iter) :
            _iter(iter) {}

    CombineSeq &begin() { return *this; }

    CombineSeq &end() { return *this; }

    CombineSeq(const CombineSeq &other) :
            _iter(other._iter), _curr(other._curr) {}


    bool operator!=(const CombineSeq &) const {
        return !_iter.empty() || _curr.has_value();
    }

    void operator++() {
        _curr = {};
        if (_iter.empty()) {
            return;
        }
        auto first = _iter.next();
        _curr = boost::lexical_cast<std::string>(first);
        if (_iter.empty()) {
            return;
        }

        auto second = _iter.next();
        if (_iter.empty() || second != first + 1) {
            _iter.push_back(second);
            return;
        }

        auto third = _iter.next();
        if (third != second + 1) {
            _iter.push_back(third);
            _iter.push_back(second);
            return;
        }
        second = third;
        for (; !_iter.empty();) {
            auto value = _iter.next();
            if (value != second + 1) {
                _iter.push_back(value);
                break;
            }
            second = value;
        }
        _curr = _curr.value() + '-' + boost::lexical_cast<std::string>(second);
    }

    auto operator*() {
        if (!_curr.has_value()) {
            this->operator++();
        }
        return _curr.value();
    }
};

template<typename Iterable>
auto PushIterable(Iterable &&iter) {
    return PushIterable_object(iter);
}


template<typename T>
auto PushIterable(std::initializer_list<T> &&iter) {
    return PushIterable_object(iter);
}


template<typename Iterable>
auto combine_seq(Iterable &&iter) {
    return CombineSeq(iter);
}


template<typename T>
auto combine_seq(std::initializer_list<T> &&iter) {
    return CombineSeq(iter);
}


template<typename T>
T dot_multiplication(const std::tuple<T, T, T> &lhs, const std::tuple<T, T, T> &rhs) {
    return std::get<0>(lhs) * std::get<0>(rhs) +
           std::get<1>(lhs) * std::get<1>(rhs) +
           std::get<2>(lhs) * std::get<2>(rhs);
}

template<typename T>
T dot_multiplication(const T &lhs, const T &rhs) {
    return lhs * rhs;
}


template<typename T>
std::tuple<T, T, T> cross_multiplication(const std::tuple<T, T, T> &lhs, const std::tuple<T, T, T> &rhs) {
    auto[u1, u2, u3] = lhs;
    auto[v1, v2, v3] = rhs;
    return {u2 * v3 - u3 * v2, u1 * v3 - u3 * v1, u1 * v2 - u2 * v1};
}

template<typename T>
T vector_norm2(const std::tuple<T, T, T> &vector1) {
    auto[u1, u2, u3] = vector1;
    return u1 * u1 + u2 * u2 + u3 * u3;
}

template<typename T>
T vector_norm(const std::tuple<T, T, T> &vector1) {
    return std::sqrt(vector_norm2(vector1));
}

template<typename T>
std::tuple<T, T, T> operator-(const std::tuple<T, T, T> &lhs, const std::tuple<T, T, T> &rhs) {
    return {std::get<0>(lhs) - std::get<0>(rhs),
            std::get<1>(lhs) - std::get<1>(rhs),
            std::get<2>(lhs) - std::get<2>(rhs)};
}

template<typename T>
std::tuple<T, T, T> operator+(const std::tuple<T, T, T> &lhs, const std::tuple<T, T, T> &rhs) {
    return {std::get<0>(lhs) + std::get<0>(rhs),
            std::get<1>(lhs) + std::get<1>(rhs),
            std::get<2>(lhs) + std::get<2>(rhs)};
}

template<typename T>
std::tuple<T, T, T> &operator+=(std::tuple<T, T, T> &lhs, const std::tuple<T, T, T> &rhs) {
    std::get<0>(lhs) += std::get<0>(rhs);
    std::get<1>(lhs) += std::get<1>(rhs);
    std::get<2>(lhs) += std::get<2>(rhs);
    return lhs;
}

template<typename T>
std::tuple<T, T, T> &operator-=(std::tuple<T, T, T> &lhs, const std::tuple<T, T, T> &rhs) {
    std::get<0>(lhs) -= std::get<0>(rhs);
    std::get<1>(lhs) -= std::get<1>(rhs);
    std::get<2>(lhs) -= std::get<2>(rhs);
    return lhs;
}

template<typename T>
std::tuple<T, T, T> operator/(const std::tuple<T, T, T> &vector1, T norm) {
    return {std::get<0>(vector1) / norm,
            std::get<1>(vector1) / norm,
            std::get<2>(vector1) / norm};
}

template<typename T>
std::tuple<T, T, T> &operator/=(std::tuple<T, T, T> &vector1, T norm) {
    std::get<0>(vector1) /= norm;
    std::get<1>(vector1) /= norm;
    std::get<2>(vector1) /= norm;
    return vector1;
}

template<typename T>
std::tuple<T, T, T> operator*(T norm, const std::tuple<T, T, T> &vector1) {
    return {std::get<0>(vector1) * norm,
            std::get<1>(vector1) * norm,
            std::get<2>(vector1) * norm};
}

template<typename T>
std::tuple<T, T, T> operator*(const std::tuple<T, T, T> &vector1, T norm) {
    return norm * vector1;
}

template<typename T>
std::tuple<T, T, T> &operator*=(std::tuple<T, T, T> &vector1, T norm) {
    std::get<0>(vector1) *= norm;
    std::get<1>(vector1) *= norm;
    std::get<2>(vector1) *= norm;
    return vector1;
}


class Grid {
public:
    int x, y, z;
};

template<typename T>
std::string chrono_cast(const T &dur) {
    auto secs = std::chrono::duration_cast<std::chrono::seconds>(dur).count();
    std::string format_string;
    auto hours = secs / 3600;
    secs %= 3600;
    if (hours > 0) {
        format_string = std::to_string(hours) + (hours > 1 ? " hours" : " hour");
    }
    auto mins = secs / 60;
    secs %= 60;
    if (mins > 0) {
        format_string += (format_string.empty() ? "" : " ") + std::to_string(mins) + (mins > 1 ? " mins" : " min");
    }
    if (secs > 0) {
        format_string += (format_string.empty() ? "" : " ") + std::to_string(secs) + (secs > 1 ? " secs" : " sec");
    }
    if (format_string.empty()) {
        format_string = "0 sec";
    }
    return format_string;
}

std::string getOutputFilename(const po::variables_map &vm);

std::string getTopologyFilename(const po::variables_map &vm);

std::string getTrajectoryFilename(const po::variables_map &vm);

std::string getPrmFilename(const po::variables_map &vm);

std::size_t getDefaultVectorReserve();

#endif //TINKER_COMMON_HPP
