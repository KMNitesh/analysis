#include "QMStructureComp.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/regex.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/irange.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/xpressive/xpressive_static.hpp>

#include "ana_module/RMSDCal.hpp"
#include "utils/common.hpp"

void QMStructureComp::process() {
    std::string _47file = choose_file("Enter .47 file > ").isExist(true).extension("47");
    std::string log_file = choose_file("Enter log fle > ").isExist(true).extension("log");

    try {
        std::ifstream ifs_47file{_47file};
        auto _47_coord = read_47_file(ifs_47file);

        std::ifstream ifs_log_file{log_file};
        auto _log_file = read_log_file(ifs_log_file);

        if (_47_coord.size() != _log_file.size()) {
            std::cerr << "Structure is not match !\n";
            std::exit(EXIT_FAILURE);
        }

        auto n = _47_coord.size();

        for (auto i : boost::irange(n)) {
            if (element_table[boost::fusion::at_c<0>(_47_coord[i]) - 1] != boost::fusion::at_c<0>(_log_file[i])) {
                std::cerr << "Structure is not match ! for atom  " << i + 1 << '\n';
                std::exit(EXIT_FAILURE);
            }
        }

        double x1[n], y1[n], z1[n];
        double x2[n], y2[n], z2[n];

        for (const auto &l : _47_coord | boost::adaptors::indexed()) {
            x1[l.index()] = boost::fusion::at_c<1>(l.value());
            y1[l.index()] = boost::fusion::at_c<2>(l.value());
            z1[l.index()] = boost::fusion::at_c<3>(l.value());
        }

        for (const auto &l : _log_file | boost::adaptors::indexed()) {
            x2[l.index()] = boost::fusion::at_c<1>(l.value());
            y2[l.index()] = boost::fusion::at_c<2>(l.value());
            z2[l.index()] = boost::fusion::at_c<3>(l.value());
        }

        double mid[3];

        RMSDCal::center(n, x1, y1, z1, mid, n);
        RMSDCal::center(n, x2, y2, z2, mid, n);

        RMSDCal::quatfit(n, x1, y1, z1, n, x2, y2, z2, n);

        auto rms_avg = RMSDCal::rmsfit(x1, y1, z1, x2, y2, z2, n);
        auto rms_max = RMSDCal::rms_max(x1, y1, z1, x2, y2, z2, n);

        std::cout << "RMS AVG = " << rms_avg << " RMS MAX = " << rms_max << '\n';
        if (rms_avg < 0.01 and rms_max < 0.02) std::cout << "Structure is Same! Congratulations @^_^@\n";

    } catch (std::exception &e) {
        std::cerr << "ERROR !! " << e.what() << '\n';
        std::exit(EXIT_FAILURE);
    }
}

std::vector<boost::fusion::vector<uint, double, double, double>> QMStructureComp::read_47_file(std::istream &is) {
    std::vector<boost::fusion::vector<uint, double, double, double>> coord;

    std::string line;
    using namespace boost::spirit::qi;
    const auto parser = copy(uint_ >> omit[uint_] >> double_ >> double_ >> double_);
    while (std::getline(is, line)) {
        if (line == " $COORD") {
            std::getline(is, line);
            while (std::getline(is, line)) {
                boost::fusion::vector<uint, double, double, double> attr;
                if (auto it = std::begin(line);
                    !(phrase_parse(it, std::end(line), parser, ascii::space, attr) and it == std::end(line))) {
                    return coord;
                }
                coord.push_back(std::move(attr));
            }
        }
    }
    throw std::runtime_error("Parse .47 file error!");
}

std::vector<boost::fusion::vector<std::string, double, double, double>> QMStructureComp::read_log_file(
    std::istream &is) {
    std::vector<std::string> line_buffer;
    using namespace boost::xpressive;
    // (Enter /public/home/xiamr/prog/g09-c01/g09/l9999.exe)
    sregex line_loc = " (Enter " >> +_ >> "l9999.exe)";
    std::string line;
    while (std::getline(is, line)) {
        if (regex_match(line, line_loc)) {
            line_buffer.clear();
            continue;
        }
        line_buffer.push_back(line);
    }

    std::string contracte_line;
    for (auto &line : line_buffer) {
        if (line.empty()) break;
        contracte_line += line.substr(1, line.size());
    }

    std::vector<std::string> fields;
    boost::split_regex(fields, contracte_line, boost::regex(R"=(\\\\)="));
    line = fields.at(3);
    using namespace boost::spirit::qi;
    const auto parser =
        copy(omit[double_ >> double_] >> +(R"=(\)=" >> lexeme[+ascii::alpha] >> double_ >> double_ >> double_));
    std::vector<boost::fusion::vector<std::string, double, double, double>> coord;

    if (auto it = std::begin(line);
        phrase_parse(it, std::end(line), parser, char_(','), coord) and it == std::end(line)) {
        return coord;
    }

    throw std::runtime_error("Parse gaussian log file error!");
}
