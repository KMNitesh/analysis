#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/xpressive/xpressive.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/phoenix/function/adapt_function.hpp>
#include <boost/process.hpp>
#include "NBOOrbitalComposition.hpp"
#include "utils/common.hpp"


void NBOOrbitalComposition::process() {
    std::string file = choose_file("NBO output file > ").isExist(true);
    std::ifstream ifstream(file);
    auto alpha_orbitals = findOccupancy(ifstream);
    std::cout << "Occupied molecular alpha orbital number : " << alpha_orbitals << '\n';
    std::cout << "Occupied molecular beta  orbital number : " << findOccupancy(ifstream) << '\n';

    driveMultiwfn(file, alpha_orbitals);
}


int NBOOrbitalComposition::findOccupancy(std::istream &is) {
    std::string line;
    int occupancy_number{};
    while (std::getline(is, line)) {
        if (boost::contains(line, "------------------ Lewis ------------------------------------------------------")) {
            while (std::getline(is, line),
                    !boost::contains(line,
                                     "---------------- non-Lewis ----------------------------------------------------")) {
                if (match(line)) ++occupancy_number;
            }
            return occupancy_number;
        }
    }
    throw std::runtime_error("NBO file sytanx errror");
}

bool NBOOrbitalComposition::match(const std::string &line) {
    using namespace boost::xpressive;
    static sregex re = bos >> *space >> +digit >> ". (" >> *_;
    return regex_match(begin(line), end(line), re);
}

namespace {
    BOOST_PHOENIX_ADAPT_FUNCTION(void, trim, boost::trim,

    1)
}


void NBOOrbitalComposition::driveMultiwfn(const std::string &file, int alpha_orbitals) {
    namespace bp = boost::process;
    bp::ipstream is;
    bp::opstream os;

    bp::child c(bp::search_path("Multiwfn"), file, bp::std_out > is, bp::std_in < os);

    os << 8 << std::endl;
    os << 7 << std::endl;

    os << 0 << std::endl;

    for (int orbital = alpha_orbitals; orbital > 0; --orbital) {
        os << orbital << std::endl;
    }

    std::string line;

    using namespace boost::spirit::qi;
    using namespace boost::phoenix;

    std::vector <boost::fusion::vector<int, std::string>> attrs;
    for (;;) {
        std::cout << "<atom, orbital> : ";
        std::getline(std::cin, line);
        if (auto it = std::begin(line);
                phrase_parse(
                        it, std::end(line),
                        +('<' >> int_ >> ',' >> as_string[lexeme[+alnum]] >> '>'),
                        ascii::space, attrs) and it == std::end(line)) {
            break;
        }
        std::cerr << "Syntax Error !\n";
    }

    int current_orbital = alpha_orbitals;
    std::map < int, std::map < std::pair < int, std::string >, double >, std::greater<>>
    contributions;
    namespace bf = boost::fusion;
    while (c.running() and std::getline(is, line)) {
        if (boost::contains(line, "Below are composition of")) {
            while (std::getline(is, line) and !boost::contains(line, "Summing up the compositions listed above:")) {
                if (auto attribute = parseLine(line); attribute.has_value()) {
                    for (auto &bf_vector : attrs) {
                        if (bf::at_c<0>(*attribute) == bf::at_c<0>(bf_vector)
                            and bf::at_c<2>(*attribute) == bf::at_c<1>(bf_vector)) {
                            contributions[current_orbital][{bf::at_c<0>(bf_vector),
                                                            bf::at_c<1>(bf_vector)}] += bf::at_c<3>(*attribute);
                        }
                    }
                }
            }
            if (--current_orbital == 0) break;
        }
    }

    std::cout << "NAO orbital contributions <percentage(%)> \n";
    std::cout << std::setw(10) << "orbital" << std::setw(5) << "n";

    for (auto &bf_vector : attrs) {
        std::cout << std::setw(15)
                  << ("<" + std::to_string(bf::at_c<0>(bf_vector)) + "," + bf::at_c<1>(bf_vector) + ">");
    }
    std::cout << '\n' << std::setprecision(6) << std::fixed;

    int current_shift{};
    for (auto[orbital, contribution] : contributions) {
        std::cout << std::setw(10) << orbital << std::setw(5) << current_shift--;
        for (auto &bf_vector : attrs) {
            std::cout << std::setw(15)
                      << contributions[orbital][{bf::at_c<0>(bf_vector), bf::at_c<1>(bf_vector)}];
        }
        std::cout << '\n';
    }
}


std::optional <boost::fusion::vector<int, std::string, std::string, double>>
NBOOrbitalComposition::parseLine(const std::string &line) {

    using namespace boost::spirit::qi;
    using namespace boost::phoenix;

    static const auto line_parser =
            copy(omit[int_] >> int_ >> '(' >> as_string[lexeme[+(char_ - ')')]][trim(_1)] >> ')'
                            >> omit[lexeme[+(char_ - ascii::space)] >> lexeme[+(char_ - '(')]]
                            >> '(' >> as_string[lexeme[+(char_ - ')')]][trim(_1)] >> ')' >> double_ >> '%');

    boost::fusion::vector<int, std::string, std::string, double> attribute;
    if (auto it = std::begin(line);
            phrase_parse(it, std::end(line), line_parser, ascii::space, attribute) and it == std::end(line)) {
        return attribute;
    }
    return {};
}




