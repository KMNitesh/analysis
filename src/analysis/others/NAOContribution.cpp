#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/xpressive/xpressive.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/phoenix/function/adapt_function.hpp>
#include <boost/process.hpp>
#include "NAOContribution.hpp"
#include "NBOOrbitalComposition.hpp"
#include "utils/common.hpp"

void NAOContribution::process() {

    std::string file = choose_file("NBO output file > ").isExist(true);

    std::ifstream ifstream(file);
    auto alpha_orbitals = NBOOrbitalComposition::findOccupancy(ifstream);
    auto beta_orbitals = NBOOrbitalComposition::findOccupancy(ifstream);
    std::cout << "Occupied molecular alpha orbital number : " << alpha_orbitals << '\n';
    std::cout << "Occupied molecular beta  orbital number : " << beta_orbitals << '\n';

    namespace bp = boost::process;
    bp::ipstream is;
    bp::opstream os;

    bp::child c(bp::search_path("Multiwfn"), file, bp::std_out > is, bp::std_in < os);

    os << 8 << std::endl; // Orbital composition analysis
    os << 7 << std::endl; // Orbital composition analysis by natural atomic orbital (NAO) method

    os << 0 << std::endl; // Show composition of an orbital

    for (int orbital = alpha_orbitals; orbital > 0; --orbital) {
        os << orbital << std::endl;
    }
    os << 0 << std::endl; // return
    os << 3 << std::endl; // Switch spin type to Beta
    os << 0 << std::endl; // Show composition of an orbital

    for (int orbital = beta_orbitals; orbital > 0; --orbital) {
        os << orbital << std::endl;
    }

    std::vector<boost::fusion::vector<int, boost::optional<int>>> raw_attr;
    for (std::string line;;) {
        std::cout << "atom groups : ";
        std::getline(std::cin, line);

        using namespace boost::spirit::qi;
        using namespace boost::phoenix;

        if (auto it = std::begin(line);
                phrase_parse(it, std::end(line), (int_ >> -('-' >> int_)) % ',',
                             ascii::space, raw_attr) and it == std::end(line)) {
            break;
        }
        std::cerr << "Syntax Error !\n";
    }

    std::vector<int> attrs;

    for (auto &v : raw_attr) {
        auto left = boost::fusion::at_c<0>(v);
        auto right = boost::fusion::at_c<1>(v);

        if (right.has_value()) {
            for (int i = left; i <= right; ++i) attrs.push_back(i);
        } else {
            attrs.push_back(left);
        }
    }

    auto alpha_contributions = read_contributions(attrs, is, alpha_orbitals);
    auto beta_contributions = read_contributions(attrs, is, beta_orbitals);

    print_contributions("Alpha", alpha_contributions);
    print_contributions("Beta", beta_contributions);
}

std::map<int, double, std::greater<>> NAOContribution::read_contributions(const std::vector<int> &attrs,
                                                                          std::istream &is,
                                                                          int orbital_number) {

    std::map<int, double, std::greater<>> contributions;

    namespace bf = boost::fusion;
    std::string line;
    while (std::getline(is, line)) {
        if (boost::contains(line, "Condensed above result to atoms:")) {
            while (std::getline(is, line) and !boost::contains(line, "Analyze which orbital?")) {
                if (auto attribute = parseLine(line); attribute.has_value()) {
                    for (auto atom : attrs) {
                        if (bf::at_c<0>(*attribute) == atom) {
                            contributions[orbital_number] += bf::at_c<1>(*attribute);
                        }
                    }
                }
            }
            if (--orbital_number == 0) break;
        }
    }
    return contributions;
}

void NAOContribution::print_contributions(
        std::string_view descriptions,
        std::map<int, double, std::greater<>> &contributions) {

    std::cout << "NAO contributions <percentage(%)> for " << descriptions << " orbital\n";
    std::cout << std::setw(10) << "orbital" << std::setw(5) << "n";

    std::cout << '\n' << std::setprecision(6) << std::fixed;

    int current_shift{};
    for (auto[orbital, contribution] : contributions) {
        std::cout << std::setw(10) << orbital << std::setw(5) << current_shift--;
        std::cout << std::setw(15) << contribution;
        std::cout << '\n';
    }
}


boost::optional<boost::fusion::vector<int, double> > NAOContribution::parseLine(const std::string &line) {

    using namespace boost::spirit::qi;
    using namespace boost::phoenix;

    static const auto line_parser = copy(int_ >> '(' >> omit[lexeme[+alpha]] >> ')' >> double_ >> '%');

    boost::fusion::vector<int, double> attribute;
    if (auto it = std::begin(line);
            phrase_parse(it, std::end(line), line_parser, ascii::space, attribute) and it == std::end(line)) {
        return attribute;
    }
    return {};
}
