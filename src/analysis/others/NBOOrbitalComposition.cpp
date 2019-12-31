#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/xpressive/xpressive.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/phoenix/function/adapt_function.hpp>
#include <boost/process.hpp>
#include "NBOOrbitalComposition.hpp"
#include "utils/common.hpp"


void NBOOrbitalComposition::process() {
    std::string file = choose_file("NBO output file > ").isExist(true);
    std::ifstream ifstream(file);
    auto alpha_orbitals = findOccupancy(ifstream);
    auto beta_orbitals = findOccupancy(ifstream);
    std::cout << "Occupied molecular alpha orbital number : " << alpha_orbitals << '\n';
    std::cout << "Occupied molecular beta  orbital number : " << beta_orbitals << '\n';

    driveMultiwfn(file, alpha_orbitals, beta_orbitals);
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
    throw std::runtime_error("NBO file syntax errror");
}

bool NBOOrbitalComposition::match(const std::string &line) {
    using namespace boost::xpressive;
    static sregex re = bos >> *space >> +digit >> ". (" >> *_;
    return regex_match(begin(line), end(line), re);
}

namespace {
    BOOST_PHOENIX_ADAPT_FUNCTION(void, trim, boost::trim, 1)
}


void NBOOrbitalComposition::driveMultiwfn(const std::string &file, int alpha_orbitals, int beta_orbitals) {
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

    std::vector<
            boost::fusion::vector<
                    std::vector<
                            boost::fusion::vector<int, boost::optional<int>>
                    >,
                    boost::variant<
                            std::vector<std::string>,
                            std::string
                    >
            >
    > attrs;

    for (std::string line;;) {
        std::cout << "<atom, orbital> : ";
        std::getline(std::cin, line);

        using namespace boost::spirit::qi;
        using namespace boost::phoenix;
        if (auto it = std::begin(line);
                phrase_parse(
                        it, std::end(line),
                        +('<' >> ((int_ >> -('-' >> int_)) % ',') >> ':'
                              >> ((as_string[lexeme[+alnum]] % ',') | string("*")) >> '>'),
                        ascii::space, attrs) and it == std::end(line)) {
            break;
        }
        std::cerr << "Syntax Error !\n";
    }

    auto alpha_contributions = read_contributions(is, alpha_orbitals);
    auto beta_contributions = read_contributions(is, beta_orbitals);


    auto[alpha_filter_contributions, alpha_column_names] = filter(attrs, alpha_contributions);
    print_contributions("Alpha", alpha_filter_contributions, alpha_column_names);

    auto[beta_filter_contributions, beta_column_names] = filter(attrs, beta_contributions);
    print_contributions("Beta", beta_filter_contributions, beta_column_names);

}

std::pair<std::map<int, std::vector<double>, std::greater<>>, std::vector<std::string>>
NBOOrbitalComposition::filter(
        const std::vector<
                boost::fusion::vector<
                        std::vector<boost::fusion::vector<int, boost::optional<int>>>,
                        boost::variant<std::vector<std::string>, std::string>
                >
        > &attrs,
        const std::map<int, std::vector<AtomComposition>, std::greater<>> &contributions) {

    std::map<int, std::vector<double>, std::greater<>> contr;
    std::vector<std::string> column_names;

    auto atom_no_matcher = [](const AtomComposition &atom, auto &bf_vector) -> bool {
        for (auto &v : boost::fusion::at_c<0>(bf_vector)) {
            auto &v1 = boost::fusion::at_c<0>(v);
            auto &v2 = boost::fusion::at_c<1>(v);
            if (v2.has_value()) {
                if (atom.atom_no >= v1 and atom.atom_no <= v2) return true;
            } else {
                if (atom.atom_no == v1) return true;
            }
        }
        return false;
    };

    auto orbital_name_matchaer = [](const AtomComposition &atom, auto &bf_vector) -> bool {

        struct name_visitor : boost::static_visitor<bool> {
            bool operator()(const std::vector<std::string> &bv) const {
                for (auto &v : bv) {
                    if (atom.orbial_name == v) return true;
                }
                return false;
            }

            bool operator()(const std::string &v) const {
                assert(v == "*");
                return true;
            }

            const AtomComposition &atom;

            explicit name_visitor(const AtomComposition &atom) : atom(atom) {}
        } matcher(atom);

        return boost::apply_visitor(matcher, boost::fusion::at_c<1>(bf_vector));
    };

    for (auto &bf_vector : attrs) {
        std::string generated;
        std::back_insert_iterator<std::string> sink{generated};

        using namespace boost::spirit::karma;
        generate(sink, '<' << ((int_ << -('-' << int_)) % ',') << ':' << ((string % ',') | string) << '>', bf_vector);
        column_names.emplace_back(std::move(generated));

        for (auto &[orbital, comp] : contributions) {
            double sum{};
            for (auto &atom : comp) {
                if (orbital_name_matchaer(atom, bf_vector) and atom_no_matcher(atom, bf_vector)) {
                    sum += atom.contribution;
                }
            }
            contr[orbital].emplace_back(sum);
        }
    }
    return {contr, column_names};
}


std::optional<boost::fusion::vector<int, std::string, std::string, double>>
NBOOrbitalComposition::parseLine(const std::string &line) {

    using namespace boost::spirit::qi;
    static const auto line_parser =
            copy(omit[int_] >> int_ >> '(' >> as_string[lexeme[+alpha]] >> ')'
                            >> omit[lexeme[+(char_ - ascii::space)] >> lexeme[+alpha]]
                            >> '(' >> as_string[lexeme[+(char_ - ')')]][trim(_1)] >> ')' >> double_ >> '%');

    boost::fusion::vector<int, std::string, std::string, double> attribute;
    if (auto it = std::begin(line);
            phrase_parse(it, std::end(line), line_parser, ascii::space, attribute) and it == std::end(line)) {
        return attribute;
    }
    return {};
}

void NBOOrbitalComposition::print_contributions(
        std::string_view descriptions,
        std::map<int, std::vector<double>, std::greater<>> &contributions,
        const std::vector<std::string> &column_names) {

    std::cout << "NAO contributions <percentage(%)> for " << descriptions << " orbital\n";
    std::cout << std::setw(10) << "orbital" << std::setw(5) << "n";

    for (auto &name : column_names) {
        std::cout << std::setw(15) << name;
    }
    std::cout << '\n' << std::setprecision(6) << std::fixed;

    int current_shift{};
    for (auto[orbital, contribution] : contributions) {
        std::cout << std::setw(10) << orbital << std::setw(5) << current_shift--;
        for (auto val : contribution) {
            std::cout << std::setw(15) << val;
        }
        std::cout << '\n';
    }
}

std::map<int, std::vector<NBOOrbitalComposition::AtomComposition>, std::greater<>>
NBOOrbitalComposition::read_contributions(std::istream &is, int orbital_number) {

    std::map<int, std::vector<AtomComposition>, std::greater<>> contributions;
    std::string line;
    while (std::getline(is, line)) {
        if (boost::contains(line, "Below are composition of")) {
            while (std::getline(is, line) and !boost::contains(line, "Summing up the compositions listed above:")) {
                if (auto attribute = parseLine(line); attribute.has_value()) {
                    contributions[orbital_number].emplace_back(
                            AtomComposition{boost::fusion::at_c<0>(*attribute),
                                            boost::fusion::at_c<1>(*attribute),
                                            boost::fusion::at_c<2>(*attribute),
                                            boost::fusion::at_c<3>(*attribute)});
                }
            }
            if (--orbital_number == 0) break;
        }
    }
    return contributions;
}




