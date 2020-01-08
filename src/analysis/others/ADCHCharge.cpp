#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/process.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include "ADCHCharge.hpp"
#include "utils/common.hpp"

boost::optional<boost::fusion::vector<std::string, double, double>> ADCHCharge::read_charge(const std::string &line) {
    using namespace boost::spirit::qi;
    //  Atom:    1Am  Corrected charge:    1.164452  Before:    1.124097
    static const auto parser = copy("Atom:" >> omit[int_] >> lexeme[+alpha]
                                            >> "Corrected charge:" >> double_ >> "Before:" >> double_);

    boost::fusion::vector<std::string, double, double> attribute;
    if (auto it = std::begin(line);
            phrase_parse(it, std::end(line), parser, ascii::space, attribute) and it == std::end(line))
        return attribute;

    return {};
}

namespace {
    std::vector<int> read_atoms(std::string_view prompt = "") {
        std::vector<int> atoms;
        for (std::string line;;) {
            using namespace boost::spirit;
            std::cout << prompt;
            std::getline(std::cin, line);
            if (auto it = std::begin(line);
                    qi::phrase_parse(it, std::end(line), +qi::int_, ascii::space, atoms) and it == std::end(line)) {
                break;
            }
            std::cerr << "Syntax Error\n";
            atoms.clear();
        }
        return atoms;
    }
}

void ADCHCharge::process() {
    std::string file;
    std::getline(std::cin, file);
    boost::trim(file);

    auto atoms = read_atoms();
    process(file, atoms);
}

void ADCHCharge::process_interactive(boost::optional<std::string> file) {
    std::string filename = file.has_value() ? *file : choose_file("Enter fchk > ").isExist(true).extension("fchk");
    auto atoms = read_atoms("Atom Numbers > ");
    process(filename, atoms);
}

void ADCHCharge::process(const std::string &file, const std::vector<int> &atoms) {
    namespace bp = boost::process;
    bp::ipstream is;
    bp::opstream os;

    bp::child c(bp::search_path("Multiwfn"), file, bp::std_out > is, bp::std_in < os);

    os << 7 << std::endl;  // Population analysis
    os << 11 << std::endl; // Atomic dipole corrected Hirshfeld population (ADCH)
    os << 1 << std::endl;  // Use build-in sphericalized atomic densities in free-states (more convenient)

    std::string line;
    std::vector<boost::fusion::vector<std::string, double, double>> charges;
    while (c.running() and std::getline(is, line)) {
        if (boost::contains(line, "======= Summary of atomic dipole moment corrected (ADC) charges =======")) {
            while (std::getline(is, line) and !boost::contains(line, "Note:")) {
                if (auto charge = read_charge(line); charge.has_value()) {
                    charges.push_back(*charge);
                }
            }
            break;
        }
    }


    boost::accumulators::accumulator_set<double,
            boost::accumulators::features<boost::accumulators::tag::mean>> adch_acc, hirshfeld_acc;

    std::cout << std::setw(22) << "ADCH" << std::setw(15) << "Hirshfeld" << '\n'
              << std::fixed << std::setprecision(6);
    for (const auto atom: atoms) {
        auto &bv = charges.at(atom - 1);
        auto adch_charge = boost::fusion::at_c<1>(bv);
        adch_acc(adch_charge);
        auto hirshfeld_charge = boost::fusion::at_c<2>(bv);
        hirshfeld_acc(hirshfeld_charge);

        std::cout << std::setw(3) << atom << '(' << std::left << std::setw(2) << boost::fusion::at_c<0>(bv) << ')'
                  << std::right << std::setw(15) << adch_charge << std::setw(15) << hirshfeld_charge << '\n';

    }

    std::cout << std::setw(7) << "Average:"
              << std::setw(14) << boost::accumulators::mean(adch_acc)
              << std::setw(15) << boost::accumulators::mean(hirshfeld_acc) << '\n';

}

