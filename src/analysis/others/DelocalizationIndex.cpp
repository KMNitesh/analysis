#include "DelocalizationIndex.hpp"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/process.hpp>
#include <boost/spirit/include/qi.hpp>
#include <iomanip>
#include <iostream>
#include <utility>

#include "MultiwfnAIMDriver.hpp"
#include "utils/common.hpp"

void DelocalizationIndex::process() {
    auto [file, center_index, ligand_atoms] = MultiwfnAIMDriver::inputParameter();
    std::vector<std::pair<int, int>> bonds;
    for (auto ligand : ligand_atoms) {
        bonds.emplace_back(center_index, ligand);
    }

    process(file, bonds);
}

void DelocalizationIndex::process(const std::string &file, const std::vector<std::pair<int, int>> &bonds) {
    namespace bp = boost::process;
    bp::ipstream is;
    bp::opstream os;

    bp::child c(bp::search_path("Multiwfn"), file, bp::std_out > is, bp::std_in < os);
    os << 9 << std::endl;  // Population analysis
    os << 7 << std::endl;  // Fuzzy bond order analysis

    using namespace boost::spirit::qi;
    using namespace boost::phoenix;

    const auto atom_parser = copy(int_ >> '(' >> as_string[lexeme[+alpha]] >> ')');

    const auto parser = copy('#' >> omit[int_] >> ':' >> atom_parser >> atom_parser >> "Alpha:" >> omit[double_] >>
                             "Beta:" >> omit[double_] >> "Total:" >> double_);

    std::string line;
    struct BondOrder {
        BondOrder(std::string symbol1, std::string symbol2, double bondorder)
            : symbol1(std::move(symbol1)), symbol2(std::move(symbol2)), bondorder(bondorder) {}

        std::string symbol1;
        std::string symbol2;
        double bondorder;
    };

    std::map<std::pair<int, int>, BondOrder> bondorders;
    while (c.running() and std::getline(is, line) and !boost::contains(line, "If output bond order matrix?")) {
        boost::fusion::vector<int, std::string, int, std::string, double> bondorder;
        if (auto it = std::begin(line);
            phrase_parse(it, std::end(line), parser, ascii::space, bondorder) and it == std::end(line)) {
            bondorders.insert(
                {{boost::fusion::at_c<0>(bondorder), boost::fusion::at_c<2>(bondorder)},
                 BondOrder{std::move(boost::fusion::at_c<1>(bondorder)), std::move(boost::fusion::at_c<3>(bondorder)),
                           boost::fusion::at_c<4>(bondorder)}});
        }
    }

    std::cout << "Delocalization Index based on Fuzzy Space\n" << std::fixed << std::setprecision(6);
    boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::mean>> acc;
    for (auto [atom1, atom2] : bonds) {
        if (atom1 > atom2) std::swap(atom1, atom2);
        const auto &bondorder = bondorders.at({atom1, atom2});
        std::cout << std::setw(5) << std::right << atom1 << '(' << std::left << std::setw(2) << bondorder.symbol1 << ')'
                  << std::setw(5) << std::right << atom2 << '(' << std::left << std::setw(2) << bondorder.symbol2 << ')'
                  << std::right << std::setw(15) << bondorder.bondorder << '\n';
        acc(bondorder.bondorder);
    }
    std::cout << "         Average  " << std::setw(15) << boost::accumulators::mean(acc) << '\n';
}

void DelocalizationIndex::process_interactive() {
    auto [file, bonds] = MultiwfnAIMDriver::user_input();
    process(file, bonds);
}
