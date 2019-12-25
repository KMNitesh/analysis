
#include <boost/process.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/format.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include "MultiwfnAIMDriver.hpp"
#include "utils/common.hpp"

void MultiwfnAIMDriver::readBCP(std::istream &is, std::vector<BCP> &bcp_vector) {
    std::string line;
    auto bcp_it = std::begin(bcp_vector);
    while (std::getline(is, line), !boost::contains(line, "============ Modify or export found CPs ============")) {
        if (boost::ends_with(line, "(3,-1)")) {
            auto fields = split(line);
            bcp_it->index = std::stoi(fields[0]);
            bcp_it->coord = {std::stoi(fields[1]), std::stoi(fields[2]), std::stoi(fields[3])};
            ++bcp_it;
        }
    }
}

std::tuple<std::string, int, std::vector<int>> MultiwfnAIMDriver::inputParameter() {
    std::string file;
    std::getline(std::cin, file);

    std::string line;

    int center_index;
    std::getline(std::cin, line);
    center_index = std::stoi(line);
    std::vector<int> ligand_atoms;
    std::getline(std::cin, line);

    using namespace boost::spirit;
    qi::phrase_parse(std::begin(line), std::end(line), +qi::int_, ascii::space, ligand_atoms);

    return {file, center_index, ligand_atoms};
}

void MultiwfnAIMDriver::process_interactive() {
    std::string file = choose_file("Fchk file > ").isExist(true).extension("fchk");
    using namespace boost::spirit;
    using boost::phoenix::construct;

    qi::rule<std::string::iterator, std::pair<int, int>(),
            qi::ascii::space_type>
            bond_parser = ('(' >> qi::int_ >> ',' >> qi::int_ >> ')')[_val = construct<std::pair<int, int>>(_1, _2)];

    auto parser = +bond_parser;
    std::vector<std::pair<int, int>> bonds;

    for (;;) {
        std::string line;
        std::cout << "Bond description [ +(a,b) ] > ";
        std::getline(std::cin, line);

        if (auto it = std::begin(line);
                qi::phrase_parse(it, end(line), parser, ascii::space, bonds) && it == end(line)) {
            process(file, bonds);
            return;
        }
        std::cerr << "Syntax Error , input again !\n";
    }
}

void MultiwfnAIMDriver::process() {

    auto[file, center_index, ligand_atoms] = inputParameter();
    std::vector<std::pair<int, int>> bonds;
    for (auto ligand : ligand_atoms) {
        bonds.emplace_back(center_index, ligand);
    }

    process(file, bonds);
}

void MultiwfnAIMDriver::process(const std::string &file, const std::vector<std::pair<int, int>> &bonds) {

    namespace bp = boost::process;
    bp::ipstream is;
    bp::opstream os;

    bp::child c(bp::search_path("Multiwfn"), file, bp::std_out > is, bp::std_in < os);
    os << "2" << std::endl;

    std::vector<BCP> bcp_vector;
    for (auto[i, j] : bonds) {
        os << "1" << std::endl;
        os << i << "," << j << std::endl;

        bcp_vector.emplace_back(i, j);
    }

    os << "-4" << std::endl;
    os << "1" << std::endl;
    os << "0" << std::endl;

    std::string line;
    // Read BCP points
    while (c.running() and std::getline(is, line)) {
        if (boost::contains(line, "Index                    Coordinate                     Type")) {
            readBCP(is, bcp_vector);
            break;
        }
    }

    while (c.running() and std::getline(is, line)) {
        if (boost::contains(line, "================ Topology analysis ===============")) break;
    }

    for (auto &bcp : bcp_vector) {
        os << "7" << std::endl;
        os << bcp.index << std::endl;
        while (c.running() and std::getline(is, line)) {
            auto fields = split(line);
            if (boost::contains(line, "Density of all electrons:")) {
                bcp.property.Density_of_all_electrons = std::stod(fields[fields.size() - 1]);
            } else if (boost::contains(line, "Energy density E(r) or H(r):")) {
                bcp.property.Hb = std::stod(fields[fields.size() - 1]);
            } else if (boost::contains(line, "Laplacian of electron density:")) {
                bcp.property.Laplacian_of_electron_density = std::stod(fields[fields.size() - 1]);
            } else if (boost::contains(line, "Ellipticity of electron density:")) {
                bcp.property.Ellipticity = std::stod(fields[fields.size() - 1]);
            } else if (boost::contains(line, "================ Topology analysis ===============")) {
                break;
            }
        }
    }

    print(bcp_vector);
}

void MultiwfnAIMDriver::printProperty(const boost::format &fmt, const std::vector<BCP> &bcp_vector,
                                      std::string_view name, double BCP_property::* field) {

    boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::mean>> acc;
    std::cout << "\n" << name << '\n' << std::scientific;
    boost::for_each(bcp_vector, [&](const BCP &bcp) {
        std::cout << boost::format(fmt)
                     % bcp.atom_pair_belonging.first % bcp.atom_pair_belonging.second
                     % bcp.property.*field;

        acc(bcp.property.*field);
    });
    std::cout << "           Average    " << boost::accumulators::mean(acc) << '\n';
}

void MultiwfnAIMDriver::print(const std::vector<BCP> &bcp_vector) {

    const boost::format fmt("        %2d  <=>  %2d   %E\n");

    printProperty(fmt, bcp_vector, "Density of all electrons:", &BCP_property::Density_of_all_electrons);
    printProperty(fmt, bcp_vector, "Laplacian of electron density:", &BCP_property::Laplacian_of_electron_density);
    printProperty(fmt, bcp_vector, "Energy density E(r) or H(r):", &BCP_property::Hb);
    printProperty(fmt, bcp_vector, "Ellipticity of electron density:", &BCP_property::Ellipticity);
}
