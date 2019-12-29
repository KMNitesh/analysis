
#include <boost/process.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/format.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/assign.hpp>
#include "MultiwfnAIMDriver.hpp"
#include "utils/common.hpp"

namespace boost {
    template<>
    struct hash<double MultiwfnAIMDriver::BCP_property::*> {
        size_t operator()(double MultiwfnAIMDriver::BCP_property::* __val) const noexcept {
            return *reinterpret_cast<int *>(&__val); // Cation !!
        }
    };
}

void MultiwfnAIMDriver::readBCP(std::istream &is, std::vector<BCP> &bcp_vector) {
    std::string line;
    auto bcp_it = std::begin(bcp_vector);
    using namespace boost::spirit;
    const auto parser = qi::copy(qi::int_ >> qi::repeat(3)[qi::double_] >> "(3,-1)");
    while (std::getline(is, line), !boost::contains(line, "============ Modify or export found CPs ============")) {
        namespace bf = boost::fusion;
        bf::vector<int, std::vector<double>> attribute;
        if (auto it = std::begin(line);
                qi::phrase_parse(it, std::end(line), parser, ascii::space, attribute)
                and it == std::end(line)) {
            bcp_it->index = bf::at_c<0>(attribute);
            auto &c = bf::at_c<1>(attribute);
            bcp_it->coord = {c[0], c[1], c[2]};
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
            bond_parser = ('(' >> qi::int_ >> ',' >> qi::int_ >> ')')
    [_val = construct<std::pair<int, int>>(boost::spirit::_1, boost::spirit::_2)];

    const auto parser = qi::copy(+bond_parser);
    std::vector<std::pair<int, int>> bonds;

    for (;;) {
        std::string line;
        std::cout << "Bond description [ +(a,b) ] > ";
        std::getline(std::cin, line);

        if (auto it = std::begin(line);
                qi::phrase_parse(it, end(line), parser, ascii::space, bonds) && it == end(line)) {

            process(file, bonds, option_menu());
            return;
        }
        std::cerr << "Syntax Error , input again !\n";
    }
}

void MultiwfnAIMDriver::process(const std::string &options) {

    auto[file, center_index, ligand_atoms] = inputParameter();
    std::vector<std::pair<int, int>> bonds;
    for (auto ligand : ligand_atoms) {
        bonds.emplace_back(center_index, ligand);
    }

    process(file, bonds, parse_options(options));
}

std::vector<double MultiwfnAIMDriver::BCP_property::*>
MultiwfnAIMDriver::parse_options(const std::string &option_string) {
    std::vector<double BCP_property::*> output_field;
    for (auto &s : split(option_string, ",")) {
        if (auto it = property_options.right.find(s); it != property_options.right.end()) {
            output_field.push_back(it->second);
        } else {
            std::cerr << "ERROR !! unknown noption `" << s << "`\n";
            std::exit(EXIT_FAILURE);
        }
    }
    return output_field;
}

void MultiwfnAIMDriver::process(const std::string &file, const std::vector<std::pair<int, int>> &bonds,
                                const std::vector<double BCP_property::*> &output_field) {

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
        while (c.running() and std::getline(is, line)
               and !boost::contains(line, "================ Topology analysis ===============")) {
            for (auto &[property, name] : property_names) {
                if (boost::contains(line, name)) {
                    auto fields = split(line);
                    bcp.property.*property = std::stod(fields[fields.size() - 1]);
                    break;
                }
            }
        }
    }

    print(bcp_vector, output_field);
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


void MultiwfnAIMDriver::print(const std::vector<BCP> &bcp_vector,
                              const std::vector<double BCP_property::*> &output_field) {

    const boost::format fmt("        %2d  <=>  %2d   %E\n");
    for (auto property : output_field) {
        printProperty(fmt, bcp_vector, property_names.left.find(property)->second, property);
    }
}


using boost::bimaps::unordered_set_of;
using boost::assign::list_of;

boost::bimap<unordered_set_of<double MultiwfnAIMDriver::BCP_property::*>,
        unordered_set_of<std::string>> MultiwfnAIMDriver::property_names
        = list_of<decltype(MultiwfnAIMDriver::property_names)::relation>
                (&MultiwfnAIMDriver::BCP_property::Density_of_all_electrons, "Density of all electrons")
                (&MultiwfnAIMDriver::BCP_property::Laplacian_of_electron_density, "Laplacian of electron density")
                (&MultiwfnAIMDriver::BCP_property::Hb, "Energy density E(r) or H(r)")
                (&MultiwfnAIMDriver::BCP_property::Ellipticity, "Ellipticity of electron density");


boost::bimap<unordered_set_of<double MultiwfnAIMDriver::BCP_property::*>,
        unordered_set_of<std::string>> MultiwfnAIMDriver::property_options
        = list_of<decltype(MultiwfnAIMDriver::property_options)::relation>
                (&MultiwfnAIMDriver::BCP_property::Density_of_all_electrons, "da")
                (&MultiwfnAIMDriver::BCP_property::Laplacian_of_electron_density, "lda")
                (&MultiwfnAIMDriver::BCP_property::Hb, "hb")
                (&MultiwfnAIMDriver::BCP_property::Ellipticity, "ell");

std::vector<double MultiwfnAIMDriver::BCP_property::*> MultiwfnAIMDriver::option_menu() {
    std::vector<double BCP_property::*> show_options{
            &BCP_property::Density_of_all_electrons,
            &BCP_property::Laplacian_of_electron_density,
            &BCP_property::Hb,
            &BCP_property::Ellipticity
    };

    using namespace boost::spirit;
    struct num_ : qi::symbols<char, double BCP_property::*> {
    } num__;

    std::cout << ">>>> Property Menu <<<<\n";
    for (const auto &element : show_options | boost::adaptors::indexed(1)) {
        std::cout << "(" << element.index() << ") " << property_names.left.find(element.value())->second << '\n';
        num__.add(std::to_string(element.index()), element.value());
    }

    for (;;) {
        std::cout << "Enter > ";

        std::string line;
        std::getline(std::cin, line);

        std::vector<double BCP_property::*> options;
        if (auto it = begin(line);
                qi::phrase_parse(it, end(line), +num__, ascii::space, options) and it == end(line)) {
            return options;
        }
        std::cerr << "Syntax Error !\n";
    }
}
