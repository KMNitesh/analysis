#include <boost/spirit/include/qi.hpp>
#include <boost/phoenix/function/adapt_function.hpp>
#include <boost/xpressive/xpressive_static.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm.hpp>
#include "GromosReader.hpp"
#include "utils/common.hpp"

namespace {
    BOOST_PHOENIX_ADAPT_FUNCTION(void, trim, boost::trim, 1)
}

void GromosReader::process() {
    std::string omd_filename = choose_file("GROMOS omd file > ").isExist(true).extension("omd");

    // mmap all content from file
    boost::iostreams::mapped_file_source file;
    file.open(omd_filename);
    if (!file) {
        std::cerr << "ERROR! GROMOS omd file cannot open \n";
        std::exit(EXIT_FAILURE);
    }

    std::vector<Energy> energies;

    using namespace boost::spirit::qi;
    using namespace boost::phoenix;

    auto timestep_parser = "TIMESTEP" >> eol >>
                                      uint_ >> double_ >> eol
                                      >> "END" >> eol;

    auto energies_parser = "ENERGIES" >> eol
                                      >> +(as_string[lexeme[+(char_ - ':')]][trim(_1)]
                                              >> ':' >> double_ >> eol)
                                      >> eol;

    auto parser = timestep_parser >> energies_parser;

    using namespace boost::xpressive;
    auto rex = "TIMESTEP" >> -*_ >> _n >> _n;

    bool bFilled = false;

    std::vector<std::string> menuStrings;
    std::size_t max_length = 0;

    for (cregex_token_iterator pos(file.begin(), file.end(), rex), end; pos != end; ++pos) {
        auto it = pos->first;
        namespace bp = boost::fusion;
        bp::vector<uint, double, std::vector<bp::vector<std::string, double>>> attribute;
        if (!(phrase_parse(it, pos->second, parser, ascii::space - boost::spirit::eol, attribute) &&
              it == pos->second)) {
            std::cerr << "omd file is ill-formed !\n";
            std::exit(EXIT_FAILURE);
        };

        Energy e;
        e.step = bp::at_c<0>(attribute);
        e.time = bp::at_c<1>(attribute);

        for (auto &bpv : bp::at_c<2>(attribute)) {
            e.energies.push_back(bp::at_c<1>(bpv));
            if (!bFilled) {
                menuStrings.push_back(bp::at_c<0>(bpv));
                max_length = std::max(max_length, bp::at_c<0>(bpv).size());
            }
        }
        bFilled = true;
        energies.push_back(std::move(e));
    }

    printMenu(menuStrings, max_length);

    auto input_parser = +uint_;

    for (;;) {

        std::cout << "Seleect Items > ";
        std::string line;
        std::getline(std::cin, line);

        auto it = begin(line);
        std::vector<uint> attribute;
        if (!(phrase_parse(it, end(line), input_parser, ascii::space, attribute) && it == end(line))) {
            std::cerr << "Paser Error !\n" << line << '\n';
            for (auto iter = line.begin(); iter != it; ++iter) std::cout << " ";
            std::cout << "^\n";
            continue;
        }
        if (boost::find_if(attribute, [&](uint i) { return i < 1 or i > menuStrings.size(); }) != std::end(attribute)) {
            std::cerr << "Number must inside 1.." << menuStrings.size() << '\n';
            continue;
        }
        if (std::unique(std::begin(attribute), std::end(attribute)) != std::end(attribute)) {
            std::cerr << "Number must inside unique\n";
            continue;
        }
        boost::sort(attribute);
        std::ofstream ofs;
        for (;;) {
            std::string outfile = choose_file("Output xvg file > ").isExist(false);
            ofs.open(outfile);
            if (ofs) break;
            std::cerr << "Cannot open file <" << outfile << ">\n";
        }

        ofs << '#' << std::setw(9) << "Time";
        boost::for_each(attribute, [&](uint i) { ofs << std::setw(15) << menuStrings[i - 1]; });
        ofs << '\n';

        for (auto &e : energies) {
            ofs << std::setw(10) << e.time;
            boost::for_each(attribute, [&](uint i) { ofs << std::setw(15) << e.energies[i - 1]; });
            ofs << '\n';
        }
        break;
    }
}

void GromosReader::printMenu(std::vector<std::string> &menuStrings, std::size_t width) {
    std::cout << ">>>  MENU  <<<<\n" << std::left;
    for (const auto &element : menuStrings | boost::adaptors::indexed(1)) {
        std::cout << '(' << element.index() << ") " << std::setw(width + 2) << element.value();
        if (element.index() % 4 == 0) std::cout << '\n';
    }
    if (menuStrings.size() % 4 != 0) std::cout << '\n';
}
