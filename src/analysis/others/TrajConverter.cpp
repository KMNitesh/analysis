
#include "TrajConverter.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/spirit/include/qi.hpp>

#include "utils/common.hpp"

int TrajConverter::parse_header(const std::string &line) {
    using namespace boost::spirit::qi;
    static const auto head_parser = copy(lit("block=") >> omit[lexeme[+(alpha | char_('_'))]] >> "records=" >> int_);

    int total_atoms;
    if (auto it = std::begin(line);
        phrase_parse(it, std::end(line), head_parser, ascii::space, total_atoms) and it == std::end(line))
        return total_atoms;

    throw std::runtime_error("file syntax error");
}

boost::fusion::vector<std::string, double, double, double> TrajConverter::parse_atom(const std::string &line) {
    using namespace boost::spirit::qi;
    static const auto atom_parser = copy(as_string[lexeme[+alpha]] >> double_ >> double_ >> double_);

    boost::fusion::vector<std::string, double, double, double> attribute;
    if (auto it = std::begin(line);
        phrase_parse(it, std::end(line), atom_parser, ascii::space, attribute) and it == std::end(line))
        return attribute;

    throw std::runtime_error("file syntax error");
}

void TrajConverter::process() {
    std::string filename = choose_file("trj file (input) > ").extension("trj").isExist(true);
    std::ifstream ifstream(filename);

    std::string pdb_filename = choose_file("pdb file (output) > ").extension("pdb").isExist(false);
    std::ofstream os(pdb_filename);

    std::string line;
    std::getline(ifstream, line);

    int mol_seq = 0;
    int total_atoms{};
    int seq = 0;
    const boost::format fmt{"HETATM%5d %4s %3s  %4d    %8.3f%8.3f%8.3f                      %2s\n"};
    while (!ifstream.eof() and std::getline(ifstream, line)) {
        if (boost::starts_with(line, "block")) {
            total_atoms = parse_header(line);
            os << "MODEL     " << std::setw(4) << ++mol_seq << '\n';
        } else {
            auto attribute = parse_atom(line);
            constexpr double bohr2Ang = 0.529177;
            os << boost::format(fmt) % ++seq % boost::fusion::at_c<0>(attribute) % "" % 0 %
                      (boost::fusion::at_c<1>(attribute) * bohr2Ang) % (boost::fusion::at_c<2>(attribute) * bohr2Ang) %
                      (boost::fusion::at_c<3>(attribute) * bohr2Ang) % boost::fusion::at_c<0>(attribute);
            if (total_atoms == seq) {
                seq = 0;
                os << "TER\nENDMDL\n";
            }
        }
    }
}
