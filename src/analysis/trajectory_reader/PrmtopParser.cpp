
#include "utils/common.hpp"
#include "trajectory_reader/PrmtopGrammar.hpp"

#include <printf.h>

#include <boost/spirit/repository/include/qi_iter_pos.hpp>
#include <boost/spirit/include/support_line_pos_iterator.hpp>
#include "PrmtopParser.hpp"


boost::optional<PrmtopStruct> PrmtopParser::parse(std::istream &is) {
    std::string content{std::istreambuf_iterator<char>(is), std::istreambuf_iterator<char>()};
    using namespace boost::spirit;
    PrmtopGrammar<line_pos_iterator<std::string::iterator>, decltype(ascii::space - qi::eol)> parser;

    line_pos_iterator<std::string::iterator> begin{std::begin(content)}, it{begin}, end;
    try {
        PrmtopStruct prmtop_struct;
        if (qi::phrase_parse(it, end, parser, ascii::space - qi::eol, prmtop_struct)) {
            return prmtop_struct;
        }
    } catch (const qi::expectation_failure<line_pos_iterator<std::string::iterator>> &x) {
        std::cerr << "Grammar Parser Failure ! Expecting : " << x.what_ << '\n';
        auto line = get_line(x.first);
        auto column = get_column(begin, x.first);
        std::string pos = " (line: " + std::to_string(line) + ", column: " + std::to_string(column) + ") ";
        std::cerr << pos << ">>>>" << get_current_line(begin, x.first, x.last) << "<<<<\n";
        std::cerr << std::string(column + pos.size() + 3, ' ') << "^~~~ here\n";
    }
    return {};
}
