#include "PrmtopParser.hpp"
#include "trajectory_reader/PrmtopGrammar.hpp"
#include "utils/common.hpp"


boost::optional<PrmtopStruct> PrmtopParser::parse(std::istream &is) {
    is.unsetf(std::ios::skipws);

    using namespace boost::spirit;
    PrmtopGrammar<istream_iterator, decltype(ascii::space - qi::eol)> parser;

    PrmtopStruct prmtop_struct;
    if (istream_iterator it(is), end;
            qi::phrase_parse(it, end, parser, ascii::space - qi::eol, prmtop_struct) and it == end) {
        return prmtop_struct;
    }
    return {};
}


