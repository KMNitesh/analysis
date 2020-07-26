
#include <boost/algorithm/string.hpp>

#include "ana_module/AbstractAnalysis.hpp"

void AbstractAnalysis::setOutFilename(std::string outfilename) {
    this->outfilename = std::move(outfilename);
    boost::trim(this->outfilename);
    if (this->outfilename.empty()) {
        throw std::runtime_error("outfilename cannot empty");
    }
}