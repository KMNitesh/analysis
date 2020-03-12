
#include "GaussianFchkReader.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/range/algorithm.hpp>

#include "utils/common.hpp"

std::vector<double> GaussianFchkReader::readCartesian(std::istream &is) {
    std::string line;
    std::vector<double> cartesians;
    while (!is.eof()) {
        std::getline(is, line);
        if (boost::starts_with(line, "Current cartesian coordinates")) {
            auto fields = split(line);
            auto atom_numbers = boost::lexical_cast<std::size_t>(fields[5]);
            while (cartesians.size() < atom_numbers) {
                std::getline(is, line);
                boost::for_each(split(line), [&](auto &field) { cartesians.push_back(std::stod(field)); });
            }
            return cartesians;
        }
    }
    return cartesians;
}

std::tuple<double, double, double> GaussianFchkReader::getCoordinateOfAtom(std::size_t index,
                                                                           const std::vector<double> &cartesians) {
    if (index == 0 or 3 * index > cartesians.size()) throw std::out_of_range("Index invaild");
    return {cartesians[3 * index - 3], cartesians[3 * index - 2], cartesians[3 * index - 1]};
}
