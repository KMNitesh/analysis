//
// Created by xiamr on 7/24/19.
//

#include "TinkerDynReader.hpp"

#include <boost/algorithm/string/predicate.hpp>

#include "utils/ThrowAssert.hpp"
#include "utils/common.hpp"

void TinkerDynReader::readContent() {
    std::string line;
    std::getline(*istream, line);
    std::getline(*istream, line);
    for (;;) {
        std::getline(*istream, line);
        if (boost::contains(line, " Periodic Box Dimensions :")) {
            std::getline(*istream, line);
            auto fields = split(line);
            throw_assert(fields.size() == 3, "Periodic Box Dimensions Information in dyn file is invalid!");
            box.xbox = std::stod(boost::replace_all_copy(fields[0], "D", "E"));
            box.ybox = std::stod(boost::replace_all_copy(fields[1], "D", "E"));
            box.zbox = std::stod(boost::replace_all_copy(fields[2], "D", "E"));
            std::getline(*istream, line);
            fields = split(line);
            throw_assert(fields.size() == 3, "Periodic Box Dimensions Information in dyn file is invalid!");
            box.alpha = std::stod(boost::replace_all_copy(fields[0], "D", "E"));
            box.beta = std::stod(boost::replace_all_copy(fields[1], "D", "E"));
            box.gamma = std::stod(boost::replace_all_copy(fields[2], "D", "E"));
        } else if (boost::contains(line, " Current Atomic Positions :")) {
            for (;;) {
                std::getline(*istream, line);
                if (boost::contains(line, " Current Atomic Velocities :")) {
                    return;
                }
                auto fields = split(line);
                throw_assert(fields.size() == 3, "Atomic Positions Information in dyn file is invalid!");
                positions.coordinates.emplace_back(std::stod(boost::replace_all_copy(fields[0], "D", "E")),
                                                   std::stod(boost::replace_all_copy(fields[1], "D", "E")),
                                                   std::stod(boost::replace_all_copy(fields[2], "D", "E")));
            }
        }
    }
}
