
#include "trajectory_reader/GroTrajectoryReader.hpp"

#include <boost/algorithm/string.hpp>
#include <cstring>

#include "data_structure/atom.hpp"
#include "data_structure/frame.hpp"
#include "utils/common.hpp"

bool GroTrajectoryReader::open(const std::string &file) {
    ifs.open(file);
    return static_cast<bool>(ifs);
}

bool GroTrajectoryReader::readOneFrameImpl(std::shared_ptr<Frame> &frame) {
    std::string line;
    std::getline(ifs, frame->title);
    std::getline(ifs, line);

    boost::trim(line);
    if (line.empty()) return false;

    auto total_atom_numbers = std::stoi(line);
    if (total_atom_numbers != frame->atom_list.size()) {
        std::cerr << "atom number from gro file do not match topology\n";
        std::exit(EXIT_FAILURE);
    }

    for (auto &atom : frame->atom_list) {
        std::getline(ifs, line);
        atom->x = 10 * std::stod(line.substr(20, 8));
        atom->y = 10 * std::stod(line.substr(28, 8));
        atom->z = 10 * std::stod(line.substr(36, 8));
    }
    std::getline(ifs, line);

    auto fields = split(line);
    gmx::matrix box;
    std::memset(box, 0, 9 * sizeof(gmx::real));

    box[0][0] = std::stod(fields[0]);
    box[1][1] = std::stod(fields[1]);
    box[2][2] = std::stod(fields[2]);

    // v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y);
    if (fields.size() > 3) {
        box[0][1] = std::stod(fields[3]);
        box[0][2] = std::stod(fields[4]);
        box[1][0] = std::stod(fields[5]);
        box[1][2] = std::stod(fields[6]);
        box[2][0] = std::stod(fields[7]);
        box[2][1] = std::stod(fields[8]);
    }
    frame->box = PBCBox(box);
    return true;
}

void GroTrajectoryReader::close() { ifs.close(); }
