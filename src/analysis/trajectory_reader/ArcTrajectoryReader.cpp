
#include "ArcTrajectoryReader.hpp"

#include <filesystem>

#include "TrajectoryInterface.hpp"
#include "data_structure/atom.hpp"
#include "data_structure/frame.hpp"
#include "topology_utils.hpp"
#include "utils/common.hpp"

bool ArcTrajectoryReader::open(const std::string &file) {
    ifs.open(file);
    if (enable_read_velocity) {
        std::filesystem::path path(file);
        auto vel = path.parent_path() / (path.stem().string() + ".vel");
        velocity_file.open(vel);
        if (!velocity_file) {
            std::cerr << "ERROR !! Could not open velocity file " << vel << '\n';
            std::exit(EXIT_FAILURE);
        }
    }
    return static_cast<bool>(ifs);
}

bool ArcTrajectoryReader::readOneFrameImpl(std::shared_ptr<Frame> &frame) {
    std::getline(ifs, line);
    if (line.empty()) return false;

    if (frame->enable_bound) {
        std::getline(ifs, line);
        field = split(line);
        if (field.empty()) return false;
        parse_box(frame);
    }

    for (auto &atom : frame->atom_list) {
        std::getline(ifs, line);
        field = split(line);
        parse_coord(atom);
    }
    if (enable_read_velocity) readOneFrameVelocity(frame);

    return true;
}

std::shared_ptr<Frame> ArcTrajectoryReader::read(const std::string &filename) {
    open(filename);

    int atom_num = 0;
    auto frame = std::make_shared<Frame>();
    std::getline(ifs, line);
    field = split(line);
    atom_num = std::stoi(field[0]);
    frame->title = line.substr(line.rfind(field[0]) + field[0].size());
    boost::trim(frame->title);

    if (frame->enable_bound) {
        std::getline(ifs, line);
        field = split(line);
        if (field.empty()) return {};
        parse_box(frame);
    }

    for (int i = 0; i < atom_num; i++) {
        std::getline(ifs, line);
        field = split(line);
        if (i == 0) {
            if (field.empty()) return {};
            try {
                parse_box(frame);
                frame->enable_bound = true;
                i--;
                continue;
            } catch (...) {
            }
        }

        auto atom = std::make_shared<Atom>();
        atom->seq = std::stoi(field[0]);
        atom->atom_name = field[1];
        parse_coord(atom);
        atom->typ = std::stoi(field[5]);
        for (size_t j = 6; j < field.size(); j++) {
            atom->con_list.push_back(std::stoi(field[j]));
        }
        frame->atom_list.push_back(atom);
        frame->atom_map[atom->seq] = atom;
    }

    topology_utils::assgin_atom_to_molecule(frame);
    apply_box(frame);
    frame->build_graph();
    return frame;
}

void ArcTrajectoryReader::parse_coord(std::shared_ptr<Atom> &atom) {
    atom->x = std::stod(field[2]);
    atom->y = std::stod(field[3]);
    atom->z = std::stod(field[4]);
}

void ArcTrajectoryReader::parse_box(std::shared_ptr<Frame> &frame) {
    frame->box = PBCBox(std::stod(field[0]), std::stod(field[1]), std::stod(field[2]), std::stod(field[3]),
                        std::stod(field[4]), std::stod(field[5]));
}

void ArcTrajectoryReader::readOneFrameVelocity(std::shared_ptr<Frame> &frame) {
    std::getline(velocity_file, line);
    for (auto &atom : frame->atom_list) {
        std::getline(velocity_file, line);
        field = split(line);
        auto len1 = field[2].length();
        auto len2 = field[3].length();
        auto len3 = field[4].length();
        atom->vx = std::stod(field[2].replace(len1 - 4, 1, "E"));
        atom->vy = std::stod(field[3].replace(len2 - 4, 1, "E"));
        atom->vz = std::stod(field[4].replace(len3 - 4, 1, "E"));
    }
    frame->has_velocity = true;
}
