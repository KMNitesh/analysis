
#include "utils/common.hpp"
#include "data_structure/atom.hpp"
#include "data_structure/frame.hpp"
#include "topology_utils.hpp"
#include "ArcTrajectoryReader.hpp"

bool ArcTrajectoryReader::open(const std::string &file) {
    ifs.open(file);
    return static_cast<bool>(ifs);
}


bool ArcTrajectoryReader::readOneFrameImpl(std::shared_ptr<Frame> &frame) {
    std::getline(ifs, line);
    if (line.empty()) return false;

    if (frame->enable_bound) {
        std::getline(ifs, line);
        field = split(line);
        if (field.empty()) return false;
        parse_box(frame, field);
    }

    for (auto &atom : frame->atom_list) {
        std::getline(ifs, line);
        field = split(line);
        parse_coord(atom, field);
    }
    apply_box(frame);

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
        parse_box(frame, field);
    }

    for (int i = 0; i < atom_num; i++) {
        std::getline(ifs, line);
        field = split(line);
        if (i == 0) {
            if (field.empty()) return {};
            try {
                parse_box(frame, field);
                frame->enable_bound = true;
                i--;
                continue;
            } catch (...) {}
        }

        auto atom = std::make_shared<Atom>();
        atom->seq = std::stoi(field[0]);
        atom->atom_name = field[1];
        parse_coord(atom, field);
        atom->typ = std::stoi(field[5]);
        for (size_t j = 6; j < field.size(); j++) {
            atom->con_list.push_back(std::stoi(field[j]));
        }
        frame->atom_list.push_back(atom);
        frame->atom_map[atom->seq] = atom;
    }

    topology_utils::assgin_atom_to_molecule(frame);
    apply_box(frame);

    return frame;
}

void ArcTrajectoryReader::apply_box(std::shared_ptr<Frame> &frame) {
    if (frame->enable_bound) {
        frame->a_axis_half = frame->a_axis / 2;
        frame->b_axis_half = frame->b_axis / 2;
        frame->c_axis_half = frame->c_axis / 2;
    }
}

void ArcTrajectoryReader::parse_coord(std::shared_ptr<Atom> &atom, const std::vector<std::string> &field) {
    atom->x = std::stod(field[2]);
    atom->y = std::stod(field[3]);
    atom->z = std::stod(field[4]);
}

void ArcTrajectoryReader::parse_box(std::shared_ptr<Frame> &frame, const std::vector<std::string> &field) {
    frame->a_axis = std::stod(field[0]);
    frame->b_axis = std::stod(field[1]);
    frame->c_axis = std::stod(field[2]);
    frame->alpha = std::stod(field[3]);
    frame->beta = std::stod(field[4]);
    frame->gamma = std::stod(field[5]);
}
