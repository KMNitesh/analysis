
#include "trajectoryreader.hpp"

#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/irange.hpp>
#include <fstream>
#include <iostream>
#include <list>

#include "ReaderFactory.hpp"
#include "data_structure/atom.hpp"
#include "data_structure/frame.hpp"
#include "data_structure/molecule.hpp"
#include "nlohmann/json.hpp"
#include "utils/PBCUtils.hpp"
#include "utils/common.hpp"

void TrajectoryReader::add_trajectoy_file(const std::string &filename) {
    if (getFileType(filename) == FileType::JSON) {
        std::ifstream i(filename);
        nlohmann::json j;
        i >> j;
        for (auto &item : j) {
            std::size_t start = item.value("start", 1);
            std::size_t end = item.value("end", 0);
            std::string name = item.at("name");
            std::string mask = item.value("mask", "");
            traj_filenames.emplace(std::move(name), start, end,
                                   mask.empty() ? boost::blank{} : parse_atoms(mask, true));
        }
    } else
        traj_filenames.push(filename);
}

void TrajectoryReader::set_topology(const std::string &filename) { topology_filename = filename; }

void TrajectoryReader::set_mask(std::string mask_string) {
    if (!mask_string.empty()) {
        mask = parse_atoms(mask_string, true);
    }
}
std::shared_ptr<Frame> TrajectoryReader::readOneFrame() {
    if (!frame) {
        readTopology();
        atoms_for_readtraj = Atom::isBlank(mask) ? frame->atom_list : PBCUtils::find_atoms(mask, frame);
    }
    for (;;) {
        if (!traj_reader) {
            if (traj_filenames.empty())
                return {};
            current_trajectory_file = traj_filenames.front();
            traj_filenames.pop();
            current_frame_pos = 1;
            traj_reader = ReaderFactory::getTrajectory(current_trajectory_file);
            traj_reader->open(current_trajectory_file);
            if (Atom::isBlank(mask)) {
                atoms_for_readtraj = Atom::isBlank(current_trajectory_file.mask)
                                         ? frame->atom_list
                                         : PBCUtils::find_atoms(current_trajectory_file.mask, frame);
            }
            while (current_frame_pos < current_trajectory_file.start) {
                traj_reader->readOneFrame(frame, atoms_for_readtraj);
                ++current_frame_pos;
            }
        }
        if (traj_reader->readOneFrame(frame, atoms_for_readtraj)) {
            if (current_trajectory_file.end == 0 or current_frame_pos <= current_trajectory_file.end) {
                ++current_frame_pos;
                return frame;
            }
        }
        traj_reader->close();
        traj_reader.reset();
    }
}

std::shared_ptr<Frame> TrajectoryReader::readTopology() {
    std::string file;
    if (topology_filename)
        file = topology_filename.get();
    else {
        if (!traj_filenames.empty() and getFileType(traj_filenames.front()) == FileType::ARC) {
            file = traj_filenames.front();
        } else {
            std::cerr << "Topology file not set !\n";
            std::exit(EXIT_FAILURE);
        }
    }

    auto reader = ReaderFactory::getTopology(file);
    frame = reader->read(file);
    frame->enable_bound = true; // TODO: detect PBC condtion
    for (const auto &element : frame->molecule_list | boost::adaptors::indexed(1)) {
        element.value()->sequence = element.index();
    }
    return frame;
}
