
#include "trajectoryreader.hpp"

#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/irange.hpp>
#include <fstream>
#include <iostream>
#include <list>

#include "ReaderFactory.hpp"
#include "config.h"
#include "data_structure/atom.hpp"
#include "data_structure/frame.hpp"
#include "data_structure/molecule.hpp"
#include "nlohmann/json.hpp"
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
            traj_filenames.emplace(std::move(name), start, end);
        }
    } else
        traj_filenames.push(filename);
}

void TrajectoryReader::set_topology(const std::string &filename) { topology_filename = filename; }

std::shared_ptr<Frame> TrajectoryReader::readOneFrame() {
    if (!frame) {
        readTopology();
        if (!mask_string.empty()) {
            auto mask_for_readtraj = parse_atoms(mask_string);
            boost::for_each(frame->atom_list, [this, &mask_for_readtraj](const std::shared_ptr<Atom> &atom) {
                if (Atom::is_match(atom, mask_for_readtraj)) {
                    atoms_for_readtraj.push_back(atom);
                }
            });
        } else {
            atoms_for_readtraj = frame->atom_list;
        }
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
