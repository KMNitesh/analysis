
#include "config.h"
#include <list>
#include <iostream>
#include <fstream>
#include <boost/range/adaptors.hpp>
#include <boost/range/irange.hpp>

#include "utils/common.hpp"
#include "trajectoryreader.hpp"
#include "data_structure/atom.hpp"
#include "data_structure/molecule.hpp"
#include "data_structure/frame.hpp"
#include "ReaderFactory.hpp"


void TrajectoryReader::add_trajectoy_file(const std::string &filename) {
    traj_filenames.push(filename);
}

void TrajectoryReader::set_topology(const std::string &filename) {
    topology_filename = filename;
}

std::shared_ptr<Frame> TrajectoryReader::readOneFrame() {
    if (!frame) readTopology();
    for (;;) {
        if (!traj_reader) {
            if (traj_filenames.empty()) return {};
            auto file = traj_filenames.front();
            traj_filenames.pop();
            traj_reader = ReaderFactory::getTrajectory(file);
            traj_reader->open(file);
        }
        if (traj_reader->readOneFrame(frame)) {
            return frame;
        }
        traj_reader->close();
        traj_reader.reset();
    }
}

std::shared_ptr<Frame> TrajectoryReader::readTopology() {

    std::string file;
    if (topology_filename) file = topology_filename.get();
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
    frame->enable_bound = true; //TODO: detect PBC condtion
    for (const auto &element : frame->molecule_list | boost::adaptors::indexed(1)) {
        element.value()->sequence = element.index();
    }
    return frame;
}
