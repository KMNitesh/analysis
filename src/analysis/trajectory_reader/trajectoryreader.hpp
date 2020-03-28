//
// Created by xiamr on 3/17/19.
//

#ifndef TINKER_TRAJECTORYREADER_HPP
#define TINKER_TRAJECTORYREADER_HPP

#include "data_structure/atom.hpp"
#include "trajectory_reader/TopologyInterface.hpp"
#include "trajectory_reader/TrajectoryInterface.hpp"

class TrajectoryReader {
public:
    void add_trajectoy_file(const std::string &filename);

    void set_topology(const std::string &filename);

    std::shared_ptr<Frame> readOneFrame();

    std::shared_ptr<Frame> readTopology();

    void set_mask(std::string mask_string);

private:
    class TrajectoryFile {
    public:
        TrajectoryFile() = default;
        TrajectoryFile(std::string name, std::size_t start = 1, std::size_t end = 0, AmberMask mask = boost::blank{})
            : name(std::move(name)), start(start), end(end), mask(std::move(mask)) {}

        operator std::string() { return name; }

        std::string name;
        std::size_t start;
        std::size_t end;
        AmberMask mask;
    };

    TrajectoryFile current_trajectory_file;
    std::size_t current_frame_pos;

    AmberMask mask;
    std::vector<std::shared_ptr<Atom>> atoms_for_readtraj;

    boost::optional<std::string> topology_filename;
    std::queue<TrajectoryFile> traj_filenames; // the continuous trajectory files

    std::shared_ptr<Frame> frame;
    std::shared_ptr<TrajectoryInterface> traj_reader;
};

#endif // TINKER_TRAJECTORYREADER_HPP
