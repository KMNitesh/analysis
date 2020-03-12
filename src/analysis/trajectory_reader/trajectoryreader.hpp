//
// Created by xiamr on 3/17/19.
//

#ifndef TINKER_TRAJECTORYREADER_HPP
#define TINKER_TRAJECTORYREADER_HPP

#include "TopologyInterface.hpp"
#include "TrajectoryInterface.hpp"

class TrajectoryReader {
public:
    void add_trajectoy_file(const std::string &filename);

    void set_topology(const std::string &filename);

    std::shared_ptr<Frame> readOneFrame();

    std::shared_ptr<Frame> readTopology();

private:
    class TrajectoryFile {
    public:
        TrajectoryFile() = default;
        TrajectoryFile(std::string name, std::size_t start = 1, std::size_t end = 0)
            : name(std::move(name)), start(start), end(end) {}

        operator std::string() { return name; }

        std::string name;
        std::size_t start;
        std::size_t end;
    };

    TrajectoryFile current_trajectory_file;
    std::size_t current_frame_pos;

    boost::optional<std::string> topology_filename;
    std::queue<TrajectoryFile> traj_filenames; // the continuous trajectory files

    std::shared_ptr<Frame> frame;
    std::shared_ptr<TrajectoryInterface> traj_reader;
};

#endif // TINKER_TRAJECTORYREADER_HPP
