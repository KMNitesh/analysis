//
// Created by xiamr on 3/17/19.
//

#ifndef TINKER_TRAJECTORYREADER_HPP
#define TINKER_TRAJECTORYREADER_HPP

#include "data_structure/atom.hpp"
#include "trajectory_reader/TopologyInterface.hpp"
#include "trajectory_reader/TrajectoryInterface.hpp"
#include <boost/fusion/sequence.hpp>
#include <boost/optional.hpp>

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
        TrajectoryFile(std::string name, AmberMask mask = boost::blank{})
            : name(std::move(name)), mask(std::move(mask)) {}

        operator std::string() { return name; }

        std::string name;
        using Ranges = std::vector<boost::fusion::vector<uint, boost::optional<boost::variant<uint, char>>>>;
        Ranges range;

        struct LessEqual : public boost::static_visitor<bool> {
            LessEqual(uint pos) : pos(pos) {}
            bool operator()(const uint end) const { return pos <= end; };

            bool operator()(const char) const { return true; };

        private:
            uint pos;
        };

        AmberMask mask;

        bool is_in_range(uint pos) const;

        bool is_end(uint pos) const;

        static Ranges parse_range(const std::string &range_string);
    };

    TrajectoryFile current_trajectory_file;
    uint current_frame_pos;

    AmberMask mask;
    std::vector<std::shared_ptr<Atom>> atoms_for_readtraj;

    boost::optional<std::string> topology_filename;
    std::queue<TrajectoryFile> traj_filenames; // the continuous trajectory files

    std::shared_ptr<Frame> frame;
    std::shared_ptr<TrajectoryInterface> traj_reader;
};

#endif // TINKER_TRAJECTORYREADER_HPP
