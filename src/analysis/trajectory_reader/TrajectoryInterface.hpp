#ifndef TINKER_TRAJECTORYINTERFACE_HPP
#define TINKER_TRAJECTORYINTERFACE_HPP

#include <filesystem>
#include <fstream>

#include "TopologyInterface.hpp"
#include "data_structure/atom.hpp"
#include "data_structure/frame.hpp"
#include "topology_utils.hpp"
#include "utils/common.hpp"

class TrajectoryInterface {
public:
    virtual bool open(const std::string &file) = 0;

    bool readOneFrame(std::shared_ptr<Frame> &frame, const std::vector<std::shared_ptr<Atom>> &atoms);

    virtual void close() = 0;

    virtual ~TrajectoryInterface() = default;

protected:
    virtual bool readOneFrameImpl(std::shared_ptr<Frame> &frame, const std::vector<std::shared_ptr<Atom>> &atoms) = 0;

    static void apply_box(std::shared_ptr<Frame> &frame);
};

#endif  // TINKER_TRAJECTORYINTERFACE_HPP
