#ifndef TINKER_TRAJECTORYINTERFACE_HPP
#define TINKER_TRAJECTORYINTERFACE_HPP

#include <fstream>
#include <filesystem>
#include "data_structure/frame.hpp"
#include "TopologyInterface.hpp"
#include "utils/common.hpp"
#include "data_structure/atom.hpp"
#include "topology_utils.hpp"

class TrajectoryInterface {

public:
    virtual bool open(const std::string &file) = 0;

    bool readOneFrame(std::shared_ptr<Frame> &frame);;

    virtual void close() = 0;

    virtual ~TrajectoryInterface() = default;

protected:
    virtual bool readOneFrameImpl(std::shared_ptr<Frame> &frame) = 0;

    static void apply_box(std::shared_ptr<Frame> &frame);
};

#endif //TINKER_TRAJECTORYINTERFACE_HPP
