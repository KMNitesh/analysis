#include "TrajectoryInterface.hpp"

#include "data_structure/frame.hpp"


bool TrajectoryInterface::readOneFrame(std::shared_ptr<Frame> &frame, const std::vector<std::shared_ptr<Atom>> &atoms) {
    frame->has_velocity = false;
    auto ok = readOneFrameImpl(frame, atoms);
    return ok;
}
