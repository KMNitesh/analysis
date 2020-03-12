#include "TrajectoryInterface.hpp"

#include "data_structure/frame.hpp"

void TrajectoryInterface::apply_box(std::shared_ptr<Frame> &frame) {}

bool TrajectoryInterface::readOneFrame(std::shared_ptr<Frame> &frame) {
    frame->has_velocity = false;
    auto ok = readOneFrameImpl(frame);
    if (ok) apply_box(frame);
    return ok;
}
