#include "data_structure/frame.hpp"
#include "TrajectoryInterface.hpp"

void TrajectoryInterface::apply_box(std::shared_ptr<Frame> &frame) {
    if (frame->enable_bound) {
        frame->a_axis_half = frame->a_axis / 2;
        frame->b_axis_half = frame->b_axis / 2;
        frame->c_axis_half = frame->c_axis / 2;
    }
}

bool TrajectoryInterface::readOneFrame(std::shared_ptr<Frame> &frame) {
    frame->has_velocity = false;
    auto ok = readOneFrameImpl(frame);
    if (ok) apply_box(frame);
    return ok;
}
