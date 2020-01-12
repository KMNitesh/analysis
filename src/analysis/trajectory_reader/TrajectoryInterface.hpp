#ifndef TINKER_TRAJECTORYINTERFACE_HPP
#define TINKER_TRAJECTORYINTERFACE_HPP

#include "data_structure/frame.hpp"

class TrajectoryInterface {
public:
    virtual bool open(const std::string &file) = 0;

    bool readOneFrame(std::shared_ptr<Frame> &frame) {
        frame->has_velocity = false;
        return readOneFrameImpl(frame);
    };

    virtual void close() = 0;

    virtual ~TrajectoryInterface() = default;

protected:
    virtual bool readOneFrameImpl(std::shared_ptr<Frame> &frame) = 0;
};

#endif //TINKER_TRAJECTORYINTERFACE_HPP
