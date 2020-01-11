#ifndef TINKER_TRAJECTORYINTERFACE_HPP
#define TINKER_TRAJECTORYINTERFACE_HPP

#include <memory>

class Frame;

class TrajectoryInterface {
public:
    virtual bool open(const std::string &file) = 0;

    virtual bool readOneFrame(std::shared_ptr<Frame> &frame) = 0;

    virtual void close() = 0;

    virtual ~TrajectoryInterface() = default;
};

#endif //TINKER_TRAJECTORYINTERFACE_HPP
