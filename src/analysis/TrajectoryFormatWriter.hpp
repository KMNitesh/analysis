//
// Created by xiamr on 6/28/19.
//

#ifndef TINKER_TRAJECTORYFORMATWRITER_HPP
#define TINKER_TRAJECTORYFORMATWRITER_HPP

#include <string>
#include <memory>

class Frame;

class TrajectoryFormatWriter {
public:

    virtual void open(const std::string &filename) = 0;

    virtual void close() = 0;

    virtual void write(const std::shared_ptr<Frame> &frame) = 0;

    virtual ~TrajectoryFormatWriter() = default;
};

#endif //TINKER_TRAJECTORYFORMATWRITER_HPP
