//
// Created by xiamr on 6/28/19.
//

#ifndef TINKER_TRAJECTORYFORMATWRITER_HPP
#define TINKER_TRAJECTORYFORMATWRITER_HPP

#include "data_structure/frame.hpp"
#include <memory>
#include <string>
#include <vector>

class Atom;

class TrajectoryFormatWriter {
public:
    virtual void open(const std::string &filename) = 0;

    virtual void close() = 0;

    virtual void write(const std::shared_ptr<Frame> &frame, const std::vector<std::shared_ptr<Atom>> &atoms) = 0;

    virtual void write(const std::shared_ptr<Frame> &frame) { write(frame, frame->atom_list); };

    virtual ~TrajectoryFormatWriter() = default;
};

#endif // TINKER_TRAJECTORYFORMATWRITER_HPP
