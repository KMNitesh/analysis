//
// Created by xiamr on 3/19/19.
//

#ifndef TINKER_TRR_WRITER_HPP
#define TINKER_TRR_WRITER_HPP

#include <memory>
#include <string>

#include "GromacsImpl.hpp"
#include "GromacsInterface.hpp"
#include "TrajectoryFormatWriter.hpp"

class Frame;

class TRRWriter : public TrajectoryFormatWriter {
    gmx::t_fileio *xd = nullptr;
    int step;
    float time;
    bool writeVelocities = false;

public:
    void open(const std::string &filename) override;

    using TrajectoryFormatWriter::write;

    void write(const std::shared_ptr<Frame> &frame, const std::vector<std::shared_ptr<Atom>> &atoms) override;

    void close() override;

    virtual void setWriteVelocities(bool writeVelocities) { this->writeVelocities = writeVelocities; }

    [[nodiscard]] bool isWriteVelocities() const { return writeVelocities; }

    virtual void setCurrentTime(gmx::real t) { this->time = t; }

protected:
    /*
     *  Interface for Mock
     */
    virtual GromacsInterface *getGromacsImpl() {
        static GromacsImpl impl;
        return &impl;
    }
};

#endif  // TINKER_TRR_WRITER_HPP
