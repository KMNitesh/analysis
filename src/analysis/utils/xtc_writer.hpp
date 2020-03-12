#ifndef TINKER_XTC_WRITER_HPP
#define TINKER_XTC_WRITER_HPP

#include <memory>
#include <string>

#include "GromacsImpl.hpp"
#include "GromacsInterface.hpp"
#include "TrajectoryFormatWriter.hpp"

class Frame;

class XTCWriter : public TrajectoryFormatWriter {
    gmx::t_fileio *xd = nullptr;
    int step;
    gmx::real prec, time;

public:
    void open(const std::string &filename) override;

    void write(const std::shared_ptr<Frame> &frame) override;

    void close() override;

protected:
    /*
     *  Interface for Mock
     */
    virtual GromacsInterface *getGromacsImpl() {
        static GromacsImpl impl;
        return &impl;
    }
};

#endif  // TINKER_XTC_WRITER_HPP
