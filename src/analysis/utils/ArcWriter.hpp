#ifndef TINKER_ARCWRITER_HPP
#define TINKER_ARCWRITER_HPP

#include "TrajectoryFormatWriter.hpp"

class ArcWriter : public TrajectoryFormatWriter {
public:
    void open(const std::string &filename) override;

    void close() override;

    void write(const std::shared_ptr<Frame> &frame) override;

private:

    std::ofstream ofs;
};


#endif //TINKER_ARCWRITER_HPP
