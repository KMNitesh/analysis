#ifndef TINKER_ARCWRITER_HPP
#define TINKER_ARCWRITER_HPP

#include "TrajectoryFormatWriter.hpp"
#include "utils/std.hpp"

class ArcWriter : public TrajectoryFormatWriter {
public:
    void open(const std::string &filename) override;

    void close() override;
    using TrajectoryFormatWriter::write;
    void write(const std::shared_ptr<Frame> &frame, const std::vector<std::shared_ptr<Atom>> &atoms) override;

private:
    std::ofstream ofs;
};

#endif  // TINKER_ARCWRITER_HPP
