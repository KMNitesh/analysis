
#ifndef TINKER_GRO_WRITER_HPP
#define TINKER_GRO_WRITER_HPP

#include <fstream>
#include <memory>
#include <string>
#include <utility>

#include "FileImpl.hpp"
#include "FileInterface.hpp"
#include "TrajectoryFormatWriter.hpp"

class Frame;

class GROWriter : public TrajectoryFormatWriter {
    std::shared_ptr<FileInterface> os;

public:
    explicit GROWriter(std::shared_ptr<FileInterface> f = std::make_shared<FileImpl>()) : os(std::move(f)) {}

    void open(const std::string &filename) override;

    void close() override;

    using TrajectoryFormatWriter::write;
    
    void write(const std::shared_ptr<Frame> &frame, const std::vector<std::shared_ptr<Atom>> &atoms) override;
};

#endif  // TINKER_GRO_WRITER_HPP
