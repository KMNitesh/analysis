#include <utility>

//
// Created by xiamr on 3/19/19.
//

#ifndef TINKER_GRO_WRITER_HPP
#define TINKER_GRO_WRITER_HPP

#include <string>
#include <memory>
#include <fstream>

#include "TrajectoryFormatWriter.hpp"
#include "FileInterface.hpp"
#include "FileImpl.hpp"

class Frame;


class GROWriter : public TrajectoryFormatWriter {

    std::shared_ptr<FileInterface> os;
public:
    explicit GROWriter(std::shared_ptr<FileInterface> f = std::make_shared<FileImpl>()) : os(std::move(f)) {}

    void open(const std::string &filename) override;

    void close() override;

    void write(const std::shared_ptr<Frame> &frame) override;

};


#endif //TINKER_GRO_WRITER_HPP
