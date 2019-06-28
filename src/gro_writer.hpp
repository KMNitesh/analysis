//
// Created by xiamr on 3/19/19.
//

#ifndef TINKER_GRO_WRITER_HPP
#define TINKER_GRO_WRITER_HPP

#include <string>
#include <memory>
#include <fstream>

#include "TrajectoryFormatWriter.hpp"


class Frame;


class GROWriter : public TrajectoryFormatWriter {

    std::fstream f;
public:

    void open(const std::string &filename) override;

    void close() override;

    void write(const std::shared_ptr<Frame> &frame) override;

};


#endif //TINKER_GRO_WRITER_HPP
