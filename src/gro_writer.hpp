//
// Created by xiamr on 3/19/19.
//

#ifndef TINKER_GRO_WRITER_HPP
#define TINKER_GRO_WRITER_HPP

#include <string>
#include <memory>

class Frame;


class GROWriter {
public:
    void write(const std::string &filename, std::shared_ptr<Frame> &frame);
};


#endif //TINKER_GRO_WRITER_HPP
