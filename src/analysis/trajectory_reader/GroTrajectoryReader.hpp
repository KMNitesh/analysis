#ifndef TINKER_GROTRAJECTORYREADER_HPP
#define TINKER_GROTRAJECTORYREADER_HPP

#include <fstream>

#include "TrajectoryInterface.hpp"

class GroTrajectoryReader : public TrajectoryInterface {
   public:
    bool open(const std::string &file) override;

    void close() override;

   protected:
    bool readOneFrameImpl(std::shared_ptr<Frame> &frame) override;

    std::ifstream ifs;
};

#endif  // TINKER_GROTRAJECTORYREADER_HPP
