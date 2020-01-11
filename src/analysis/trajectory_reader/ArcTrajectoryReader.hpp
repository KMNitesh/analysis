#ifndef TINKER_ARCTRAJECTORYREADER_HPP
#define TINKER_ARCTRAJECTORYREADER_HPP

#include <fstream>
#include "TrajectoryInterface.hpp"
#include "TopologyInterface.hpp"

class Frame;

class ArcTrajectoryReader : public TrajectoryInterface, public TopologyInterface {
public:
    bool open(const std::string &file) override;

    std::shared_ptr<Frame> read(const std::string &filename) override;

    bool readOneFrame(std::shared_ptr<Frame> &frame) override;

    void close() override {};

private:

    std::ifstream ifs;

    std::string line;
    std::vector<std::string> field;
};


#endif //TINKER_ARCTRAJECTORYREADER_HPP
