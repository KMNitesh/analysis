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

    void close() override {};

protected:
    bool readOneFrameImpl(std::shared_ptr<Frame> &frame) override;

private:

    static void apply_box(std::shared_ptr<Frame> &frame);

    void parse_box(std::shared_ptr<Frame> &frame);

    void parse_coord(std::shared_ptr<Atom> &atom);

    void readOneFrameVelocity(std::shared_ptr<Frame> &frame);

    std::ifstream ifs, velocity_file;

    std::string line;
    std::vector<std::string> field;


};


#endif //TINKER_ARCTRAJECTORYREADER_HPP
