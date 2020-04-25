#ifndef TINKER_NETCDFTRAJECTORYREADER_HPP
#define TINKER_NETCDFTRAJECTORYREADER_HPP

#include "TrajectoryInterface.hpp"
#include "utils/AmberNetcdf.h"

class NetcdfTrajectoryReader : public TrajectoryInterface {
public:
    bool open(const std::string &file) override;

    void close() override;

    ~NetcdfTrajectoryReader() override;

protected:
    bool readOneFrameImpl(std::shared_ptr<Frame> &frame, const std::vector<std::shared_ptr<Atom>> &atoms) override;

private:
    std::unique_ptr<struct AmberNetcdf> NC;
};

#endif // TINKER_NETCDFTRAJECTORYREADER_HPP
