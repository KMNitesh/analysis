#ifndef TINKER_NETCDFTRAJECTORYREADER_HPP
#define TINKER_NETCDFTRAJECTORYREADER_HPP

#include "TrajectoryInterface.hpp"

#include "utils/AmberNetcdf.h"
#include "utils/amber_netcdf.h"

class NetcdfTrajectoryReader : public TrajectoryInterface {
public:
    bool open(const std::string &file) override;

    bool readOneFrame(std::shared_ptr<Frame> &frame) override;

    void close() override;

    ~NetcdfTrajectoryReader() override;

private:
    std::unique_ptr<struct AmberNetcdf> NC;
};


#endif //TINKER_NETCDFTRAJECTORYREADER_HPP
