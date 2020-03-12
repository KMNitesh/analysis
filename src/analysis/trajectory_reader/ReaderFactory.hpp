#ifndef TINKER_READERFACTORY_HPP
#define TINKER_READERFACTORY_HPP

#include "TopologyInterface.hpp"
#include "TrajectoryInterface.hpp"

class ReaderFactory {
public:
    [[nodiscard]] static std::shared_ptr<TopologyInterface> getTopology(const std::string &filename);

    [[nodiscard]] static std::shared_ptr<TrajectoryInterface> getTrajectory(const std::string &filename);
};

#endif  // TINKER_READERFACTORY_HPP
