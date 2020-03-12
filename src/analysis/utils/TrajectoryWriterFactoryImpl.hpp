//
// Created by xiamr on 6/28/19.
//

#ifndef TINKER_TRAJECTORYWRITERFACTORYIMPL_HPP
#define TINKER_TRAJECTORYWRITERFACTORYIMPL_HPP

#include <memory>

#include "TrajectoryWriterFactoryInterface.hpp"
#include "utils/common.hpp"

class TrajectoryWriterFactoryImpl : public TrajectoryWriterFactoryInterface {
public:
    std::shared_ptr<TrajectoryFormatWriter> make_instance(FileType t) override;
};

#endif  // TINKER_TRAJECTORYWRITERFACTORYIMPL_HPP
