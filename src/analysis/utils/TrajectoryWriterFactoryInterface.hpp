//
// Created by xiamr on 6/28/19.
//

#ifndef TINKER_TRAJECTORYWRITERFACTORYINTERFACE_HPP
#define TINKER_TRAJECTORYWRITERFACTORYINTERFACE_HPP

#include <memory>
#include "utils/common.hpp"
#include "TrajectoryFormatWriter.hpp"

class TrajectoryWriterFactoryInterface {
public:
    virtual std::shared_ptr<TrajectoryFormatWriter> make_instance(FileType t) = 0;
};


#endif //TINKER_TRAJECTORYWRITERFACTORYINTERFACE_HPP
