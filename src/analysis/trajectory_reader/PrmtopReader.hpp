#ifndef TINKER_PRMTOPREADER_HPP
#define TINKER_PRMTOPREADER_HPP

#include "TopologyInterface.hpp"

class PrmtopReader : public TopologyInterface {
public:
    std::shared_ptr<Frame> read(const std::string &filename) override;
};

#endif  // TINKER_PRMTOPREADER_HPP
