#ifndef TINKER_MOL2READER_HPP
#define TINKER_MOL2READER_HPP

#include "TopologyInterface.hpp"

class Mol2Reader : public TopologyInterface {
public:
    std::shared_ptr<Frame> read(const std::string &filename) override;
};


#endif //TINKER_MOL2READER_HPP
