#ifndef TINKER_TRRTRAJECTORYREADER_HPP
#define TINKER_TRRTRAJECTORYREADER_HPP

#include "TrajectoryInterface.hpp"


namespace gmx {

#include "gromacs/fileio/xtcio.h"
#include "gromacs/utility/smalloc.h"

}

class TrrTrajectoryReader : public TrajectoryInterface {
public:
    bool open(const std::string &file) override;

    bool readOneFrame(std::shared_ptr<Frame> &frame) override;

    void close() override;

    ~TrrTrajectoryReader() override;

private:
    gmx::t_fileio *fio = nullptr;
};


#endif //TINKER_TRRTRAJECTORYREADER_HPP
