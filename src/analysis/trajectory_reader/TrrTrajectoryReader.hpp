#ifndef TINKER_TRRTRAJECTORYREADER_HPP
#define TINKER_TRRTRAJECTORYREADER_HPP

#include "TrajectoryInterface.hpp"

namespace gmx {

#include "gromacs/fileio/xtcio.h"
#include "gromacs/utility/smalloc.h"

}  // namespace gmx

class TrrTrajectoryReader : public TrajectoryInterface {
public:
    bool open(const std::string &file) override;

    void close() override;

    ~TrrTrajectoryReader() override;

protected:
    bool readOneFrameImpl(std::shared_ptr<Frame> &frame, const std::vector<std::shared_ptr<Atom>> &atoms) override;

private:
    gmx::t_fileio *fio = nullptr;
};

#endif  // TINKER_TRRTRAJECTORYREADER_HPP
