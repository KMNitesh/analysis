#ifndef TINKER_XTCTRAJECTORYREADER_HPP
#define TINKER_XTCTRAJECTORYREADER_HPP

#include "TrajectoryInterface.hpp"

namespace gmx {

#include "gromacs/fileio/xtcio.h"
#include "gromacs/utility/smalloc.h"

}


class XtcTrajectoryReader : public TrajectoryInterface {
public:
    bool open(const std::string &file) override;

    void close() override;

    ~XtcTrajectoryReader() override;

protected:
    bool readOneFrameImpl(std::shared_ptr<Frame> &frame) override;

private:

    gmx::t_fileio *fio = nullptr;

    int natoms, step;
    gmx::rvec *x = nullptr;
    gmx::real prec, time;

};


#endif //TINKER_XTCTRAJECTORYREADER_HPP
