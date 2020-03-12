#ifndef TINKER_TPRREADER_HPP
#define TINKER_TPRREADER_HPP

#include "TopologyInterface.hpp"

namespace gmx {

#include "gromacs/fileio/trnio.h"
#include "gromacs/fileio/xtcio.h"
#include "gromacs/utility/smalloc.h"

}  // namespace gmx

class TprReader : public TopologyInterface {
public:
    std::shared_ptr<Frame> read(const std::string &filename) override;

private:
    gmx::t_fileio *fio = nullptr;
    int natoms, step;
    gmx::rvec *x = nullptr;
    gmx::real prec, time;
};

#endif  // TINKER_TPRREADER_HPP
