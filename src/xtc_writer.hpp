//
// Created by xiamr on 3/19/19.
//

#ifndef TINKER_XTC_WRITER_HPP
#define TINKER_XTC_WRITER_HPP

#include <string>
#include <memory>
namespace gmx {

#include "gromacs/fileio/xtcio.h"
#include "gromacs/fileio/trnio.h"
#include "gromacs/utility/smalloc.h"

}

class Frame;

class XTCWriter {
    gmx::t_fileio *xd = nullptr;
    gmx::rvec *x = nullptr;
    int step;
    gmx::real prec, time;
public:
    void open(const std::string &filename);

    void write(std::shared_ptr<Frame> &frame);

    void close();

};



#endif //TINKER_XTC_WRITER_HPP
