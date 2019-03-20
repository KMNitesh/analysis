//
// Created by xiamr on 3/19/19.
//

#ifndef TINKER_NETCDF_WRITER_HPP
#define TINKER_NETCDF_WRITER_HPP

#include <string>
#include <memory>
namespace gmx {

#include "gromacs/fileio/xtcio.h"
#include "gromacs/fileio/trnio.h"
#include "gromacs/utility/smalloc.h"

}

#include "AmberNetcdf.h"

class Frame;


class NetCDFWriter {
    struct AmberNetcdf NC;
    double *x = nullptr;
    int step = 0;
public:
    void open(const std::string &filename, int natom);

    void close();

    void write(const std::shared_ptr<Frame> &frame);
};


#endif //TINKER_NETCDF_WRITER_HPP
