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

#include "boost/optional.hpp"
#include "AmberNetcdf.h"

#include "TrajectoryFormatWriter.hpp"
#include "NetcdfInterface.hpp"
#include "NetcdfImpl.hpp"

class Frame;


class NetCDFWriter : public TrajectoryFormatWriter {
    struct AmberNetcdf NC;
    bool _is_open = false;
    int step = 0;
    std::string filename;
public:
    void open(const std::string &filename) override;

    void close() override;

    void write(const std::shared_ptr<Frame> &frame) override;

protected:
    virtual NetcdfInterface *getNetcdfImpl() {
        static NetcdfImpl impl;
        return &impl;
    }
};


#endif //TINKER_NETCDF_WRITER_HPP
