//
// Created by xiamr on 6/28/19.
//

#ifndef TINKER_NETCDFINTERFACE_HPP
#define TINKER_NETCDFINTERFACE_HPP

#include "AmberNetcdf.h"

class NetcdfInterface {
public:
    virtual int netcdfClose(struct AmberNetcdf *A) = 0;

    virtual int netcdfCreate(struct AmberNetcdf *A, const char *filename, int natom, int isBox) = 0;

    virtual int netcdfWriteNextFrame(struct AmberNetcdf *A, double *X, double *box) = 0;

    virtual ~NetcdfInterface() = default;
};

#endif //TINKER_NETCDFINTERFACE_HPP
