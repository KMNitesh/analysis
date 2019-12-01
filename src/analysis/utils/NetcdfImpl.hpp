//
// Created by xiamr on 6/28/19.
//

#ifndef TINKER_NETCDFIMPL_HPP
#define TINKER_NETCDFIMPL_HPP

#include "NetcdfInterface.hpp"

class NetcdfImpl : public NetcdfInterface {
    int netcdfClose(struct AmberNetcdf *A) override {
        return ::netcdfClose(A);
    }

    int netcdfCreate(struct AmberNetcdf *A, const char *filename, int natom, int isBox) override {
        return ::netcdfCreate(A, filename, natom, isBox);
    }

    int netcdfWriteNextFrame(struct AmberNetcdf *A, double *X, double *box) override {
        return ::netcdfWriteNextFrame(A, X, box);
    }

};


#endif //TINKER_NETCDFIMPL_HPP
