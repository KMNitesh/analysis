//
// Created by xiamr on 4/13/18.
//

#ifndef TINKER_AMBER_NETCDF_H
#define TINKER_AMBER_NETCDF_H

extern "C" {


void netcdf_read_next_(struct AmberNetcdf **handle,
                       double *xbox, double *ybox, double *zbox, double *alpha, double *beta, double *gamma,
                       double *x, double *y, double *z, int *ret);
};


#endif //TINKER_AMBER_NETCDF_H
