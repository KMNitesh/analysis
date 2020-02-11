#include <fstream>
#include <cstring>
#include <cmath>
#include <iostream>

#include "AmberNetcdf.h"

using namespace std;

extern "C" {

void netcdf_rewind_(struct AmberNetcdf **handle) {
    auto nc = *handle;
    nc->currentFrame = 0;
}
void netcdf_close_(struct AmberNetcdf **handle) {
    if (netcdfClose(*handle)) {
        cerr << "Error close amber netcdf file" << endl;
        exit(1);
    }
    delete *handle;
}

void netcdf_open_(const char *fn, int *leng, struct AmberNetcdf **handle, int *ret) {
    char filename[256];
    strncpy(filename, fn, *leng);
    filename[*leng] = '\0';
    auto nc = new struct AmberNetcdf();
    *ret = netcdfLoad(nc, filename);
    *handle = nc;
}

void netcdf_read_next_(struct AmberNetcdf **handle,
                       double *xbox, double *ybox, double *zbox, double *alpha, double *beta, double *gamma,
                       double *x, double *y, double *z, int *ret) {
    auto nc = *handle;
    auto coord = new double[nc->ncatom3];
    double box[6];
    *ret = netcdfGetNextFrame(nc, coord, box);
    if (*ret == 0) return;
    for (int i = 0; i < nc->ncatom; i++) {
        x[i] = coord[3 * i];
        y[i] = coord[3 * i + 1];
        z[i] = coord[3 * i + 2];
    }
    *xbox = box[0];
    *ybox = box[1];
    *zbox = box[2];
    *alpha = box[3];
    *beta = box[4];
    *gamma = box[5];

    delete[] coord;
}

void netcdf_append_box_x_(const char *fn, int *leng,
                          double *xbox, double *ybox, double *zbox, double *alpha, double *beta, double *gamma,
                          int *natoms, double *x, double *y, double *z) {
    char filename[256];
    strncpy(filename, fn, *leng);
    filename[*leng] = '\0';
    fstream ncfile(filename);
    struct AmberNetcdf nc;
    if (ncfile.good()) {
        if (netcdfOpen(&nc, filename)) {
            std::cerr << "Error open " << filename << endl;
            exit(1);
        };
    } else {
        if (netcdfCreate(&nc, filename, *natoms, 1)) {
            std::cerr << "Error create " << filename << endl;
            exit(1);
        }
    }
    double box[6];
    box[0] = *xbox;
    box[1] = *ybox;
    box[2] = *zbox;
    box[3] = *alpha;
    box[4] = *beta;
    box[5] = *gamma;

    auto coord = new double[nc.ncatom3];
    for (int i = 0; i < *natoms; i++) {
        coord[3 * i] = x[i];
        coord[3 * i + 1] = y[i];
        coord[3 * i + 2] = z[i];
    }

    int frame_set = nc.ncframe < 0 ? 0 : nc.ncframe;
    if (netcdfWriteFrame(&nc, frame_set, coord, box)) {
        cerr << "Error write  amber netcdf frame " << endl;
        exit(1);
    }

    if (netcdfClose(&nc)) {
        cerr << "Error close amber netcdf file" << endl;
        exit(1);
    }
    delete[] coord;

}
void netcdf_append_x_(const char *fn, int *leng, int *natoms, double *x, double *y, double *z) {
    char filename[256];
    strncpy(filename, fn, *leng);
    filename[*leng] = '\0';
    fstream ncfile(filename);
    struct AmberNetcdf nc;
    if (ncfile.good()) {
        if (netcdfOpen(&nc, filename)) {
            std::cerr << "Error open " << filename << endl;
            exit(1);
        }
    } else {
        if (netcdfCreate(&nc, filename, *natoms, 1)) {
            std::cerr << "Error create " << filename << endl;
            exit(1);
        }
    }

    auto coord = new double[nc.ncatom3];
    for (int i = 0; i < *natoms; i++) {
        coord[3 * i] = x[i];
        coord[3 * i + 1] = y[i];
        coord[3 * i + 2] = z[i];
    }
    int frame_set = nc.ncframe < 0 ? 0 : nc.ncframe;
    if (netcdfWriteFrame(&nc, frame_set, coord, NULL)) {
        cerr << "Error write  amber netcdf frame " << endl;
        exit(1);
    }

    if (netcdfClose(&nc)) {
        cerr << "Error close amber netcdf file" << endl;
        exit(1);
    }
    delete[] coord;
}
}
