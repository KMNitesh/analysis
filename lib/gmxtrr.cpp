#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>

#include "gmxtrr.h"
#include "../src/common.hpp"
using namespace std;



extern "C" {

void trr_close_(gmx::t_fileio **handle){
    gmx::t_fileio *fio = *handle;
    gmx::close_trn(fio);
}
void trr_open_(const char *fn, int *leng, gmx::t_fileio **handle){
    char filename[256];
    strncpy(filename,fn,*leng);
    filename[*leng]='\0';
    gmx::t_fileio *fio = gmx::open_trn(filename,"r");
    *handle = fio;
}

void trr_read_next_(gmx::t_fileio **handle,
                    double *xbox, double *ybox, double *zbox, double *alpha, double *beta, double *gamma,
                    double *x, double *y, double *z,int *ret){

    gmx::t_fileio *fio = *handle;
    gmx::t_trnheader trnheader;
    gmx::gmx_bool bOK;
    gmx::rvec box[3];
    gmx::rvec *coord = nullptr;
    if (gmx::fread_trnheader(fio,&trnheader,&bOK)){
        if (bOK){
            coord = new gmx::rvec[trnheader.x_size];
            if (trnheader.box_size){
                gmx::fread_htrn(fio,&trnheader,box,coord,NULL,NULL);
                translate(box,xbox,ybox,zbox,alpha,beta,gamma);
            }else{
                *xbox = 0.0;
                *ybox = 0.0;
                *zbox = 0.0;
                *alpha = 0.0;
                *beta = 0.0;
                *gamma = 0.0;
                gmx::fread_htrn(fio,&trnheader,NULL,coord,NULL,NULL);
            }
            for( int i = 0; i < trnheader.x_size; i++){
                x[i] = coord[i][0] *10;
                y[i] = coord[i][1] *10;
                z[i] = coord[i][2] *10;
            }
            *ret = 1;
        } else{
            *ret = 0;
        }
    } else{ *ret = 0;}

    delete [] coord;
}

void append_frame_x_box_(const char *fn,int *leng, int *step, double *time,
                         double *xbox, double *ybox,double *zbox,double *alpha, double *beta, double *gamma,
                         int *natoms, double *x, double *y, double *z ){
    char filename[256];
    strncpy(filename,fn,*leng);
    filename[*leng] = '\0';
    fstream trrfile(filename);
    gmx::t_fileio *fio;
    if (trrfile.good()){
        fio = gmx::open_trn(filename, "a");
    }else{
        fio = gmx::open_trn(filename, "w");
    }
    gmx::rvec box[3];
    translate(*xbox/10.0,*ybox/10.0,*zbox/10.0,*alpha,*beta,*gamma, box);

    auto coord = new gmx::rvec[(*natoms)];
    for (int i = 0; i < *natoms; i++){
        coord[i][0] = x[i] /10.0;
        coord[i][1] = y[i] /10.0;
        coord[i][2] = z[i] /10.0;
    }
    gmx::fwrite_trn(fio, *step, *time, 0.0, box, *natoms , coord, NULL,NULL);
    gmx::close_trn(fio);
    delete [] coord;
}
void append_frame_x_(const char *fn,int *leng, int *step, double *time,
                     int *natoms, double *x, double *y, double *z ){
    char filename[256];
    strncpy(filename,fn,*leng);
    filename[*leng] = '\0';
    fstream trrfile(filename);
    gmx::t_fileio *fio;
    if (trrfile.good()){
        fio = gmx::open_trn(filename, "a");
    }else{
        fio = gmx::open_trn(filename, "w");
    }

    auto coord = new gmx::rvec[(*natoms)];
    for (int i = 0; i < *natoms; i++){
        coord[i][0] = x[i] /10.0;
        coord[i][1] = y[i] /10.0;
        coord[i][2] = z[i] /10.0;
    }
    gmx::fwrite_trn(fio, *step, *time, 0.0, NULL, *natoms , coord, NULL,NULL);
    gmx::close_trn(fio);
    delete [] coord;
}
}


void translate( gmx::real xbox, gmx::real ybox, gmx::real zbox,
                gmx::real alpha, gmx::real beta, gmx::real gamma, gmx::rvec *box){
    bool triclinic = false;
    bool orthogonal = false;
    bool octahedron = false;
    bool monoclinic = false;
    if (alpha == 90.0 and  beta == 90.0 and gamma == 90.0)
        orthogonal = true;
    else if (alpha == 90.0 and gamma == 90.0)
        monoclinic = true;
    else
        triclinic = true;


    if (octahedron){
        if (xbox == ybox and xbox == zbox and orthogonal){
            orthogonal = false;
            monoclinic = false;
            triclinic = false;
        }
    }
    double alpha_cos,beta_sin,beta_cos,gamma_sin,gamma_cos,beta_term,gamma_term;
    if (orthogonal or octahedron) {
        alpha_cos = 0.0;
        beta_sin = 1.0;
        beta_cos = 0.0;
        gamma_sin = 1.0;
        gamma_cos = 0.0;
        beta_term = 0.0;
        gamma_term = 1.0;
    } else if (monoclinic) {
        alpha_cos = 0.0;
        beta_sin = sin(beta/radian);
        beta_cos = cos(beta/radian);
        gamma_sin = 1.0;
        gamma_cos = 0.0;
        beta_term = 0.0;
        gamma_term = beta_sin;
    } else if (triclinic) {
        alpha_cos = cos(alpha/radian);
        beta_sin = sin(beta/radian);
        beta_cos = cos(beta/radian);
        gamma_sin = sin(gamma/radian);
        gamma_cos = cos(gamma/radian);
        beta_term = (alpha_cos - beta_cos*gamma_cos) / gamma_sin;
        gamma_term = sqrt(beta_sin*beta_sin - beta_term*beta_term);
    }
    box[0][0] = xbox;
    box[0][1] = 0.0;
    box[0][2] = 0.0;
    box[1][0] = ybox * gamma_cos;
    box[1][1] = ybox * gamma_sin;
    box[1][2] = 0.0;
    box[2][0] = zbox * beta_cos;
    box[2][1] = zbox * beta_term;
    box[2][2] = zbox * gamma_term;
}

void translate(gmx::rvec *box, double *xbox, double *ybox, double *zbox, double *alpha, double *beta, double *gamma) {

    if (box[0][1] != 0.0 or box[0][2] != 0.0 or box[1][0] != 0.0 or box[1][2] != 0.0 or box[2][0] != 0.0 or box[2][1] != 0.0){
        std::cerr << "box is not orthogonal !" << endl;
        exit(1);
    }
    // bugbug only work in orthogonal case
    *xbox = box[0][0] *10;
    *ybox = sqrt(box[1][0]*box[1][0]+box[1][1]*box[1][1]) *10;
//    *gamma = radian * acos(box[1][0]/(*ybox));
     *zbox = box[2][2] *10;
     *alpha = 90.00;
     *beta  = 90.00;
     *gamma = 90.00;
}



