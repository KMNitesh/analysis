//
// Created by xiamr on 6/28/19.
//

#ifndef TINKER_GROMACSINTERFACE_HPP
#define TINKER_GROMACSINTERFACE_HPP

namespace gmx {

#include "gromacs/fileio/xtcio.h"
#include "gromacs/fileio/trnio.h"
#include "gromacs/utility/smalloc.h"

}

class GromacsInterface {
public:
    /* Open a trr / trr file */
    virtual gmx::t_fileio *open_trn(const char *fn, const char *mode) = 0;

    /* Close it */
    virtual void close_trn(gmx::t_fileio *fio) = 0;


    /* Write a trn frame to file fp, box, x, v, f may be NULL */
    virtual void fwrite_trn(gmx::t_fileio *fio, int step, gmx::real t, gmx::real lambda,
                            gmx::rvec *box, int natoms, gmx::rvec *x, gmx::rvec *v, gmx::rvec *f) = 0;

    /* Open a file for xdr I/O */
    virtual gmx::t_fileio *open_xtc(const char *filename, const char *mode) = 0;

    /* Close the file for xdr I/O */
    virtual void close_xtc(gmx::t_fileio *fio) = 0;

    /* Write a frame to xtc file */
    virtual int write_xtc(gmx::t_fileio *fio,
                          int natoms, int step, gmx::real time,
                          gmx::matrix box, gmx::rvec *x, gmx::real prec) = 0;


    virtual ~GromacsInterface() = default;
};


#endif //TINKER_GROMACSINTERFACE_HPP
