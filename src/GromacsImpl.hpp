//
// Created by xiamr on 6/28/19.
//

#ifndef TINKER_GROMACSIMPL_HPP
#define TINKER_GROMACSIMPL_HPP

#include "GromacsInterface.hpp"

class GromacsImpl : public GromacsInterface {
public:
    gmx::t_fileio *open_trn(const char *fn, const char *mode) override { return gmx::open_trn(fn, mode); }

    void close_trn(gmx::t_fileio *fio) override { gmx::close_trn(fio); }

    void
    fwrite_trn(gmx::t_fileio *fio, int step, gmx::real t, gmx::real lambda, gmx::rvec *box, int natoms, gmx::rvec *x,
               gmx::rvec *v, gmx::rvec *f) override {
        gmx::fwrite_trn(fio, step, t, lambda, box, natoms, x, v, f);
    }

    gmx::t_fileio *open_xtc(const char *filename, const char *mode) override { return gmx::open_xtc(filename, mode); }

    void close_xtc(gmx::t_fileio *fio) override { gmx::close_xtc(fio); }

    int write_xtc(gmx::t_fileio *fio, int natoms, int step, gmx::real time, gmx::matrix box, gmx::rvec *x,
                  gmx::real prec) override {
        return gmx::write_xtc(fio, natoms, step, time, box, x, prec);
    }
};

#endif // TINKER_GROMACSIMPL_HPP
