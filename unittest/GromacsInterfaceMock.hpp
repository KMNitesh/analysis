//
// Created by xiamr on 6/28/19.
//

#ifndef TINKER_GROMACSINTERFACEMOCK_HPP
#define TINKER_GROMACSINTERFACEMOCK_HPP

#include <gmock/gmock.h>
#include "GromacsInterface.hpp"

using namespace testing;

class GromacsMock : public GromacsInterface {
public:
    MOCK_METHOD2(open_trn, gmx::t_fileio *(
            const char *fn,
            const char *mode));

    MOCK_METHOD1(close_trn, void(gmx::t_fileio
            *fio));

    MOCK_METHOD9(fwrite_trn, void(gmx::t_fileio
            *fio, int
            step, gmx::real
            t, gmx::real
            lambda, gmx::rvec * box, int
            natoms,
                    gmx::rvec * x, gmx::rvec * v, gmx::rvec * f));

    MOCK_METHOD2(open_xtc, gmx::t_fileio *(
            const char *filename,
            const char *mode));

    MOCK_METHOD1(close_xtc, void(gmx::t_fileio
            *fio));

    MOCK_METHOD7(
            write_xtc, int(gmx::t_fileio
            *fio, int
            natoms, int
            step, gmx::real
            time, gmx::matrix
            box, gmx::rvec * x, gmx::real
            prec));
};


#endif //TINKER_GROMACSINTERFACEMOCK_HPP
