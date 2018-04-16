//
// Created by xiamr on 4/13/18.
//

#ifndef TINKER_GMXTRR_H
#define TINKER_GMXTRR_H

namespace gmx{
#include "gromacs/fileio/trnio.h"
}
const double radian = 57.29577951308232088;

void translate( gmx::real xbox, gmx::real ybox, gmx::real zbox,
                gmx::real alpha, gmx::real beta, gmx::real gamma, gmx::rvec *box);

void translate( gmx::rvec *box, double *xbox, double *ybox, double *zbox,
                double *alpha, double *beta, double *gamma);

#endif //TINKER_GMXTRR_H
