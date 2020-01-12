
#include "utils/common.hpp"
#include "data_structure/atom.hpp"
#include "data_structure/frame.hpp"
#include "XtcTrajectoryReader.hpp"
#include "utils/gmxtrr.h"

bool XtcTrajectoryReader::open(const std::string &file) {
    fio = gmx::open_xtc(file.c_str(), "r");
    return fio != nullptr;
}

bool XtcTrajectoryReader::readOneFrameImpl(std::shared_ptr<Frame> &frame) {
    gmx::matrix box;
    gmx::gmx_bool bOK;

    int ret;
    if (!x) {
        ret = gmx::read_first_xtc(fio, &natoms, &step, &time, box, &x, &prec, &bOK);
    } else {
        ret = gmx::read_next_xtc(fio, natoms, &step, &time, box, x, &prec, &bOK);
    }
    if (!ret) return false;

    if (bOK) {
        if (natoms != static_cast<int>(frame->atom_list.size())) {
            std::cerr << "ERROR! the atom number do not match" << std::endl;
            exit(1);
        }
        if (frame->enable_bound) {
            translate(box, &(frame->a_axis), &(frame->b_axis), &(frame->c_axis),
                      &(frame->alpha), &(frame->beta), &(frame->gamma));
        }
        int i = 0;
        for (auto &atom : frame->atom_list) {
            atom->x = x[i][0] * 10;
            atom->y = x[i][1] * 10;
            atom->z = x[i][2] * 10;
            i++;
        }
        return true;
    }
    std::cerr << "\nWARNING: Incomplete frame at time "
              << std::scientific << time << std::defaultfloat << "\n";
    return false;
}

void XtcTrajectoryReader::close() {
    if (x) {
        gmx::sfree(x);
        x = nullptr;
    }
    if (fio) {
        gmx::close_xtc(fio);
        fio = nullptr;
    }
}

XtcTrajectoryReader::~XtcTrajectoryReader() {
    XtcTrajectoryReader::close();
}
