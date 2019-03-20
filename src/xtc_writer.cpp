//
// Created by xiamr on 3/19/19.
//

#include "xtc_writer.hpp"
#include "frame.hpp"
#include "atom.hpp"
#include "gmxtrr.h"

void XTCWriter::open(const std::string &filename) {
    // First Time
    xd = gmx::open_xtc(filename.c_str(), "w");
    step = 0;
    time = 0.0;
    prec = 1000.0;

}

void XTCWriter::write(std::shared_ptr<Frame> &frame) {
    if (!x) {
        x = new gmx::rvec[frame->atom_list.size()];
    }
    int i = 0;
    for (auto &atom : frame->atom_list) {
        x[i][0] = atom->x / 10.0;
        x[i][1] = atom->y / 10.0;
        x[i][2] = atom->z / 10.0;
        i++;
    }

    gmx::rvec box[3];
    translate(frame->a_axis / 10.0, frame->b_axis / 10.0, frame->c_axis / 10.0,
              frame->alpha, frame->beta, frame->gamma, box);
    gmx::write_xtc(xd, i, step, time, box, x, prec);

    step++;
    time++;
}

void XTCWriter::close() {
    gmx::close_xtc(xd);
    delete[] x;
}
