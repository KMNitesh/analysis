
#include "xtc_writer.hpp"
#include "data_structure/frame.hpp"
#include "data_structure/atom.hpp"
#include "ThrowAssert.hpp"

void XTCWriter::open(const std::string &filename) {
    // First Time
    throw_assert(xd == nullptr, "Gromacs Xtc File handle not close");
    xd = getGromacsImpl()->open_xtc(filename.c_str(), "w");
    throw_assert(xd != nullptr, "Gromacs Xtc File handle open error");
    step = 0;
    time = 0.0;
    prec = 1000.0;

}

void XTCWriter::write(const std::shared_ptr<Frame> &frame) {
    throw_assert(xd != nullptr, "Gromacs Xtc File handle invalid");
    auto x = std::make_unique<gmx::rvec[]>(frame->atom_list.size());
    int i = 0;
    for (auto &atom : frame->atom_list) {
        x[i][0] = atom->x / 10.0;
        x[i][1] = atom->y / 10.0;
        x[i][2] = atom->z / 10.0;
        i++;
    }

    gmx::matrix box;
    frame->box.getBoxParameter(box);
    getGromacsImpl()->write_xtc(xd, i, step, time, box, x.get(), prec);

    step++;
    time++;
}

void XTCWriter::close() {
    throw_assert(xd != nullptr, "Gromacs Xtc File handle invalid");
    getGromacsImpl()->close_xtc(xd);
    xd = nullptr;
}
