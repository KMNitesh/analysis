//
// Created by xiamr on 3/19/19.
//

#include <boost/range/adaptors.hpp>

#include "trr_writer.hpp"
#include "data_structure/frame.hpp"
#include "data_structure/atom.hpp"
#include "ThrowAssert.hpp"

void TRRWriter::open(const std::string &filename) {
    // First Time
    throw_assert(xd == nullptr, "Gromacs Trr File handle not close");
    xd = getGromacsImpl()->open_trn(filename.c_str(), "w");
    throw_assert(xd != nullptr, "Gromacs Trr File handle open error");
    step = 0;
    time = 0.0;
}

void TRRWriter::write(const std::shared_ptr<Frame> &frame) {
    throw_assert(xd != nullptr, "Gromacs Trr File handle invalid");
    auto x = std::make_unique<gmx::rvec[]>(frame->atom_list.size());
    for (const auto &ele : frame->atom_list | boost::adaptors::indexed()) {
        auto i = ele.index();
        auto &atom = ele.value();
        x[i][0] = atom->x / 10.0;
        x[i][1] = atom->y / 10.0;
        x[i][2] = atom->z / 10.0;
    }
    std::unique_ptr<gmx::rvec[]> v;
    if (writeVelocities) {
        v = std::make_unique<gmx::rvec[]>(frame->atom_list.size());
        for (const auto &ele : frame->atom_list | boost::adaptors::indexed()) {
            auto i = ele.index();
            auto &atom = ele.value();
            v[i][0] = atom->vx / 10.0;
            v[i][1] = atom->vy / 10.0;
            v[i][2] = atom->vz / 10.0;
        }
    }

    gmx::matrix box;
    frame->box.getBoxParameter(box);
    getGromacsImpl()->fwrite_trn(xd, step, time, 0.0, box, frame->atom_list.size(), x.get(), v.get(), NULL);
    step++;
    time++;
}

void TRRWriter::close() {
    throw_assert(xd != nullptr, "Gromacs Trr File handle invalid");
    getGromacsImpl()->close_trn(xd);
    xd = nullptr;
}
