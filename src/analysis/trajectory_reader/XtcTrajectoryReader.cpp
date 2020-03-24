
#include "XtcTrajectoryReader.hpp"

#include "data_structure/atom.hpp"
#include "data_structure/frame.hpp"
#include "utils/common.hpp"

bool XtcTrajectoryReader::open(const std::string &file) {
    fio = gmx::open_xtc(file.c_str(), "r");
    return fio != nullptr;
}

bool XtcTrajectoryReader::readOneFrameImpl(std::shared_ptr<Frame> &frame,
                                           const std::vector<std::shared_ptr<Atom>> &atoms) {
    gmx::matrix box;
    gmx::gmx_bool bOK;

    int ret;
    if (!x) {
        ret = gmx::read_first_xtc(fio, &natoms, &step, &time, box, &x, &prec, &bOK);
    } else {
        ret = gmx::read_next_xtc(fio, natoms, &step, &time, box, x, &prec, &bOK);
    }
    if (!ret)
        return false;

    if (bOK) {
        static bool has_Warning_d = false;
        if (!has_Warning_d and natoms != static_cast<int>(atoms.size())) {
            std::cerr << boost::format("WARNING: topology has %d atoms, whereas trajectory has %d\n") % atoms.size() %
                             natoms;
            has_Warning_d = true;
        }
        if (frame->enable_bound) {
            frame->box = PBCBox(box);
        }
        for (auto i : boost::irange(natoms)) {
            auto &atom = atoms[i];
            atom->x = x[i][0] * 10;
            atom->y = x[i][1] * 10;
            atom->z = x[i][2] * 10;
        }
        return true;
    }
    std::cerr << "\nWARNING: Incomplete frame at time " << std::scientific << time << std::defaultfloat << "\n";
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

XtcTrajectoryReader::~XtcTrajectoryReader() { XtcTrajectoryReader::close(); }
