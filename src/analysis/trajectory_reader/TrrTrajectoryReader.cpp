
#include "utils/common.hpp"
#include "data_structure/atom.hpp"
#include "data_structure/frame.hpp"
#include "TrrTrajectoryReader.hpp"
#include "utils/gmxtrr.h"

bool TrrTrajectoryReader::open(const std::string &file) {
    fio = gmx::open_trn(file.c_str(), "r");
    return fio != nullptr;
}

bool TrrTrajectoryReader::readOneFrameImpl(std::shared_ptr<Frame> &frame) {
    gmx::t_trnheader trnheader;
    gmx::gmx_bool bOK;
    gmx::rvec box[3];
    std::unique_ptr<gmx::rvec[]> coord;
    std::unique_ptr<gmx::rvec[]> velocities;
    if (gmx::fread_trnheader(fio, &trnheader, &bOK)) {
        frame->setCurrentTime(trnheader.t);
        if (bOK) {
            coord = std::make_unique<gmx::rvec[]>(trnheader.natoms);
            if (trnheader.v_size) {
                velocities = std::make_unique<gmx::rvec[]>(trnheader.natoms);
                frame->has_velocity = true;
            } else if (enable_read_velocity) {
                std::cerr << "ERROR!! Gromacs TRR Trajectory file does not have velocities !\n";
                exit(4);
            }
            if (trnheader.box_size) {
                gmx::fread_htrn(fio, &trnheader, box, coord.get(), velocities.get(), nullptr);
                translate(box, &(frame->a_axis), &(frame->b_axis), &(frame->c_axis),
                          &(frame->alpha), &(frame->beta), &(frame->gamma));
                if (frame->enable_bound) {
                    frame->a_axis_half = frame->a_axis / 2;
                    frame->b_axis_half = frame->b_axis / 2;
                    frame->c_axis_half = frame->c_axis / 2;
                }
            } else {
                if (frame->enable_bound) {
                    std::cerr << "ERROR !! trr trajectory does not  have PBC enabled" << std::endl;
                    exit(1);
                }
                frame->a_axis = 0.0;
                frame->b_axis = 0.0;
                frame->c_axis = 0.0;
                frame->alpha = 0.0;
                frame->beta = 0.0;
                frame->gamma = 0.0;
                frame->enable_bound = false;
                gmx::fread_htrn(fio, &trnheader, nullptr, coord.get(), velocities.get(), nullptr);
            }
            if (static_cast<int>(frame->atom_list.size()) != trnheader.natoms) {
                std::cerr << "ERROR! the atom number do not match" << std::endl;
                exit(1);
            }
            int i = 0;
            for (auto &atom : frame->atom_list) {
                atom->x = coord[i][0] * 10;
                atom->y = coord[i][1] * 10;
                atom->z = coord[i][2] * 10;
                if (velocities) {
                    atom->vx = velocities[i][0] * 10;
                    atom->vy = velocities[i][1] * 10;
                    atom->vz = velocities[i][2] * 10;
                }
                i++;
            }

            return true;
        }
    }
    return false;

}

void TrrTrajectoryReader::close() {
    if (fio) {
        gmx::close_trn(fio);
        fio = nullptr;
    }
}


TrrTrajectoryReader::~TrrTrajectoryReader() {
    TrrTrajectoryReader::close();
}
