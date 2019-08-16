//
// Created by xiamr on 8/16/19.
//

#include "ConvertVelocityToVelocityCharge.hpp"
#include "common.hpp"
#include "frame.hpp"
#include "trr_writer.hpp"
#include "atom.hpp"

ConvertVelocityToVelocityCharge::ConvertVelocityToVelocityCharge() {
    enable_read_velocity = true;
}

void ConvertVelocityToVelocityCharge::processFirstFrame(std::shared_ptr<Frame> &frame) {
    writer = std::make_unique<TRRWriter>();
    writer->open(trr_vq_outfilename);
    writer->setWriteVelocities(true);
}

void ConvertVelocityToVelocityCharge::process(std::shared_ptr<Frame> &frame) {
    std::tuple<double, double, double> vq;
    for (auto &atom : frame->atom_list) {
        vq += atom->charge.value() * atom->getVelocities();
        atom->vx = atom->vy = atom->vz = 0.0;
    }
    auto &last_atom = frame->atom_list.back();
    last_atom->vx = std::get<0>(vq);
    last_atom->vy = std::get<1>(vq);
    last_atom->vz = std::get<2>(vq);

    writer->setCurrentTime(frame->getCurrentTime().value());
    writer->write(frame);
}

void ConvertVelocityToVelocityCharge::print(std::ostream &os) {
    writer->close();
}

void ConvertVelocityToVelocityCharge::readInfo() {
    trr_vq_outfilename = choose_file("Enter trr output filename > ", false, "trr");
}


