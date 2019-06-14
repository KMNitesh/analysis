//
// Created by xiamr on 6/14/19.
//

#include "CoordinateNumPerFrame.hpp"
#include "frame.hpp"

void CoordinateNumPerFrame::process(std::shared_ptr<Frame> &frame) {

    int cn_sum = 0;
    for (auto &atom1 : group1) {
        for (auto &atom2 : group2) {
            if (atom_distance(atom1, atom2, frame) <= this->dist_cutoff) {
                cn_sum++;
            }
        }
    }
    cn_list.push_back(cn_sum);

}

void CoordinateNumPerFrame::print() {
    outfile << "***************************" << '\n';
    outfile << "******* CN per Frame ******" << '\n';
    outfile << "type 1 :" << ids1 << '\n';
    outfile << "type 2 :" << ids2 << '\n';
    outfile << "cutoff :" << dist_cutoff << '\n';
    outfile << "***************************" << '\n';
    for (auto[frame, coordination_number] : enumerate(cn_list).start(1)) {
        outfile << frame << "   " << coordination_number << '\n';
    }
    outfile << "***************************" << '\n';
}

void CoordinateNumPerFrame::readInfo() {
    Atom::select2group(ids1, ids2);
    dist_cutoff = choose(0.0, std::numeric_limits<double>::max(), "Please enter distance cutoff:");
}

void CoordinateNumPerFrame::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](std::shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->ids1)) this->group1.insert(atom);
                      if (Atom::is_match(atom, this->ids2)) this->group2.insert(atom);
                  });
}

