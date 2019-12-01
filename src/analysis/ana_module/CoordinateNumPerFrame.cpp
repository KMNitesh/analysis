//
// Created by xiamr on 6/14/19.
//

#include <boost/range/adaptors.hpp>
#include "CoordinateNumPerFrame.hpp"
#include "data_structure/frame.hpp"

CoordinateNumPerFrame::CoordinateNumPerFrame() {
    enable_outfile = true;
}

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

void CoordinateNumPerFrame::print(std::ostream &os) {
    os << "***************************" << '\n';
    os << "******* CN per Frame ******" << '\n';
    os << "type 1 :" << ids1 << '\n';
    os << "type 2 :" << ids2 << '\n';
    os << "cutoff :" << dist_cutoff << '\n';
    os << "***************************" << '\n';
    for (const auto &element: cn_list | boost::adaptors::indexed(1)) {
        os << element.index() << "   " << element.value() << '\n';
    }
    os << "***************************" << '\n';
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

