//
// Created by xiamr on 6/14/19.
//
#include <boost/range/adaptors.hpp>
#include "Distance.hpp"
#include "frame.hpp"


void Distance::process(std::shared_ptr<Frame> &frame) {

    double x1, y1, z1, x2, y2, z2;
    x1 = y1 = z1 = x2 = y2 = z2 = 0.0;
    double weigh1, weigh2;
    weigh1 = weigh2 = 0.0;
    for (auto &atom1 : group1) {
        assert(atom1->mass);
        double mass = atom1->mass.get();
        x1 += atom1->x * mass;
        y1 += atom1->y * mass;
        z1 += atom1->z * mass;
        weigh1 += mass;
    }
    x1 /= weigh1;
    y1 /= weigh1;
    z1 /= weigh1;

    for (auto &atom2 : group2) {
        assert(atom2->mass);
        double mass = atom2->mass.get();
        x2 += atom2->x * mass;
        y2 += atom2->y * mass;
        z2 += atom2->z * mass;
        weigh2 += mass;
    }
    x2 /= weigh2;
    y2 /= weigh2;
    z2 /= weigh2;

    double xr = x2 - x1;
    double yr = y2 - y1;
    double zr = z2 - z1;

    frame->image(xr, yr, zr);

    double dist = sqrt(xr * xr + yr * yr + zr * zr);
    group_dist_list.push_back(dist);
}


void Distance::print(std::ostream &os) {
    os << "****************************************\n";
    os << "Distance between \n";
    os << "group 1 :" << ids1 << std::endl;
    os << "group 2 :" << ids2 << std::endl;
    os << "****************************************\n";

    for (const auto &element : group_dist_list | boost::adaptors::indexed(1)) {
        os << element.index() << "    " << element.value() << '\n';
    }
}

void Distance::readInfo() {
    enable_forcefield = true;
    Atom::select2group(ids1, ids2);
}

void Distance::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](std::shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->ids1)) this->group1.insert(atom);
                      if (Atom::is_match(atom, this->ids2)) this->group2.insert(atom);
                  });
}
// End of Distance
