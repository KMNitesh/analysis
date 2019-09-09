//
// Created by xiamr on 9/9/19.
//

#include <boost/range/algorithm.hpp>
#include <boost/range/adaptors.hpp>
#include "RadiusOfGyration.hpp"
#include "common.hpp"
#include "frame.hpp"
#include "molecule.hpp"

RadiusOfGyration::RadiusOfGyration() {
    enable_outfile = true;
    enable_forcefield = true;
}

void RadiusOfGyration::processFirstFrame(std::shared_ptr<Frame> &frame) {
    boost::for_each(frame->atom_list, [this](std::shared_ptr<Atom> &atom) {
        if (Atom::is_match(atom, atomMask)) {
            moles.insert(atom->molecule.lock());
        }
    });

    assert(!moles.empty());
}

void RadiusOfGyration::process(std::shared_ptr<Frame> &frame) {

    double total_radius{};
    int ntime = 0;
    for (auto &mol : moles) {
        auto mass_center = mol->calc_weigh_center(frame);
        double radius{};
        for (auto &atom : mol->atom_list) {
            auto v = atom->getCoordinate() - mass_center;
            frame->image(v);
            radius += vector_norm(v) * atom->mass.value();
        }
        mol->calc_mass();
        radius /= mol->mass;
        ++ntime;
        total_radius += radius;
    }

    series.push_back(total_radius / ntime);
}

void RadiusOfGyration::print(std::ostream &os) {
    os << std::string(50, '#') << '\n';
    os << "# " << title() << '\n';
    os << "# AtomMask for selected molecules > " << atomMask << '\n';
    os << std::string(50, '#') << '\n';

    os << format("#%15s %15s\n", "Frame", "Rg(Ang)");

    for (const auto rg : series | boost::adaptors::indexed(1)) {
        os << boost::format("%15.5f %15.5f\n") % rg.index() % rg.value();
    }
}

void RadiusOfGyration::readInfo() {
    Atom::select1group(atomMask, "Enter mask for selected molecules : ");
}
