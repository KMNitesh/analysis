//
// Created by xiamr on 9/9/19.
//

#include <boost/range/algorithm.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include "RadiusOfGyration.hpp"
#include "utils/common.hpp"
#include "data_structure/frame.hpp"
#include "data_structure/molecule.hpp"
#include "HBond.hpp"

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

    using namespace boost::accumulators;
    accumulator_set<double, features<tag::mean, tag::variance>> acc;

    for (auto &mol : moles) {
        double radius2{};
        double mol_mass{};
        auto mass_center = mol->calc_weigh_center(frame, bIncludeHydrogen);
        for (auto &atom : mol->atom_list) {
            if (!bIncludeHydrogen and which(atom) == Symbol::Hydrogen) continue;
            auto v = atom->getCoordinate() - mass_center;
            frame->image(v);
            radius2 += vector_norm2(v) * atom->mass.value();
            mol_mass += atom->mass.value();
        }

        acc(std::sqrt(radius2 / mol_mass));
    }

    series.emplace_back(mean(acc), std::sqrt(variance(acc)));
}

void RadiusOfGyration::print(std::ostream &os) {
    os << std::string(50, '#') << '\n';
    os << "# " << title() << '\n';
    os << "# AtomMask for selected molecules > " << atomMask << '\n';
    os << "# Include Hydrogen Atom > " << (bIncludeHydrogen ? 'Y' : 'N') << '\n';
    os << std::string(50, '#') << '\n';

    os << boost::format("#%15s %15s %15s\n") % "Frame" % "Rg(Ang)" % "STD(Ang)";

    const boost::format fmt{"%15.5f %15.5f %15.5f\n"};
    for (const auto rg : series | boost::adaptors::indexed(1)) {
        os << boost::format(fmt) % rg.index() % std::get<0>(rg.value()) % std::get<1>(rg.value());
    }
}

void RadiusOfGyration::readInfo() {
    Atom::select1group(atomMask, "Enter mask for selected molecules : ");
    bIncludeHydrogen = choose_bool("Include Hydrogen Atom [N] : ", Default(false));
}
