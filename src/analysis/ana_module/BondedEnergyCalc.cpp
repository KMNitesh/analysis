
#include "BondedEnergyCalc.hpp"
#include "data_structure/frame.hpp"
#include "utils/common.hpp"

BondedEnergyCalc::BondedEnergyCalc() {
    enable_forcefield = true;
    enable_outfile = true;
}

void BondedEnergyCalc::processFirstFrame([[maybe_unused]] std::shared_ptr<Frame> &frame) {
    calculator = std::make_unique<BondEnergyCalculator>(mask, frame);
}

void BondedEnergyCalc::process(std::shared_ptr<Frame> &frame) { energy_terms.push_back(calculator->energy(frame)); }

void BondedEnergyCalc::print(std::ostream &os) {

    os << '#' << title() << '\n';
    os << "# mask > " << mask << '\n';
    os << std::string(50, '#') << '\n';

    os << boost::format("#%8s %15s %15s %15s %15s %15s\n") % "Frame" % "Bond" % "Angle" % "Proper Dih." %
              "Improper Dih." % "Total";

    const boost::format fmt("%9s %15.6f %15.6f %15.6f %15.6f %15.6f\n");

    std::size_t index = 0;
    for (const auto &[bond, angle, dihedral, improper] : energy_terms) {
        os << boost::format(fmt) % ++index % bond % angle % dihedral % improper % (bond + angle + dihedral + improper);
    }
}

void BondedEnergyCalc::readInfo() { Atom::select1group(mask); }