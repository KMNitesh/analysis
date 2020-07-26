
#include <nlohmann/json.hpp>

#include <boost/range/algorithm.hpp>
#include <boost/range/numeric.hpp>

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

void BondedEnergyCalc::process(std::shared_ptr<Frame> &frame) {
    auto map = calculator->energy_with_residue_rank(frame);
    energy_terms.push_back(std::move(map));
}

void BondedEnergyCalc::print(std::ostream &os) {
    os << '#' << title() << '\n';
    os << "# mask > " << mask << '\n';
    os << std::string(50, '#') << '\n';
    os << "# Unit : kcal/mol\n";

    os << boost::format("#%8s %15s %15s %15s %15s %15s\n") % "Frame" % "Bond" % "Angle" % "Proper Dih." %
              "Improper Dih." % "Total";

    const boost::format fmt("%9d %15.6f %15.6f %15.6f %15.6f %15.6f\n");

    std::size_t index = 0;
    for (auto &frame : energy_terms) {
        const auto &[bond, angle, dihedral, improper] = boost::accumulate(
            frame, BondEnergyCalculator::Term{},
            [](const BondEnergyCalculator::Term &sum, const std::pair<int, BondEnergyCalculator::Term> &term) {
                return BondEnergyCalculator::Term{.bond = sum.bond + term.second.bond,
                                                  .angle = sum.angle + term.second.angle,
                                                  .dihedral = sum.dihedral + term.second.dihedral,
                                                  .improper = sum.improper + term.second.improper};
            });

        os << boost::format(fmt) % ++index % bond % angle % dihedral % improper % (bond + angle + dihedral + improper);
    }
    os << std::string(50, '#') << '\n';
    index = 0;
    os << boost::format("#%8s %5s %15s %15s %15s %15s %15s\n") % "Frame" % "Resid" % "Bond" % "Angle" % "Proper Dih." %
              "Improper Dih." % "Total";
    const boost::format fmt2("%9d %5d %15.6f %15.6f %15.6f %15.6f %15.6f\n");
    for (auto &frame : energy_terms) {
        auto it = boost::max_element(
            frame, [](const auto &lhs, const auto &rhs) { return lhs.second.total() < rhs.second.total(); });
        const auto [resid, term] = *it;
        const auto &[bond, angle, dihedral, improper] = term;
        os << boost::format(fmt2) % ++index % resid % bond % angle % dihedral % improper %
                  (bond + angle + dihedral + improper);
    }

    nlohmann::json j{{"title", to_string(mask)}, {"data", energy_terms}};

    std::ofstream file("bonded_energy-dump.json");
    file << std::setw(2) << j;
}

void BondedEnergyCalc::readInfo() { select1group(mask); }