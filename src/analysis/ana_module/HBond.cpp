//
// Created by xiamr on 6/14/19.
//

#include "HBond.hpp"
#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/algorithm.hpp>

#include "data_structure/atom.hpp"
#include "data_structure/forcefield.hpp"
#include "data_structure/frame.hpp"
#include "utils/common.hpp"

using namespace std;

HBond::HBond() {
    enable_outfile = true;
    enable_forcefield = true;
}

Symbol which(const std::shared_ptr<Atom> &atom) {
    double mass = forcefield.find_mass(atom);
    if (mass < 2.0)
        return Symbol::Hydrogen;
    else if (mass > 11.5 and mass < 13.0)
        return Symbol::Carbon;
    else if (mass > 13.0 and mass < 15.0)
        return Symbol::Nitrogen;
    else if (mass >= 15.0 and mass < 17.0)
        return Symbol::Oxygen;
    else if (mass >= 30.0 and mass < 31.5)
        return Symbol::Phosphorus;
    else if (mass >= 31.5 and mass < 33.0)
        return Symbol::Sulfur;
    else if (mass >= 22.5 and mass < 23.5)
        return Symbol::Sodium;
    else
        return Symbol::Unknown;
}

void HBond::print(std::ostream &os) {
    os << std::string(20, '#') << '\n';
    os << '#' << title() << '\n';
    if (mode == Selector::Both) {
        os << "#  Donor and Acceptor > " << mask1 << '\n';
    } else {
        os << "#    Donor > " << mask1 << '\n';
        os << "# Acceptor > " << mask2 << '\n';
    }
    switch (hbond_type) {
    case HBondType::VMDVerion:
        os << "#HBond criteria : VMD version\n";
        break;
    case HBondType::GMXVersion:
        os << "#HBond criteria : GMX version\n";
        break;
    }
    os << "#distance cutoff : " << donor_acceptor_dist_cutoff << '\n';
    os << "#angle cutoff : " << angle_cutoff << '\n';
    os << std::string(20, '#') << '\n';
    // os << boost::format("#%14s %14s\n") % "Frame" % "Number";
    // for (const auto &element : hbdist | boost::adaptors::indexed(1)) {
    //     os << boost::format("%15d %14d\n") % element.index() % element.value().size();
    // }
    // os << std::string(20, '#') << '\n';
    output_distance_statistics(os);
}

void HBond::process(std::shared_ptr<Frame> &frame) {
    // hbdist.emplace_back();
    residue_map.emplace_back();
    if (mode == Selector::Both)
        Selector_Both(frame);
    else
        Selector_Donor_Acceptor(frame);
}

void HBond::readInfo() {
    select1group(mask1, "Please enter group:");

    mode = static_cast<Selector>(choose(1, 3, "Which selector: [(1)Acceptor | (2)Donor | (3)Both]:") - 1);

    if (this->mode not_eq Selector::Both) {
        select1group(mask2, "Please enter group2:");
    }
    cout << "HBond criteria:\n(1)VMD version\n(2)GMX version\n";
    hbond_type = static_cast<HBondType>(choose(1, 2, "choose:") - 1);
    if (mode == Selector::Acceptor) {
        mode = Selector::Donor;
        std::swap(mask1, mask2);
    }
    donor_acceptor_dist_cutoff = choose(0.0, "Donor-Acceptor Distance:");
    angle_cutoff = choose(0.0, "Angle cutoff:");
}

void HBond::Selector_Donor_Acceptor(const std::shared_ptr<Frame> &frame) {
    for (const auto &[donor, hydrogen] : donor_hydrogens) {
        for (const auto &acceptor : acceptors) {
            check_hbond(donor, hydrogen, acceptor, frame);
        }
    }
}

void HBond::Selector_Both(const std::shared_ptr<Frame> &frame) {
    for (const auto &[donor, hydrogen] : donor_hydrogens) {
        for (const auto &acceptor : acceptors) {
            if (donor not_eq acceptor)
                check_hbond(donor, hydrogen, acceptor, frame);
        }
    }
}

void HBond::processFirstFrame(std::shared_ptr<Frame> &frame) {
    for (const auto &atom : frame->atom_list) {
        if (is_match(atom, mask1)) {
            auto symbol = which(atom);
            if (symbol == Symbol::Hydrogen) {
                const auto &donor = frame->atom_map[atom->con_list.front()];
                if (donor_acceptor_symbols.contains(which(donor))) {
                    donor_hydrogens.emplace_back(std::array{donor, atom});
                }
            } else if (mode == Selector::Both and donor_acceptor_symbols.contains(symbol))
                acceptors.push_back(atom);
        }
        if (mode not_eq Selector::Both and is_match(atom, mask2) and donor_acceptor_symbols.contains(which(atom)))
            acceptors.push_back(atom);
    };
}

bool HBond::check_hbond(const std::shared_ptr<Atom> &donor, const std::shared_ptr<Atom> &hydrogen,
                        const std::shared_ptr<Atom> &acceptor, const std::shared_ptr<Frame> &frame) {
    if (atom_distance(donor, acceptor, frame) <= donor_acceptor_dist_cutoff) {
        auto v1 = hydrogen->getCoordinate() - donor->getCoordinate();
        frame->image(v1);
        auto len1 = vector_norm2(v1);

        auto v2 = acceptor->getCoordinate() -
                  (hbond_type == HBondType::GMXVersion ? donor->getCoordinate() : hydrogen->getCoordinate());
        frame->image(v2);
        auto len2 = vector_norm2(v2);

        auto cosine = dot_multiplication(v1, v2) / std::sqrt(len1 * len2);
        if (std::abs(radian * std::acos(cosine)) <= angle_cutoff) {
            record_hbond(donor, acceptor, frame);
            return true;
        }
    }
    return false;
}

void HBond::record_hbond(const std::shared_ptr<Atom> &donor, const std::shared_ptr<Atom> &acceptor,
                         const std::shared_ptr<Frame> &frame) {

    auto distance = atom_distance(donor, acceptor, frame);
    // hbdist.back()[std::array{donor, acceptor}] = distance;

    auto key = donor->residue_name.value() + std::to_string(donor->residue_num.value());
    auto &last = residue_map.back();
    if (auto it = last.find(key); it != std::end(last)) {
        it->second = std::min(it->second, distance);
    } else {
        last.emplace(key, distance);
    }
}

void HBond::output_distance_statistics(std::ostream &os) {
    std::map<std::string, std::size_t> countor_map;
    for (auto &seq : residue_map) {
        for (auto &[key, value] : seq) {
            ++countor_map[key];
        }
    }
    std::vector<std::pair<std::string, std::size_t>> countor(begin(countor_map), end(countor_map));

    boost::sort(countor, [](const auto &lhs, const auto &rhs) { return lhs.second > rhs.second; });

    os << "#Frame ";
    for (const auto &p : countor) {
        os << std::setw(15) << p.first;
    }
    os << '\n';
    os << std::string(7 + 15 * countor.size(), '#') << '\n';
    os << "# Count";
    for (const auto &p : countor) {
        os << std::setw(15) << p.second;
    }
    os << '\n';
    os << std::string(7 + 15 * countor.size(), '#') << '\n';

    for (const auto &element : residue_map | boost::adaptors::indexed(1)) {
        os << std::setw(7) << element.index();
        const auto &hb_distances = element.value();
        for (const auto &[key, value] : countor) {
            if (auto it = hb_distances.find(key); it != hb_distances.end()) {
                os << std::setw(15) << it->second;
            } else {
                os << std::setw(15) << "";
            }
        }
        os << '\n';
    }
}
