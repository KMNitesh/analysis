
#include "KeyInteractionResidueFinder.hpp"
#include "data_structure/frame.hpp"
#include "utils/PBCUtils.hpp"
#include "utils/common.hpp"

#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm.hpp>

KeyInteractionResidueFinder::KeyInteractionResidueFinder() { enable_outfile = true; }

void KeyInteractionResidueFinder::processFirstFrame(std::shared_ptr<Frame> &frame) {
    sidechain_N_O_atoms = PBCUtils::find_atoms(sidechain_N_O_mask, frame);
    dna_P_atoms = PBCUtils::find_atoms(dna_P_mask, frame);
}

void KeyInteractionResidueFinder::process(std::shared_ptr<Frame> &frame) {
    std::map<std::string, double> residue_map;
    for (auto &atom1 : sidechain_N_O_atoms) {
        for (auto &atom2 : dna_P_atoms) {
            auto dist = atom_distance(atom1, atom2, frame);
            if (dist < cutoff) {
                auto key = atom1->residue_name.value() + std::to_string(atom1->residue_num.value());
                if (auto it = residue_map.find(key); it != std::end(residue_map)) {
                    it->second = std::min(it->second, dist);
                } else {
                    residue_map.emplace(key, dist);
                }
            }
        }
    }
    residue_distance_evolution.push_back(std::move(residue_map));
}

void KeyInteractionResidueFinder::print(std::ostream &os) {

    os << std::string(50, '#') << '\n';
    os << "# " << title() << " # \n";
    os << "# sidechain(N or O) atom mask > " << sidechain_N_O_mask << '\n';
    os << "# DNA (P) atom mask > " << dna_P_mask << '\n';
    os << "# cutoff(Ang) = " << cutoff << '\n';
    os << std::string(50, '#') << '\n';

    std::map<std::string, std::size_t> counter_map;

    for (const auto &m : residue_distance_evolution) {
        for (const auto &[k, v] : m) {
            ++counter_map[k];
        }
    }
    std::vector<std::pair<std::string, std::size_t>> counter(std::begin(counter_map), std::end(counter_map));

    boost::sort(counter, [](const auto &lhs, const auto &rhs) { return lhs.second > rhs.second; });

    os << "#Frame ";
    for (const auto &p : counter) {
        os << std::setw(15) << p.first;
    }
    os << '\n';
    os << std::string(7 + 15 * counter.size(), '#') << '\n';
    os << "# Count";
    for (const auto &p : counter) {
        os << std::setw(15) << p.second;
    }
    os << '\n';
    os << std::string(7 + 15 * counter.size(), '#') << '\n';

    for (const auto &element : residue_distance_evolution | boost::adaptors::indexed(1)) {
        os << std::setw(7) << element.index();
        const auto &dists = element.value();
        for (const auto &[key, value] : counter) {
            if (auto it = dists.find(key); it != dists.end()) {
                os << std::setw(15) << it->second;
            } else {
                os << std::setw(15) << "";
            }
        }
        os << '\n';
    }
}

void KeyInteractionResidueFinder::readInfo() {

    select1group(sidechain_N_O_mask, "Sidechain(N or O) > ");
    select1group(dna_P_mask, "DNA (P) > ");

    cutoff = choose(0.0, "Cutoff between sidechain(N or O) and DNA(P) [4.7Ang] >", Default(4.7));
}