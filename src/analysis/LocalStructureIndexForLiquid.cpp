//
// Created by xiamr on 9/24/19.
//

#include <boost/range/algorithm.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/numeric.hpp>
#include <boost/range/irange.hpp>
#include "LocalStructureIndexForLiquid.hpp"
#include "frame.hpp"
#include "common.hpp"
#include "LocalStructureIndex.hpp"
#include "Histrogram2D.hpp"

LocalStructureIndexForLiquid::LocalStructureIndexForLiquid() {
    enable_outfile = true;
}

void LocalStructureIndexForLiquid::processFirstFrame(std::shared_ptr<Frame> &frame) {
    boost::for_each(frame->atom_list, [this](std::shared_ptr<Atom> &atom) {
        if (Atom::is_match(atom, center_atom_mask)) center_atoms.push_back(atom);
    });
}

void LocalStructureIndexForLiquid::process(std::shared_ptr<Frame> &frame) {
    std::deque<double> rs;
    std::vector<double> r2_total;
    for (auto &atom1 : center_atoms) {
        rs.clear();
        r2_total.clear();
        for (auto &atom2 : center_atoms) {
            if (atom1 != atom2) {
                auto r2 = atom_distance(atom1, atom2, frame);
                r2_total.push_back(r2);
                if (r2 < cutoff2) {
                    rs.push_back(std::sqrt(r2));
                }
            }
        }
        auto lsi = LocalStructureIndex::calculateLSI(rs);
        if (r2_total.size() >= r_index) {
            std::nth_element(std::begin(r2_total), std::begin(r2_total) + r_index - 1, std::end(r2_total));
            localStructureIndices.emplace_back(std::sqrt(r2_total[r_index - 1]), lsi);
        } else {
            std::cerr << "WARNING !! no r(" << r_index << ") for r distance vector \n";
        }
    }
}

void LocalStructureIndexForLiquid::print(std::ostream &os) {
    os << std::string(50, '#') << '\n';
    os << "# " << title() << " # \n";
    os << "# center atom mask > " << center_atom_mask << '\n';
    os << "# distance cutoff(Ang) = " << std::sqrt(cutoff2) << '\n';
    os << "# index for ri = " << r_index << '\n';
    os << std::string(50, '#') << '\n';

    auto[Ri_min_iter, Ri_max_iter] = std::minmax_element(
            std::begin(localStructureIndices), std::end(localStructureIndices),
            [](auto &lhs, auto &rhs) { return std::get<0>(lhs) < std::get<0>(rhs); });
    auto[LSI_min_iter, LSI_max_iter] = std::minmax_element(
            std::begin(localStructureIndices), std::end(localStructureIndices),
            [](auto &lhs, auto &rhs) { return std::get<1>(lhs) < std::get<1>(rhs); });

    Histrogram2D histrogram2D(
            {LSI_min_iter->second, LSI_max_iter->second}, (LSI_max_iter->second - LSI_min_iter->second) / 100,
            {Ri_min_iter->first, Ri_max_iter->first}, (Ri_max_iter->first - Ri_min_iter->first) / 100);

    for (auto[ri, lsi] : localStructureIndices) {
        histrogram2D.update(lsi, ri);
    }

    os << boost::format("%15s %15s %15s\n")
          % "LSI(Ang^2)"
          % ("R" + std::to_string(r_index) + "(Ang)") %
          "Normalized Probability Density";

    auto distribution = histrogram2D.getDistribution();
    auto max_population = std::get<2>(*boost::max_element(
            distribution, [](auto &lhs, auto &rhs) { return std::get<2>(lhs) < std::get<2>(rhs); }));

    for (auto[ri, lsi, population] : distribution) {
        os << boost::format("%15.6d %15.6f %15.6f\n") % ri % lsi % (population / max_population);
    }
}

void LocalStructureIndexForLiquid::readInfo() {
    Atom::select1group(center_atom_mask, "Enter mask for center atom > ");
    auto cutoff = choose(0.0, "Enter cutoff for local distance(Ang) [3.7] > ", Default(3.7));
    cutoff2 = cutoff * cutoff;
    r_index = choose(1, "Enter i for ri [5] > ", Default(5));
}



