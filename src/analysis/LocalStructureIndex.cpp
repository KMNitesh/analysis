//
// Created by xiamr on 9/22/19.
//

#include <boost/range/algorithm.hpp>
#include <boost/range/numeric.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/irange.hpp>
#include "LocalStructureIndex.hpp"
#include "frame.hpp"
#include "common.hpp"
#include "Histrogram2D.hpp"

LocalStructureIndex::LocalStructureIndex() {
    enable_outfile = true;
}

void LocalStructureIndex::processFirstFrame(std::shared_ptr<Frame> &frame) {
    boost::for_each(frame->atom_list,
                    [this](std::shared_ptr<Atom> &atom) {
                        if (Atom::is_match(atom, metal_mask)) metal = atom;
                        else if (Atom::is_match(atom, Ow_atom_mask)) Ow_atoms.push_back(atom);
                    });
    assert(metal);
    assert(Ow_atoms.size() > 1);
}

void LocalStructureIndex::process(std::shared_ptr<Frame> &frame) {
    std::deque<double> rs;
    for (auto &atom : Ow_atoms) {
        auto r2 = atom_distance2(metal, atom, frame);
        if (r2 < cutoff2) {
            rs.push_back(std::sqrt(r2));
        }
    }
    std::vector<double> delta_r;
    delta_r.reserve(rs.size() - 1);
    for (std::size_t i = 0; i < rs.size() - 1; ++i) {
        delta_r.push_back(rs[i + 1] - rs[i]);
    }
    auto delta_r_ = boost::accumulate(delta_r, 0.0) / delta_r.size();
    double sum = 0;
    for (auto val : delta_r) {
        sum += std::pow(val - delta_r_, 2);
    }
    localStructureIndices.emplace_back(rs.back(), sum / delta_r.size());
}

void LocalStructureIndex::print(std::ostream &os) {
    os << std::string(50, '#') << '\n';
    os << "# " << title() << " # \n";
    os << "# metal atom mask > " << metal_mask << '\n';
    os << "# coodination atom mask > " << Ow_atom_mask << '\n';
    os << "# first hydration shell cutoff(Ang) = " << std::sqrt(cutoff2) << '\n';
    os << std::string(50, '#') << '\n';
    os << boost::format("%15s %15s %15s\n") % "Frame" % "RMax(Ang)" % "LSI(Ang^2)";
    for (const auto &element: localStructureIndices | boost::adaptors::indexed(1)) {
        os << boost::format("%15d %15.6f %15.6f\n") % element.index() % element.value().first % element.value().second;
    }
    os << std::string(50, '#') << '\n';
    auto[RMax_min_iter, RMax_max_iter] = std::minmax_element(
            std::begin(localStructureIndices), std::end(localStructureIndices),
            [](auto &lhs, auto &rhs) { return std::get<0>(lhs) < std::get<0>(rhs); });
    auto[LSI_min_iter, LSI_max_iter] = std::minmax_element(
            std::begin(localStructureIndices), std::end(localStructureIndices),
            [](auto &lhs, auto &rhs) { return std::get<1>(lhs) < std::get<1>(rhs); });

    Histrogram2D histrogram2D(
            {LSI_min_iter->second, LSI_max_iter->second}, (LSI_max_iter->second - LSI_min_iter->second) / 100,
            {RMax_min_iter->first, RMax_max_iter->first}, (RMax_max_iter->first - RMax_min_iter->first) / 100);

    for (auto[rmax, lsi] : localStructureIndices) {
        histrogram2D.update(lsi, rmax);
    }

    os << boost::format("%15s %15s %15s\n") % "RMax(Ang)" % "LSI(Ang^2)" % "Normalized Probability Density";
    auto distribution = histrogram2D.getDistribution();
    auto max_population = std::get<2>(*boost::max_element(
            distribution, [](auto &lhs, auto &rhs) { return std::get<2>(lhs) < std::get<2>(rhs); }));

    for (auto[rmax, lsi, population] : distribution) {
        os << boost::format("%15.6d %15.6f %15.6f\n") % rmax % lsi % (population / max_population);
    }
}

void LocalStructureIndex::readInfo() {
    Atom::select2group(metal_mask, Ow_atom_mask, "Enter mask for center metal > ",
                       "Enter mask for coordination atom > ");
    auto cutoff = choose(0.0, "Enter cutoff for first hydration shell (Ang) [ 3.0 ] > ", Default(3.0));
    cutoff2 = cutoff * cutoff;
}
