//
// Created by xiamr on 9/22/19.
//

#include "LocalStructureIndex.hpp"

#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/numeric.hpp>

#include "data_structure/frame.hpp"
#include "utils/Histogram2D.hpp"
#include "utils/common.hpp"

LocalStructureIndex::LocalStructureIndex() { enable_outfile = true; }

void LocalStructureIndex::processFirstFrame(std::shared_ptr<Frame> &frame) {
    boost::for_each(frame->atom_list, [this](std::shared_ptr<Atom> &atom) {
        if (Atom::is_match(atom, metal_mask))
            metal = atom;
        else if (Atom::is_match(atom, Ow_atom_mask))
            Ow_atoms.push_back(atom);
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

    auto lsi = calculateLSI(rs);
    localStructureIndices.emplace_back(rs.back(), lsi);
}

template <typename RandomAccessRange>
double LocalStructureIndex::calculateLSI(RandomAccessRange &distance_within_cutoff_range) {
    boost::stable_sort(distance_within_cutoff_range);
    std::vector<double> distance_differences;
    distance_differences.reserve(distance_within_cutoff_range.size() - 1);

    for (std::size_t i = 0; i < distance_within_cutoff_range.size() - 1; ++i) {
        distance_differences.push_back(distance_within_cutoff_range[i + 1] - distance_within_cutoff_range[i]);
    }
    auto distance_difference_average = boost::accumulate(distance_differences, 0.0) / distance_differences.size();
    double sum = 0;
    for (auto val : distance_differences) {
        sum += std::pow(val - distance_difference_average, 2);
    }
    return sum / distance_differences.size();
}

template double LocalStructureIndex::calculateLSI(std::deque<double> &);

template double LocalStructureIndex::calculateLSI(std::vector<double> &);

void LocalStructureIndex::print(std::ostream &os) {
    os << std::string(50, '#') << '\n';
    os << "# " << title() << " # \n";
    os << "# metal atom mask > " << metal_mask << '\n';
    os << "# coodination atom mask > " << Ow_atom_mask << '\n';
    os << "# first hydration shell cutoff(Ang) = " << std::sqrt(cutoff2) << '\n';
    os << std::string(50, '#') << '\n';
    os << boost::format("%15s %15s %15s\n") % "Frame" % "RMax(Ang)" % "LSI(Ang^2)";
    for (const auto &element : localStructureIndices | boost::adaptors::indexed(1)) {
        os << boost::format("%15d %15.6f %15.6f\n") % element.index() % element.value().first % element.value().second;
    }
    os << std::string(50, '#') << '\n';
    auto [RMax_min_iter, RMax_max_iter] =
        std::minmax_element(std::begin(localStructureIndices), std::end(localStructureIndices),
                            [](auto &lhs, auto &rhs) { return std::get<0>(lhs) < std::get<0>(rhs); });
    auto [LSI_min_iter, LSI_max_iter] =
        std::minmax_element(std::begin(localStructureIndices), std::end(localStructureIndices),
                            [](auto &lhs, auto &rhs) { return std::get<1>(lhs) < std::get<1>(rhs); });

    Histogram2D histrogram2D(
        {LSI_min_iter->second, LSI_max_iter->second}, (LSI_max_iter->second - LSI_min_iter->second) / 100,
        {RMax_min_iter->first, RMax_max_iter->first}, (RMax_max_iter->first - RMax_min_iter->first) / 100);

    for (auto [rmax, lsi] : localStructureIndices) {
        histrogram2D.update(lsi, rmax);
    }

    os << boost::format("%15s %15s %15s\n") % "LSI(Ang^2)" % "RMax(Ang)" % "Normalized Probability Density";
    auto distribution = histrogram2D.getDistribution();
    auto max_population = std::get<2>(
        *boost::max_element(distribution, [](auto &lhs, auto &rhs) { return std::get<2>(lhs) < std::get<2>(rhs); }));

    for (auto [rmax, lsi, population] : distribution) {
        os << boost::format("%15.6d %15.6f %15.6f\n") % rmax % lsi % (population / max_population);
    }
}

void LocalStructureIndex::readInfo() {
    Atom::select2group(metal_mask, Ow_atom_mask, "Enter mask for center metal > ",
                       "Enter mask for coordination atom > ");
    auto cutoff = choose(0.0, "Enter cutoff for first hydration shell (Ang) [ 3.0 ] > ", Default(3.0));
    cutoff2 = cutoff * cutoff;
}
