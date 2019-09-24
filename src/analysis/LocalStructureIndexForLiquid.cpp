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
#include "Histrogram2D.hpp"
#include "molecule.hpp"

LocalStructureIndexForLiquid::LocalStructureIndexForLiquid() {
    enable_outfile = true;
    enable_forcefield = true;
}


void LocalStructureIndexForLiquid::process(std::shared_ptr<Frame> &frame) {
    std::vector<std::pair<double, std::vector<double>>> distance_within_cutoff_vector_pair;
    std::vector<double> distance_square_vector;
    std::vector<double> distance_within_cutoff_range;
    for (auto &mol1 : frame->molecule_list) {
        distance_square_vector.clear();
        distance_within_cutoff_range.clear();
        auto mol1_coord = mol1->calc_weigh_center(frame);
        for (auto &mol2 : frame->molecule_list) {
            if (mol1 != mol2) {
                auto mol2_coord = mol2->calc_weigh_center(frame);
                auto v = mol2_coord - mol1_coord;
                frame->image(v);
                auto r2 = vector_norm2(v);
                distance_square_vector.push_back(r2);
                if (r2 < cutoff2) {
                    distance_within_cutoff_range.push_back(std::sqrt(r2));
                }
            }
        }
        boost::stable_sort(distance_within_cutoff_range);
        std::vector<double> distance_differences;
        distance_differences.reserve(distance_within_cutoff_range.size() - 1);

        for (std::size_t i = 0; i < distance_within_cutoff_range.size() - 1; ++i) {
            distance_differences.push_back(distance_within_cutoff_range[i + 1] - distance_within_cutoff_range[i]);
        }
        if (distance_square_vector.size() >= r_index) {
            std::nth_element(std::begin(distance_square_vector), std::begin(distance_square_vector) + r_index - 1,
                             std::end(distance_square_vector));
            distance_within_cutoff_vector_pair.emplace_back(std::sqrt(distance_square_vector[r_index - 1]),
                                                            std::move(distance_differences));
        } else {
            std::cerr << "WARNING !! no r(" << r_index << ") for r distance vector \n";
        }
    }
    double distance_difference_average = 0.0;
    std::size_t total_count = 0;
    for (auto &elemenet : distance_within_cutoff_vector_pair) {
        distance_difference_average = boost::accumulate(elemenet.second, distance_difference_average);
        total_count += elemenet.second.size();
    }

    distance_difference_average /= total_count;

    for (auto &[ri, distance_differences] : distance_within_cutoff_vector_pair) {
        double lsi_sum = 0;
        for (auto val : distance_differences) {
            lsi_sum += std::pow(val - distance_difference_average, 2);
        }
        localStructureIndices.emplace_back(ri, lsi_sum / distance_differences.size());
    }
}

void LocalStructureIndexForLiquid::print(std::ostream &os) {
    os << std::string(50, '#') << '\n';
    os << "# " << title() << " # \n";
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
          % ("R" + std::to_string(r_index) + "(Ang)")
          % "Normalized Probability Density";

    auto distribution = histrogram2D.getDistribution();
    auto max_population = std::get<2>(*boost::max_element(
            distribution, [](auto &lhs, auto &rhs) { return std::get<2>(lhs) < std::get<2>(rhs); }));

    for (auto[ri, lsi, population] : distribution) {
        os << boost::format("%15.6d %15.6f %15.6f\n") % ri % lsi % (population / max_population);
    }
}

void LocalStructureIndexForLiquid::readInfo() {
    auto cutoff = choose(0.0, "Enter cutoff for local distance(Ang) [3.7] > ", Default(3.7));
    cutoff2 = cutoff * cutoff;
    r_index = choose(1, "Enter i for ri [5] > ", Default(5));
}



