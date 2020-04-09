
#include "others/MMPBSA.hpp"
#include "trajectory_reader/trajectoryreader.hpp"
#include "utils/common.hpp"
#include <boost/algorithm/algorithm.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/tokenizer.hpp>

std::istream &operator>>(std::istream &is,
                         std::vector<std::pair<MMPBSA::Residue, MMPBSA::ResidueComponent>> &all_residues_total) {
    std::string line;
    bool current_is_total_componet = false;
    int lengths[] = {1, 4, 4};
    boost::offset_separator ofs(std::begin(lengths), std::end(lengths));
    while (!is.eof()) {
        std::getline(is, line);
        boost::trim(line);
        if (line.empty())
            break;
        if (current_is_total_componet) {
            auto fields = split(line, ",");
            boost::tokenizer tokenizer(fields[1], ofs);
            auto it = std::begin(tokenizer);
            MMPBSA::Residue r;
            r.location = *it++ == "R" ? MMPBSA::Residue::R : MMPBSA::Residue::L;
            r.name = boost::trim_copy(*it++);
            r.no = std::stoi(*it);

            all_residues_total.emplace_back(std::move(r), MMPBSA::ResidueComponent{});
            all_residues_total.back().second.internal = {fields[2], fields[3], fields[4]};
            all_residues_total.back().second.vdW = {fields[5], fields[6], fields[7]};
            all_residues_total.back().second.electrostatic = {fields[8], fields[9], fields[10]};
            all_residues_total.back().second.polar_solv = {fields[11], fields[12], fields[13]};
            all_residues_total.back().second.nonpolar_solv = {fields[14], fields[15], fields[16]};
            all_residues_total.back().second.total = {fields[17], fields[18], fields[19]};
        } else {
            if (boost::starts_with(line, "Total Energy Decomposition:")) {
                std::getline(is, line);
                std::getline(is, line);
                current_is_total_componet = true;
            }
        }
    }
    return is;
}

std::ostream &operator<<(std::ostream &os,
                         const std::vector<std::pair<MMPBSA::Residue, MMPBSA::ResidueComponent>> all_residues_total) {
    os << boost::format("%-10s  %10s %10s %10s %10s %10s\n") % "Residue" % "vdW" % "Ele" % "PolarSolv" % "NonPol_Sov" %
              "Total";
    const boost::format fmt{"%-10s %10.3f %10.3f %10.3f %10.3f %10.3f\n"};
    for (const auto &[res, e] : all_residues_total) {
        os << boost::format(fmt) % (res.name + std::to_string(res.no)) % e.vdW.avg % e.electrostatic.avg %
                  e.polar_solv.avg % e.nonpolar_solv.avg % e.total.avg;
    }
    return os;
}

std::ostream &operator<<(std::ostream &os,
                         const std::vector<std::pair<std::string, MMPBSA::ResidueComponent>> all_residues_total) {
    os << boost::format("%-10s  %10s %10s %10s %10s %10s\n") % "Residue" % "vdW" % "Ele" % "PolarSolv" % "NonPol_Sov" %
              "Total";
    const boost::format fmt{"%-10s %10.3f %10.3f %10.3f %10.3f %10.3f\n"};
    for (const auto &[res, e] : all_residues_total) {
        os << boost::format(fmt) % res % e.vdW.avg % e.electrostatic.avg % e.polar_solv.avg % e.nonpolar_solv.avg %
                  e.total.avg;
    }
    return os;
}

void MMPBSA::process(const std::string &topology_filename) {

    std::vector<std::pair<Residue, ResidueComponent>> all_residues_total;

    std::string filename = choose_file("MMPBSA Decomposition data > ").isExist(true).extension("dat");

    std::ifstream ifs(filename);

    ifs >> all_residues_total;

    boost::sort(all_residues_total,
                [](const auto &lhs, const auto &rhs) { return lhs.second.total.avg < rhs.second.total.avg; });

    std::cout << all_residues_total;

    std::cout << std::string(50, '#') << '\n';
    std::cout << "Result Combines\n";
    std::cout << std::string(50, '#') << '\n';

    std::map<std::string, ResidueComponent> combine_map;

    TrajectoryReader reader;
    reader.set_topology(topology_filename);
    auto frame = reader.readTopology();

    for (auto &[residue, e] : all_residues_total) {
        if (residue.location == Residue::R) {
            if (auto it = frame->real_residue_num_map.find(residue.no); it != end(frame->real_residue_num_map)) {
                auto key = it->second->residue_name.get() + std::to_string(it->second->residue_num.get());
                if (auto it2 = combine_map.find(key); it2 != end(combine_map)) {
                    it2->second += e;
                } else {
                    combine_map.emplace(key, e);
                }
            }
        }
    }
    std::vector<std::pair<std::string, ResidueComponent>> combine_vector(begin(combine_map), end(combine_map));
    boost::sort(combine_vector,
                [](const auto &lhs, const auto &rhs) { return lhs.second.total.avg < rhs.second.total.avg; });

    std::cout << combine_vector;
}