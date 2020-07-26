

#include <cmath>

#include "HBondLifeTime.hpp"

#include "HBond.hpp"
#include "data_structure/atom.hpp"
#include "data_structure/frame.hpp"
#include "data_structure/molecule.hpp"
#include "utils/ThrowAssert.hpp"
#include "utils/common.hpp"

HBondLifeTime::HBondLifeTime() { enable_outfile = true; }

void HBondLifeTime::processFirstFrame(std::shared_ptr<Frame> &frame) {
    for (auto &atom : frame->atom_list) {
        if (is_match(atom, water_mask)) {
            switch (which(atom)) {
                case Symbol::Hydrogen:
                    hydrogens.push_back(atom);
                    break;
                case Symbol::Oxygen:
                    oxygens.push_back(atom);
                    break;
                default:
                    throw std::runtime_error("water_mask error!");
            }
        }
    }
    hb_histroy.resize(hydrogens.size());
}

void HBondLifeTime::process(std::shared_ptr<Frame> &frame) {
    int vector_shift = 0;

    for (auto &hydrogen : hydrogens) {
        auto o1 = frame->atom_map[hydrogen->con_list.front()];
        int hbond_dest_oxygen_num = 0;
        for (auto &o2 : oxygens) {
            if (o2 != o1 and atom_distance(o2, o1, frame) <= dist_R_cutoff) {
                auto o1_h_vector = hydrogen->getCoordinate() - o1->getCoordinate();
                auto o1_o2_vector = o2->getCoordinate() - o1->getCoordinate();

                frame->image(o1_h_vector);
                frame->image(o1_o2_vector);

                o1_h_vector /= vector_norm(o1_h_vector);
                o1_o2_vector /= vector_norm(o1_o2_vector);

                auto angle = radian * acos(dot_multiplication(o1_h_vector, o1_o2_vector));

                if (angle <= angle_HOO_cutoff) {
                    hbond_dest_oxygen_num = o2->seq;
                    break;
                }
            }
        }
        hb_histroy[vector_shift].push_back(hbond_dest_oxygen_num);
        ++vector_shift;
    }
}

void HBondLifeTime::print(std::ostream &os) {
    auto acf = calculateAcf();
    printData(os, acf, title());
}

void HBondLifeTime::printData(std::ostream &os, const std::vector<double> &acf, std::string_view title) const {
    os << std::string(50, '#') << '\n';
    os << "# " << title << '\n';
    os << "# dist_R_cutoff     > " << dist_R_cutoff << '\n';
    os << "# angle_HOO_cutoff  > " << angle_HOO_cutoff << '\n';
    os << "# time_increment_ps > " << time_increment_ps << '\n';
    os << "# max_time_grap_ps  > " << max_time_grap_ps << '\n';
    os << "# water_mask > " << water_mask << '\n';
    os << std::string(50, '#') << '\n';

    os << boost::format("%15s %15s\n") % "Time(ps)" % "ACF";

    for (std::size_t i = 0; i < acf.size(); ++i) {
        os << boost::format("%15.3f %15.3f\n") % (time_increment_ps * i) % acf[i];
    }

    os << std::string(50, '#') << '\n';
}

std::vector<double> HBondLifeTime::calculateAcf() const {
    auto max_time_grap_frame = std::ceil(max_time_grap_ps / time_increment_ps);
    std::vector<long> acf(std::min<int>(hb_histroy.at(0).size(), max_time_grap_frame + 1), 0);
    std::vector<long> ntime(std::min<int>(hb_histroy.at(0).size(), max_time_grap_frame + 1), 0);

    for (auto &sample : hb_histroy) {
        for (std::size_t i = 0; i < sample.size() - 1; ++i) {
            for (std::size_t j = i; j < std::min<std::size_t>(sample.size(), max_time_grap_frame + 1 + i); ++j) {
                auto n = j - i;
                assert(n < ntime.size());
                ++ntime[n];
                if (sample[i] != 0 and sample[i] == sample[j]) {
                    ++acf[n];
                }
            }
        }
    }

    std::vector<double> acff(acf.size(), 0);
    acff[0] = double(acf[0]) / ntime[0];
    for (std::size_t i = 1; i < acf.size(); ++i) {
        acff[i] = acf[i] / (ntime[i] * acff[0]);
    }
    acff[0] = 1.0;
    return acff;
}

void HBondLifeTime::readInfo() {
    dist_R_cutoff = choose(0.0, 100.0, "Distance Cutoff(O-O) for Hydogen Bond [3.5 Ang] :", Default(3.5));
    angle_HOO_cutoff = choose(0.0, 100.0, "Angle Cutoff(H-O-O) for Hydogen Bond [30 degree] :", Default(30.0));
    time_increment_ps = choose(0.0, 100.0, "time_increment_ps [0.1 ps] :", Default(0.1));
    max_time_grap_ps = choose(0.0, 100.0, "max_time_grap_ps [100 ps] :", Default(100.0));
    select1group(water_mask, "Enter mask for water atoms in system > ");
}
