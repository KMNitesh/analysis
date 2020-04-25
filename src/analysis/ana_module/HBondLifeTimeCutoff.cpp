//
// Created by xiamr on 8/9/19.
//

#include "HBondLifeTimeCutoff.hpp"

#include "HBond.hpp"
#include "data_structure/frame.hpp"
#include "data_structure/molecule.hpp"
#include "utils/ThrowAssert.hpp"
#include "utils/common.hpp"

bool operator==(const HBondLifeTimeCutoff::InnerAtom &i1, const HBondLifeTimeCutoff::InnerAtom &i2) {
    return i1.index == i2.index and i1.list_ptr1 == i2.list_ptr1 and i1.list_ptr2 == i2.list_ptr2;
}

auto HBondLifeTimeCutoff::find_in(int seq) {
    return std::find_if(inner_atoms.begin(), inner_atoms.end(), [seq](auto &i) { return i.index == seq; });
}

HBondLifeTimeCutoff::HBondLifeTimeCutoff() { enable_outfile = true; }

void HBondLifeTimeCutoff::processFirstFrame(std::shared_ptr<Frame> &frame) {
    for (auto &mol : frame->molecule_list) {
        for (auto &atom : mol->atom_list) {
            if (is_match(atom, Ow_atom_mask)) {
                std::deque<std::shared_ptr<Atom>> hydrogens;
                for (auto &hydro : mol->atom_list) {
                    if (which(hydro) == Symbol::Hydrogen) {
                        hydrogens.push_back(hydro);
                    }
                }
                water_struct.emplace_back(atom, hydrogens);
            } else if (is_match(atom, center_Metal_atom_mask)) {
                metal.emplace(atom);
            }
        }
    }
    throw_assert(metal.size() == 1, "Only one metal atom is allowed!");
}

void HBondLifeTimeCutoff::process(std::shared_ptr<Frame> &frame) {
    for (auto &ref : metal) {
        for (auto &[Ow, hydrogens] : water_struct) {
            auto it = find_in(Ow->molecule.lock()->seq());
            if (atom_distance2(ref, Ow, frame) < cutoff2) {
                // in the shell
                auto &o1 = Ow;
                if (it != inner_atoms.end()) {
                    it->list_ptr1->push_back(find_hbond(o1, hydrogens[0], frame));
                    it->list_ptr2->push_back(find_hbond(o1, hydrogens[1], frame));
                } else {
                    auto list_ptr1 = std::make_shared< std::deque<int>>();
                    auto list_ptr2 = std::make_shared< std::deque<int>>();
                    list_ptr1->push_back(find_hbond(o1, hydrogens[0], frame));
                    list_ptr2->push_back(find_hbond(o1, hydrogens[1], frame));

                    inner_atoms.insert(InnerAtom(Ow->molecule.lock()->seq(), list_ptr1, list_ptr2));
                    hb_histroy.emplace_back(list_ptr1);
                    hb_histroy.emplace_back(list_ptr2);
                }
            } else {
                if (it != inner_atoms.end()) {
                    inner_atoms.erase(it);
                }
            }
        }
    }
}

int HBondLifeTimeCutoff::find_hbond(const std::shared_ptr<Atom> &o1, const std::shared_ptr<Atom> &hydrogen,
                                    std::shared_ptr<Frame> &frame) const {
    int hbond_dest_oxygen_num = 0;
    for (auto &ele : water_struct) {
        auto &o2 = ele.first;
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
    return hbond_dest_oxygen_num;
}

void HBondLifeTimeCutoff::printData(std::ostream &os, const std::vector<double> &acf, std::string_view title) const {
    os << std::string(50, '#') << '\n';
    os << "# " << title << '\n';
    os << "# dist_R_cutoff    (Ang) > " << dist_R_cutoff << '\n';
    os << "# angle_HOO_cutoff (Ang) > " << angle_HOO_cutoff << '\n';
    os << "# time_increment_ps (ps) > " << time_increment_ps << '\n';
    os << "# max_time_grap_ps  (ps) > " << max_time_grap_ps << '\n';
    os << "# Ow_mask                > " << Ow_atom_mask << '\n';
    os << "# metal_mask             > " << center_Metal_atom_mask << '\n';
    os << "# cutoff          (Ang)  > " << std::sqrt(cutoff2) << '\n';
    os << std::string(50, '#') << '\n';

    os << boost::format("%15s %15s\n") % "Time(ps)" % "ACF";

    for (std::size_t i = 0; i < acf.size(); ++i) {
        os << boost::format("%15.3f %15.3f\n") % (time_increment_ps * i) % acf[i];
    }

    os << std::string(50, '#') << '\n';
}

void HBondLifeTimeCutoff::print(std::ostream &os) {
    auto acf = calculateAcf();
    printData(os, acf, title());
}

std::vector<double> HBondLifeTimeCutoff::calculateAcf() const {
    auto max_time_grap_frame = ceil(max_time_grap_ps / time_increment_ps);
    std::vector<long> acf(1, 0);
    std::vector<long> ntime(1, 0);

    for (auto &sample : hb_histroy) {
        for (std::size_t i = 0; i < sample->size() - 1; ++i) {
            for (std::size_t j = i; j < std::min<std::size_t>(sample->size(), max_time_grap_frame + 1 + i); ++j) {
                auto n = j - i;
                if (ntime.size() <= n) {
                    ntime.push_back(0);
                    acf.push_back(0);
                }
                ++ntime[n];
                if ((*sample)[i] != 0 and (*sample)[i] == (*sample)[j]) {
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

void HBondLifeTimeCutoff::readInfo() {
    dist_R_cutoff = choose(0.0, 100.0, "Distance Cutoff(O-O) for Hydogen Bond [3.5 Ang] :", Default(3.5));
    angle_HOO_cutoff = choose(0.0, 100.0, "Angle Cutoff(H-O-O) for Hydogen Bond [30 degree] :", Default(30.0));
    time_increment_ps = choose(0.0, 100.0, "time_increment_ps [0.1 ps] :", Default(0.1));
    max_time_grap_ps = choose(0.0, 100.0, "max_time_grap_ps [100 ps] :", Default(100.0));
    auto cutoff = choose(0.0, 100.0, "cutoff :");
    cutoff2 = cutoff * cutoff;
    select1group(center_Metal_atom_mask, "Enter mask for Metal > ");
    select1group(Ow_atom_mask, "Enter mask for Ow > ");
}
