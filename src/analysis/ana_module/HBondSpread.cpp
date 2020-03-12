//
// Created by xiamr on 8/11/19.
//

#include "HBondSpread.hpp"

#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm.hpp>

#include "HBond.hpp"
#include "data_structure/frame.hpp"
#include "data_structure/molecule.hpp"
#include "utils/ThrowAssert.hpp"
#include "utils/common.hpp"

HBondSpread::HBondSpread() { enable_outfile = true; }

void HBondSpread::processFirstFrame(std::shared_ptr<Frame> &frame) {
    for (auto &mol : frame->molecule_list) {
        for (auto &atom : mol->atom_list) {
            if (Atom::is_match(atom, Ow_atom_mask)) {
                Ow_mapping[atom->seq] = boost::add_vertex(atom->seq, g);
                Ow.emplace_back(atom);
                for (auto &hydro : mol->atom_list) {
                    if (which(hydro) == Symbol::Hydrogen) {
                        hydrogens.emplace_back(hydro);
                    }
                }
            } else if (Atom::is_match(atom, center_Metal_atom_mask)) {
                metal.emplace(atom);
            }
        }
    }
    throw_assert(metal.size() == 1, "Only one metal atom is allowed!");
}

void HBondSpread::process(std::shared_ptr<Frame> &frame) {
    for (auto &hydrogen : hydrogens) {
        auto o1 = frame->atom_map[hydrogen->con_list.front()];
        for (auto &o2 : Ow) {
            if (o2 != o1 and atom_distance(o2, o1, frame) <= dist_R_cutoff) {
                auto o1_h_vector = hydrogen->getCoordinate() - o1->getCoordinate();
                auto o1_o2_vector = o2->getCoordinate() - o1->getCoordinate();

                frame->image(o1_h_vector);
                frame->image(o1_o2_vector);

                o1_h_vector /= vector_norm(o1_h_vector);
                o1_o2_vector /= vector_norm(o1_o2_vector);

                auto angle = radian * acos(dot_multiplication(o1_h_vector, o1_o2_vector));

                if (angle <= angle_HOO_cutoff) {
                    boost::add_edge(Ow_mapping[o1->seq], Ow_mapping[o2->seq], g);
                    //     std::cout << o1->seq << " <-> " << o2->seq << '\n';
                    break;
                }
            }
        }
    }
    auto &m = *metal.begin();
    std::unordered_set<int> hbond_connected;

    for (auto &atom : Ow) {
        if (atom_distance2(m, atom, frame) < cutoff2) {
            MyVisitor vis(hbond_connected);
            boost::depth_first_visit(
                g, Ow_mapping[atom->seq], vis,
                boost::make_vector_property_map<boost::default_color_type>(boost::get(boost::vertex_index, g)));
            // std::cout << atom->seq << '\n';
        }
    }

    double max_distance2 = 0;

    for (auto &atom : hbond_connected) {
        max_distance2 = std::max(max_distance2, atom_distance2(m, frame->atom_map[atom], frame));
    }

    hbond_connected_num_and_max_distance.emplace_back(hbond_connected.size(), std::sqrt(max_distance2));

    boost::remove_edge_if([](auto) { return true; }, g);
}

void HBondSpread::print(std::ostream &os) {
    os << std::string(50, '#') << '\n';
    os << "# " << title() << '\n';
    os << "# dist_R_cutoff    (Ang) > " << dist_R_cutoff << '\n';
    os << "# angle_HOO_cutoff (Ang) > " << angle_HOO_cutoff << '\n';
    os << "# Ow_mask                > " << Ow_atom_mask << '\n';
    os << "# metal_mask             > " << center_Metal_atom_mask << '\n';
    os << "# cutoff          (Ang)  > " << std::sqrt(cutoff2) << '\n';
    os << std::string(50, '#') << '\n';

    os << boost::format("%15s %15s %15s\n") % "Frame" % "HBondClusterAtomNum" % "MaxOwDistance(Ang)";

    for (const auto &ele : hbond_connected_num_and_max_distance | boost::adaptors::indexed(1)) {
        os << boost::format("%15d %15.3f %15.3f\n") % ele.index() % (ele.value().first) % (ele.value().second);
    }

    os << std::string(50, '#') << '\n';
}

void HBondSpread::readInfo() {
    dist_R_cutoff = choose(0.0, 100.0, "Distance Cutoff(O-O) for Hydogen Bond [3.5 Ang] :", Default(3.5));
    angle_HOO_cutoff = choose(0.0, 100.0, "Angle Cutoff(H-O-O) for Hydogen Bond [30 degree] :", Default(30.0));
    auto cutoff = choose(0.0, 100.0, "cutoff :");
    cutoff2 = cutoff * cutoff;
    Atom::select1group(center_Metal_atom_mask, "Enter mask for Metal > ");
    Atom::select1group(Ow_atom_mask, "Enter mask for Ow > ");
}
