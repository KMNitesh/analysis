
#include "PBCUtils.hpp"

#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/range/numeric.hpp>

#include "data_structure/atom.hpp"
#include "data_structure/frame.hpp"
#include "data_structure/molecule.hpp"

std::shared_ptr<Atom> PBCUtils::find_atom(AmberMask &mask, std::shared_ptr<Frame> &frame) {
    std::shared_ptr<Atom> ret;
    for (auto &atom : frame->atom_list) {
        if (Atom::is_match(atom, mask)) {
            if (ret) {
                throw std::runtime_error("More then one atom selected");
            }
            ret = atom;
        }
    }
    return ret;
}

std::vector<std::shared_ptr<Atom>> PBCUtils::find_atoms(AmberMask &mask, std::shared_ptr<Frame> &frame) {
    std::vector<std::shared_ptr<Atom>> ret;
    for (auto &atom : frame->atom_list) {
        if (Atom::is_match(atom, mask)) {
            ret.push_back(atom);
        }
    }
    return ret;
}

std::shared_ptr<Molecule> PBCUtils::find_molecule(AmberMask &mask, std::shared_ptr<Frame> &frame) {
    std::shared_ptr<Molecule> ret;
    for (auto &atom : frame->atom_list) {
        if (Atom::is_match(atom, mask)) {
            if (ret and ret != atom->molecule.lock()) {
                throw std::runtime_error("More then one moleule selected");
            }
            ret = atom->molecule.lock();
        }
    }
    return ret;
}

void PBCUtils::do_move_center_basedto_atom(AmberMask &mask, std::shared_ptr<Frame> &frame) const {
    do_molecule_aggregate(frame);
    if (selected_atoms.empty()) selected_atoms.push_back(find_atom(mask, frame));
    auto center = selected_atoms.front()->getCoordinate();
    for (auto &mol : frame->molecule_list) {
        auto r = mol->calc_geom_center_inplace() - center;
        auto old = r;
        frame->image(r);
        r -= old;
        for (auto &atom : mol->atom_list) {
            std::tie(atom->x, atom->y, atom->z) += r;
        }
    }
    for (auto &atom : frame->atom_list) {
        std::tie(atom->x, atom->y, atom->z) -=
            center - 0.5 * (std::make_tuple(frame->box.box[0][0], frame->box.box[0][1], frame->box.box[0][2]) +
                            std::make_tuple(frame->box.box[1][0], frame->box.box[1][1], frame->box.box[1][2]) +
                            std::make_tuple(frame->box.box[2][0], frame->box.box[2][1], frame->box.box[2][2]));
    }
}

void PBCUtils::do_move_center_basedto_atom_group(AmberMask &mask, std::shared_ptr<Frame> &frame) const {
    if (selected_atoms.empty()) selected_atoms = find_atoms(mask, frame);
    auto mols = calculate_intermol(selected_atoms, frame);
    move(mols, frame);
    auto center = boost::accumulate(selected_atoms, std::tuple<double, double, double>{},
                                    [](auto sum, auto &atom) { return sum + atom->getCoordinate(); }) /
                  selected_atoms.size();

    for (auto &mol : frame->molecule_list) {
        if (mols.first.contains(mol)) continue;
        mol->do_aggregate(frame);
        auto r = mol->calc_geom_center_inplace() - center;
        auto old = r;
        frame->image(r);
        r -= old;
        for (auto &atom : mol->atom_list) {
            std::tie(atom->x, atom->y, atom->z) += r;
        }
    }

    for (auto &atom : frame->atom_list) {
        std::tie(atom->x, atom->y, atom->z) -=
            center - 0.5 * (std::make_tuple(frame->box.box[0][0], frame->box.box[0][1], frame->box.box[0][2]) +
                            std::make_tuple(frame->box.box[1][0], frame->box.box[1][1], frame->box.box[1][2]) +
                            std::make_tuple(frame->box.box[2][0], frame->box.box[2][1], frame->box.box[2][2]));
    }
}

void PBCUtils::do_move_center_basedto_molecule(AmberMask &mask, std::shared_ptr<Frame> &frame) const {
    do_molecule_aggregate(frame);
    if (!selected_molecule) selected_molecule = find_molecule(mask, frame);
    auto center = selected_molecule->calc_geom_center_inplace();
    for (auto &mol : frame->molecule_list) {
        auto r = mol->calc_geom_center_inplace() - center;
        auto old = r;
        frame->image(r);
        r -= old;
        for (auto &atom : mol->atom_list) {
            std::tie(atom->x, atom->y, atom->z) += r;
        }
    }
    for (auto &atom : frame->atom_list) {
        std::tie(atom->x, atom->y, atom->z) -=
            center - 0.5 * (std::make_tuple(frame->box.box[0][0], frame->box.box[0][1], frame->box.box[0][2]) +
                            std::make_tuple(frame->box.box[1][0], frame->box.box[1][1], frame->box.box[1][2]) +
                            std::make_tuple(frame->box.box[2][0], frame->box.box[2][1], frame->box.box[2][2]));
    }
}

void PBCUtils::do_molecule_aggregate(std::shared_ptr<Frame> &frame) const {
    for (auto &mol : frame->molecule_list) {
        mol->do_aggregate(frame);
    }
}

void PBCUtils::do_move_all_atom_into_box(std::shared_ptr<Frame> &frame) const {
    for (auto &atom : frame->atom_list) {
        frame->image(atom->x, atom->y, atom->z);
    }
}

void PBCUtils::doPBC(Trajconv::PBCType pbc_mode, AmberMask &mask, std::shared_ptr<Frame> &frame) const {
    if (pbc_mode == Trajconv::PBCType::AllIntoBox) {
        do_move_all_atom_into_box(frame);
    } else if (pbc_mode == Trajconv::PBCType::OneAtom) {
        do_move_center_basedto_atom(mask, frame);
    } else if (pbc_mode == Trajconv::PBCType::OneMol) {
        do_move_center_basedto_molecule(mask, frame);
    } else if (pbc_mode == Trajconv::PBCType::AtomGroup) {
        do_move_center_basedto_atom_group(mask, frame);
    }
}

void PBCUtils::move(std::set<std::shared_ptr<Molecule>> &mols_set,
                    std::vector<std::pair<std::shared_ptr<Molecule>, std::shared_ptr<Molecule>>> &mols_seq,
                    const std::shared_ptr<Frame> &frame) {
    for (auto &mol : mols_set) {
        mol->do_aggregate(frame);
        mol->calc_geom_center_inplace();
    }
    for (auto &[source, target] : mols_seq) {
        auto r = target->inplace_geometry_center - source->inplace_geometry_center;
        auto old = r;
        frame->image(r);
        r -= old;
        for (auto &atom : target->atom_list) {
            std::tie(atom->x, atom->y, atom->z) += r;
        }
        target->inplace_geometry_center += r;
    }
}

PBCUtils::MolPair PBCUtils::calculate_intermol_imp(const std::set<std::shared_ptr<Molecule>> &mols_set,
                                                   const std::shared_ptr<Frame> &frame) {
    for (auto &mol : mols_set) {
        mol->do_aggregate(frame);
        mol->calc_geom_center_inplace();
    }
    std::vector<std::shared_ptr<Molecule>> mols(std::begin(mols_set), std::end(mols_set));
    Graph g;
    for (auto &mol : mols) mol->vertex_descriptor = boost::add_vertex(g);

    using namespace boost;
    for (auto it1 = mols.begin(); it1 != --mols.end(); ++it1) {
        auto center1 = (*it1)->inplace_geometry_center;
        auto it2 = it1;
        for (++it2; it2 != mols.end(); ++it2) {
            auto center2 = (*it2)->inplace_geometry_center;
            auto r = center2 - center1;
            frame->image(r);
            auto d2 = vector_norm2(r);
            add_edge((*it1)->vertex_descriptor, (*it2)->vertex_descriptor, d2, g);
        }
    }

    std::vector<graph_traits<Graph>::vertex_descriptor> p(num_vertices(g));
    prim_minimum_spanning_tree(g, p.data());

    std::vector<std::pair<std::size_t, std::size_t>> edges;
    std::size_t root = 0;
    for (std::size_t i = 0; i != p.size(); ++i)
        if (p[i] != i)
            edges.emplace_back(p[i], i);
        else
            root = i;

    adjacency_list<vecS, vecS, directedS> tg(edges.begin(), edges.end(), p.size());
    std::vector<std::pair<std::shared_ptr<Molecule>, std::shared_ptr<Molecule>>> mole_seq;
    Visitor vis(mols, mole_seq);
    breadth_first_search(tg, root, boost::visitor(vis));
    return {mols_set, mole_seq};
}
