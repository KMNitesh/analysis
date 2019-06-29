//
// Created by xiamr on 6/14/19.
//

#include "RotAcfCutoff.hpp"

#include "frame.hpp"
#include "atom.hpp"

using namespace std;

bool operator==(const RotAcfCutoff::InnerAtom &i1, const RotAcfCutoff::InnerAtom &i2) {
    return i1.index == i2.index and i1.list_ptr == i2.list_ptr;
}

auto RotAcfCutoff::find_in(int seq) {
    for (auto iter = inner_atoms.begin(); iter != inner_atoms.end(); ++iter) {
        if (seq == iter->index) return iter;
    }
    return inner_atoms.end();
}

void RotAcfCutoff::process(std::shared_ptr<Frame> &frame) {

    shared_ptr<Atom> ref;

    for (auto &mol : frame->molecule_list) {
        mol->bExculde = true;
    }

    for (auto &atom : group1) {
        ref = atom;
        if (!atom->molecule.expired()) atom->molecule.lock()->calc_mass();
    }

    for (auto &atom : group2) {
        if (!atom->molecule.expired()) {
            auto mol = atom->molecule.lock();

            mol->calc_mass();
            mol->bExculde = false;
        }
    }

    if (!ref) {
        std::cerr << "reference atom not found" << std::endl;
        exit(5);
    }
    double ref_x = ref->x;
    double ref_y = ref->y;
    double ref_z = ref->z;
    for (auto &mol: frame->molecule_list) {
        if (!mol->bExculde) {
            auto coord = mol->calc_weigh_center(frame);
            double x1 = get<0>(coord);
            double y1 = get<1>(coord);
            double z1 = get<2>(coord);

            double xr = x1 - ref_x;
            double yr = y1 - ref_y;
            double zr = z1 - ref_z;

            frame->image(xr, yr, zr);

            auto it = find_in(mol->seq());
            if (xr * xr + yr * yr + zr * zr < cutoff2) {
                // in the shell
                if (it != inner_atoms.end()) {
                    it->list_ptr->push_back(calVector(mol, frame));
                } else {
                    auto list_ptr = new std::list<std::tuple<double, double, double>>();
                    list_ptr->push_back(calVector(mol, frame));
                    inner_atoms.insert(InnerAtom(mol->seq(), list_ptr));
                    rots.emplace_back(list_ptr);
                }
            } else {
                if (it != inner_atoms.end()) {
                    inner_atoms.erase(it);
                }
            }
        }
    }
}

void RotAcfCutoff::readInfo() {

    Atom::select2group(ids1, ids2);

    double cutoff = choose(0.0, std::numeric_limits<double>::max(), "Please enter distance cutoff:");
    this->cutoff2 = cutoff * cutoff;

    this->time_increment_ps = choose(0.0, std::numeric_limits<double>::max(),
                                     "Enter the Time Increment in Picoseconds [0.1]:", true, 0.1);

}

tuple<double, double, double> RotAcfCutoff::calVector(shared_ptr<Molecule> &mol, shared_ptr<Frame> &frame) {
    auto iter = mol->atom_list.begin();

    auto atom_i = *iter;
    iter++;
    auto atom_j = *iter;
    iter++;
    auto atom_k = *iter;

    auto u1 = atom_i->x - atom_j->x;
    auto u2 = atom_i->y - atom_j->y;
    auto u3 = atom_i->z - atom_j->z;

    auto v1 = atom_k->x - atom_j->x;
    auto v2 = atom_k->y - atom_j->y;
    auto v3 = atom_k->z - atom_j->z;

    frame->image(u1, u2, u3);
    frame->image(v1, v2, v3);

    auto xv3 = u2 * v3 - u3 * v2;
    auto yv3 = u3 * v1 - u1 * v3;
    auto zv3 = u1 * v2 - u2 * v1;

    return std::make_tuple(xv3, yv3, zv3);
}

void RotAcfCutoff::print(std::ostream &os) {
    std::vector<std::pair<int, double>> acf;
    acf.emplace_back(0, 0.0);
    for (auto list_ptr : this->rots) {
        size_t i = 0;
        for (auto it1 = list_ptr->begin(); it1 != --list_ptr->end(); it1++) {
            i++;
            auto j = i;
            auto it2 = it1;
            for (it2++; it2 != list_ptr->end(); it2++) {
                j++;
                auto m = j - i;
                if (m >= acf.size()) acf.emplace_back(0, 0.0);

                double xr1 = get<0>(*it1);
                double yr1 = get<1>(*it1);
                double zr1 = get<2>(*it1);

                double xr2 = get<0>(*it2);
                double yr2 = get<1>(*it2);
                double zr2 = get<2>(*it2);


                double r1_2 = xr1 * xr1 + yr1 * yr1 + zr1 * zr1;
                double r2_2 = xr2 * xr2 + yr2 * yr2 + zr2 * zr2;

                double dot = xr1 * xr2 + yr1 * yr2 + zr1 * zr2;

                double cos = dot / std::sqrt(r1_2 * r2_2);

                acf[m].second += cos;
                acf[m].first++;
            }
        }
    }

    for (size_t i = 1; i < acf.size(); i++) {
        double counts = acf[i].first;
        acf[i].second /= counts;
    }

    acf[0].second = 1.0;
    // intergrate;

    std::vector<double> integrate(acf.size());
    integrate[0] = 0.0;

    for (std::size_t i = 1; i < integrate.size(); i++) {
        integrate[i] = integrate[i - 1] + 0.5 * (acf[i - 1].second + acf[i].second) * time_increment_ps;
    }

    os << "*********************************************************" << endl;
    os << "cutoff : " << std::sqrt(cutoff2) << std::endl;
    os << "First Type : " << ids1 << " Second Type : " << ids2 << endl;
    os << " rotational autocorrelation function" << endl;

    os << "    Time Gap      ACF       intergrate" << endl;
    os << "      (ps)                    (ps)" << endl;

    for (std::size_t t = 0; t < acf.size(); t++) {
        os << boost::format("%12.2f%18.14f%15.5f") % (t * time_increment_ps) % acf[t].second % integrate[t]
           << endl;
    }
    os << "*********************************************************" << endl;

}

void RotAcfCutoff::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->ids1)) this->group1.insert(atom);
                      if (Atom::is_match(atom, this->ids2)) this->group2.insert(atom);
                  });
}
