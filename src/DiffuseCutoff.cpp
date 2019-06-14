//
// Created by xiamr on 6/14/19.
//

#include "DiffuseCutoff.hpp"
#include "frame.hpp"


using namespace std;

bool operator==(const DiffuseCutoff::InnerAtom &i1, const DiffuseCutoff::InnerAtom &i2) {
    return i1.index == i2.index and i1.list_ptr == i2.list_ptr;
}

auto DiffuseCutoff::find_in(int seq) {
    for (auto iter = inner_atoms.begin(); iter != inner_atoms.end(); ++iter) {
        if (seq == iter->index) return iter;
    }
    return inner_atoms.end();
}

void DiffuseCutoff::process(std::shared_ptr<Frame> &frame) {


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
                    auto &old = it->list_ptr->back();

                    double xold = get<0>(old);
                    double yold = get<1>(old);
                    double zold = get<2>(old);

                    xr = x1 - xold;
                    yr = y1 - yold;
                    zr = z1 - zold;

                    frame->image(xr, yr, zr);

                    it->list_ptr->emplace_back(xr + xold, yr + yold, zr + zold);
                } else {
                    auto list_ptr = new std::list<std::tuple<double, double, double>>();
                    list_ptr->emplace_back(x1, y1, z1);
                    inner_atoms.insert(InnerAtom(mol->seq(), list_ptr));
                    rcm.emplace_back(list_ptr);
                }
            } else {
                if (it != inner_atoms.end()) {
                    inner_atoms.erase(it);
                }
            }

        }
    }
}

void DiffuseCutoff::print() {
    std::vector<std::pair<int, tuple<double, double, double>>> msd;
    for (auto list_ptr : this->rcm) {
        size_t i = 0;
        for (auto it1 = list_ptr->begin(); it1 != --list_ptr->end(); it1++) {
            i++;
            auto j = i;
            auto it2 = it1;
            for (it2++; it2 != list_ptr->end(); it2++) {
                j++;
                auto m = j - i - 1;
                if (m >= msd.size()) msd.emplace_back(0, std::make_tuple(0.0, 0.0, 0.0));
                double xdiff = get<0>(*it2) - get<0>(*it1);
                double ydiff = get<1>(*it2) - get<1>(*it1);
                double zdiff = get<2>(*it2) - get<2>(*it1);
                get<0>(msd[m].second) += xdiff * xdiff;
                get<1>(msd[m].second) += ydiff * ydiff;
                get<2>(msd[m].second) += zdiff * zdiff;
                msd[m].first++;
            }
        }
    }

    const double dunits = 10.0;

    for (auto &i : msd) {
        double counts = i.first;
        get<0>(i.second) /= counts;
        get<1>(i.second) /= counts;
        get<2>(i.second) /= counts;
    }

    outfile << "*********************************************************" << endl;
    outfile << "cutoff : " << std::sqrt(cutoff2) << std::endl;
    outfile << "First Type : " << ids1 << " Second Type : " << ids2 << endl;
    outfile << "Mean Squared Displacements and Self-Diffusion Constant" << endl;
    outfile << "    Time Gap      X MSD       Y MSD       Z MSD       R MSD       Diff Const" << endl;
    outfile << "      (ps)       (Ang^2)     (Ang^2)     (Ang^2)     (Ang^2)    (x 10^-5 cm**2/sec)" << endl;

    for (size_t i = 0; i < msd.size(); i++) {
        double delta = time_increment_ps * (i + 1);
        double xvalue = get<0>(msd[i].second);
        double yvalue = get<1>(msd[i].second);
        double zvalue = get<2>(msd[i].second);
        double rvalue = xvalue + yvalue + zvalue;
        double dvalue = dunits * rvalue / delta / 6.0;
        outfile << boost::format("%12.2f%12.2f%12.2f%12.2f%12.2f%12.4f\n") %
                   delta % xvalue % yvalue % zvalue % rvalue % dvalue;
    }
    outfile << "*********************************************************" << endl;
}

void DiffuseCutoff::readInfo() {

    Atom::select2group(ids1, ids2);

    double cutoff = choose(0.0, std::numeric_limits<double>::max(), "Please enter distance cutoff:");
    this->cutoff2 = cutoff * cutoff;

    while (true) {
        time_increment_ps = 0.1;
        string input_line = input(" Enter the Time Increment in Picoseconds [0.1]: ");
        boost::trim(input_line);
        if (!input_line.empty()) {
            time_increment_ps = stod(input_line);
            if (time_increment_ps <= 0.0) {
                cout << "error time increment " << time_increment_ps << endl;
                continue;
            }
        }
        break;
    }

}

void DiffuseCutoff::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->ids1)) this->group1.insert(atom);
                      if (Atom::is_match(atom, this->ids2)) this->group2.insert(atom);
                  });
}
