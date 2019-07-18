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
    for (auto &ref : group1) {
        for (auto &atom2 : group2) {
            auto it = find_in(atom2->molecule.lock()->seq());
            if (atom_distance2(ref, atom2, frame) < cutoff2) {
                // in the shell
                auto coord = atom2->molecule.lock()->calc_weigh_center(frame);
                if (it != inner_atoms.end()) {
                    auto &old = it->list_ptr->back();
                    auto shift = coord - old;
                    frame->image(shift);
                    it->list_ptr->push_back(shift + old);
                } else {
                    auto list_ptr = new std::list<std::tuple<double, double, double>>();
                    list_ptr->push_back(coord);
                    inner_atoms.insert(InnerAtom(atom2->molecule.lock()->seq(), list_ptr));
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

void DiffuseCutoff::print(std::ostream &os) {
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

    os << "*********************************************************" << endl;
    os << "cutoff : " << std::sqrt(cutoff2) << std::endl;
    os << "First Type : " << ids1 << " Second Type : " << ids2 << endl;
    os << "Mean Squared Displacements and Self-Diffusion Constant" << endl;
    os << "    Time Gap      X MSD       Y MSD       Z MSD       R MSD       Diff Const" << endl;
    os << "      (ps)       (Ang^2)     (Ang^2)     (Ang^2)     (Ang^2)    (x 10^-5 cm**2/sec)" << endl;

    for (size_t i = 0; i < msd.size(); i++) {
        double delta = time_increment_ps * (i + 1);
        double xvalue = get<0>(msd[i].second);
        double yvalue = get<1>(msd[i].second);
        double zvalue = get<2>(msd[i].second);
        double rvalue = xvalue + yvalue + zvalue;
        double dvalue = dunits * rvalue / delta / 6.0;
        os << boost::format("%12.2f%12.2f%12.2f%12.2f%12.2f%12.4f\n") %
              delta % xvalue % yvalue % zvalue % rvalue % dvalue;
    }
    os << "*********************************************************" << endl;
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

void DiffuseCutoff::setParameters(const Atom::Node &M, const Atom::Node &L,
                                  double cutoff, double time_increment_ps, const std::string &outfilename) {
    this->ids1 = M;
    this->ids2 = L;

    if (cutoff <= 0) {
        throw runtime_error("`cutoff` must large than zero");
    }
    this->cutoff2 = cutoff * cutoff;

    if (time_increment_ps <= 0) {
        throw runtime_error("`time_increment_ps` must large than zero");
    }
    this->time_increment_ps = time_increment_ps;

    this->outfilename = outfilename;
    boost::trim(this->outfilename);
    if (this->outfilename.empty()) {
        throw runtime_error("outfilename cannot empty");
    }
}

void DiffuseCutoff::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->ids1)) this->group1.insert(atom);
                      if (Atom::is_match(atom, this->ids2)) this->group2.insert(atom);
                  });
    if (group1.size() > 1) {
        cerr << "the reference(metal cation) atom for DiffuseCutoff function can only have one\n";
        exit(EXIT_FAILURE);
    }
}
