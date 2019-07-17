//
// Created by xiamr on 6/14/19.
//

#include "RotAcfCutoff.hpp"

#include "frame.hpp"
#include "atom.hpp"
#include "VectorSelectorFactory.hpp"

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

    if (!ref) {
        std::cerr << "reference atom not found" << std::endl;
        exit(5);
    }
    double ref_x = ref->x;
    double ref_y = ref->y;
    double ref_z = ref->z;
    for (auto &atom2: group2) {
        auto mol = atom2->molecule.lock();
        auto coord = atom2->getCoordinate();
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

void RotAcfCutoff::readInfo() {

    Atom::select2group(ids1, ids2);

    vectorSelector = VectorSelectorFactory::getVectorSelector();
    vectorSelector->readInfo();

    std::cout << "Legendre Polynomial\n";
    std::cout << "1. P1 = x\n";
    std::cout << "2. P2 = (1/2)(3x^2 -1)\n";
    LegendrePolynomial = choose(1, 2, "select > ");
    double cutoff = choose(0.0, std::numeric_limits<double>::max(), "Please enter distance cutoff:");
    this->cutoff2 = cutoff * cutoff;

    this->time_increment_ps = choose(0.0, std::numeric_limits<double>::max(),
                                     "Enter the Time Increment in Picoseconds [0.1]:", true, 0.1);
    this->max_time_grap = choose(0.0, std::numeric_limits<double>::max(),
                                 "Enter the Max Time Grap in Picoseconds :");

}

void RotAcfCutoff::setParameters(const Atom::Node &M, const Atom::Node &L, std::shared_ptr<VectorSelector> vector,
                                 int LegendrePolynomial, double cutoff, double time_increment_ps,
                                 double max_time_grap_ps, const std::string &outfilename) {
    this->ids1 = M;
    this->ids2 = L;

    if (!vector) {
        throw runtime_error("vector not vaild");
    } else {
        this->vectorSelector = vector;
    }
    if (!(LegendrePolynomial == 1 or LegendrePolynomial == 2)) {
        throw runtime_error("legendre polynomial must be 1 or 2");
    }

    this->LegendrePolynomial = LegendrePolynomial;

    if (cutoff <= 0) {
        throw runtime_error("`cutoff` must be postive");
    } else {
        cutoff2 = cutoff * cutoff;
    }

    if (time_increment_ps <= 0) {
        throw runtime_error("`time_increment_ps` must be postive");
    } else {
        this->time_increment_ps = time_increment_ps;
    }

    if (max_time_grap_ps <= 0) {
        throw runtime_error("`max_time_grap_ps` must be postive");
    } else if (max_time_grap_ps <= this->time_increment_ps) {
        throw runtime_error("`max_time_grap_ps` must be larger than `time_increment_ps`");
    } else {
        this->max_time_grap = max_time_grap_ps;
    }
    this->outfilename = outfilename;
    boost::trim(this->outfilename);
    if (this->outfilename.empty()) {
        throw runtime_error("outfilename cannot empty");
    }
}


tuple<double, double, double> RotAcfCutoff::calVector(shared_ptr<Molecule> &mol, shared_ptr<Frame> &frame) {
    return vectorSelector->calculateVector(mol, frame);
}

void RotAcfCutoff::print(std::ostream &os) {
    std::vector<std::pair<unsigned long long, double>> acf;
    acf.emplace_back(0, 0.0);
    switch (LegendrePolynomial) {
        case 1:
            calculateAutocorrelaionFunction(acf, [](auto x) { return x; });
            break;
        case 2:
            calculateAutocorrelaionFunction(acf, [](auto x) { return 0.5 * (3 * x * x - 1); });
            break;
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
    vectorSelector->print(os);
    os << "cutoff : " << std::sqrt(cutoff2) << std::endl;
    os << "First Type : " << ids1 << " Second Type : " << ids2 << endl;

    os << " rotational autocorrelation function" << endl;
    os << "Legendre Polynomial : ";
    switch (LegendrePolynomial) {
        case 1:
            os << "P1 = x\n";
            break;
        case 2:
            os << "P2 = (1/2)(3x^2 -1)\n";
            break;
    }
    os << "    Time Gap      ACF       integrate" << endl;
    os << "      (ps)                    (ps)" << endl;

    for (std::size_t t = 0; t < acf.size(); t++) {
        os << boost::format("%12.2f%18.14f%15.5f") % (t * time_increment_ps) % acf[t].second % integrate[t]
           << endl;
    }
    os << "*********************************************************" << endl;

}

template<typename Function>
void RotAcfCutoff::calculateAutocorrelaionFunction(vector<pair<unsigned long long, double>> &acf, Function f) const {
    size_t max_time_grap_step = ceil(max_time_grap / time_increment_ps);
    for (auto list_ptr : rots) {
        size_t i = 0;
        for (auto it1 = list_ptr->begin(); it1 != --list_ptr->end(); it1++) {
            i++;
            auto j = i;
            auto it2 = it1;
            for (it2++; it2 != list_ptr->end(); it2++) {
                j++;
                auto m = j - i;
                if (m > max_time_grap_step) break;
                if (m >= acf.size()) acf.emplace_back(0, 0.0);

                double cos = dot_multiplication(*it1, *it2);

                acf[m].second += f(cos);
                acf[m].first++;
            }
        }
    }
}

void RotAcfCutoff::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->ids1)) this->group1.insert(atom);
                      if (Atom::is_match(atom, this->ids2)) this->group2.insert(atom);
                  });
}

