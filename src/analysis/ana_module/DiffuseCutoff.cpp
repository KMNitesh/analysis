
#include "DiffuseCutoff.hpp"

#include <boost/range/algorithm.hpp>

#include "data_structure/frame.hpp"
#include "data_structure/molecule.hpp"
#include "utils/common.hpp"

DiffuseCutoff::DiffuseCutoff() {
    enable_outfile = true;
    enable_forcefield = true;
}

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
                    auto list_ptr = new std::deque<std::tuple<double, double, double>>();
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
    std::vector<std::pair<int, std::tuple<double, double, double>>> msd;
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
                double xdiff = std::get<0>(*it2) - std::get<0>(*it1);
                double ydiff = std::get<1>(*it2) - std::get<1>(*it1);
                double zdiff = std::get<2>(*it2) - std::get<2>(*it1);
                std::get<0>(msd[m].second) += xdiff * xdiff;
                std::get<1>(msd[m].second) += ydiff * ydiff;
                std::get<2>(msd[m].second) += zdiff * zdiff;
                msd[m].first++;
            }
        }
    }

    const double dunits = 10.0;

    for (auto &i : msd) {
        double counts = i.first;
        std::get<0>(i.second) /= counts;
        std::get<1>(i.second) /= counts;
        std::get<2>(i.second) /= counts;
    }

    os << "*********************************************************\n";
    os << "# cutoff : " << std::sqrt(cutoff2) << '\n';
    os << "# First Type : " << mask1 << " Second Type : " << mask2 << '\n';
    os << "# Mean Squared Displacements and Self-Diffusion Constant\n";
    os << "#    Time Gap      X MSD       Y MSD       Z MSD       R MSD       Diff Const\n";
    os << "#      (ps)       (Ang^2)     (Ang^2)     (Ang^2)     (Ang^2)    (x 10^-5 cm**2/sec)\n";

    const boost::format fmt("%12.2f%12.2f%12.2f%12.2f%12.2f%12.4f\n");
    for (size_t i = 0; i < msd.size(); i++) {
        double delta = time_increment_ps * static_cast<double>(i + 1);
        double xvalue = std::get<0>(msd[i].second);
        double yvalue = std::get<1>(msd[i].second);
        double zvalue = std::get<2>(msd[i].second);
        double rvalue = xvalue + yvalue + zvalue;
        double dvalue = dunits * rvalue / delta / 6.0;
        os << boost::format(fmt) % delta % xvalue % yvalue % zvalue % rvalue % dvalue;
    }
    os << "*********************************************************\n";
}

void DiffuseCutoff::readInfo() {
    select2group(mask1, mask2);

    double cutoff = choose(0.0, std::numeric_limits<double>::max(), "Please enter distance cutoff:");
    this->cutoff2 = cutoff * cutoff;

    while (true) {
        time_increment_ps = 0.1;
        std::string input_line = input(" Enter the Time Increment in Picoseconds [0.1]: ");
        boost::trim(input_line);
        if (!input_line.empty()) {
            time_increment_ps = stod(input_line);
            if (time_increment_ps <= 0.0) {
                std::cout << "error time increment " << time_increment_ps << '\n';
                continue;
            }
        }
        break;
    }
}

void DiffuseCutoff::setParameters(const AmberMask &M, const AmberMask &L, double cutoff, double time_increment_ps,
                                  const std::string &outfilename) {
    this->mask1 = M;
    this->mask2 = L;

    if (cutoff <= 0) {
        throw std::runtime_error("`cutoff` must large than zero");
    }
    this->cutoff2 = cutoff * cutoff;

    if (time_increment_ps <= 0) {
        throw std::runtime_error("`time_increment_ps` must large than zero");
    }
    this->time_increment_ps = time_increment_ps;

    this->outfilename = outfilename;
    boost::trim(this->outfilename);
    if (this->outfilename.empty()) {
        throw std::runtime_error("outfilename cannot empty");
    }
}

void DiffuseCutoff::processFirstFrame(std::shared_ptr<Frame> &frame) {
    boost::for_each(frame->atom_list, [this](std::shared_ptr<Atom> &atom) {
        if (is_match(atom, this->mask1)) this->group1.insert(atom);
        if (is_match(atom, this->mask2)) this->group2.insert(atom);
    });
    if (group1.size() > 1) {
        std::cerr << "the reference(metal cation) atom for DiffuseCutoff function can only have one\n";
        exit(EXIT_FAILURE);
    }
}

std::string DiffuseCutoff::description() {
    std::stringstream ss;
    std::string title_line = "------ " + std::string(title()) + " ------";
    ss << title_line << "\n";
    ss << " M                 = [ " << mask1 << " ]\n";
    ss << " L                 = [ " << mask2 << " ]\n";
    ss << " cutoff            = " << sqrt(cutoff2) << " (Ang)\n";
    ss << " time_increment_ps = " << time_increment_ps << " (ps)\n";
    ss << " outfilename       = " << outfilename << "\n";
    ss << std::string(title_line.size(), '-') << '\n';
    return ss.str();
}

DiffuseCutoff::~DiffuseCutoff() {
    boost::for_each(rcm, std::default_delete<std::deque<std::tuple<double, double, double>>>());
}
