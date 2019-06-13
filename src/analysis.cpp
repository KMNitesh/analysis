#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <unordered_set>
#include <memory>
#include <cmath>
#include <tbb/tbb.h>
#include <readline/readline.h>

#include <alloca.h>
#include <stdexcept>
#include <algorithm>
#include <cstdlib>
#include <fmt/format.h>
#include <fmt/printf.h>
#include <string>
#include <unistd.h>
#include <functional>

#include <Eigen/Eigen>

#define BOOST_RESULT_OF_USE_DECLTYPE
#define BOOST_SPIRIT_USE_PHOENIX_V3

#define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
#define BOOST_MPL_LIMIT_VECTOR_SIZE 30


#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>

#include <boost/filesystem.hpp>

// Boost metaprogramming library
#include <boost/mpl/vector.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/unique.hpp>
#include <boost/mpl/string.hpp>


namespace mpl = boost::mpl;

#include <boost/format.hpp>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/phoenix.hpp>
#include <boost/variant.hpp>
#include <boost/optional.hpp>
#include <boost/fusion/sequence/intrinsic/at_c.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/phoenix/function/adapt_function.hpp>

#include <boost/assert.hpp>

#include <boost/bimap.hpp>
#include <boost/bimap/unordered_set_of.hpp>
#include <boost/bimap/set_of.hpp>
#include <boost/bimap/list_of.hpp>
#include <boost/assign.hpp>

#include "common.hpp"
#include "grammar.hpp"
#include "molecule.hpp"
#include "frame.hpp"
#include "trajectoryreader.hpp"
#include "forcefield.hpp"
#include "gro_writer.hpp"
#include "trr_writer.hpp"
#include "xtc_writer.hpp"
#include "netcdf_writer.hpp"
#include "center_selection_grammar.hpp"


using namespace std;

constexpr int ATOM_MAX = 10000;

std::fstream outfile;

using boost::variant;


namespace qi = boost::spirit::qi;
namespace fusion = boost::fusion;
namespace phoenix = boost::phoenix;


class BasicAnalysis {
public:
    virtual void processFirstFrame(std::shared_ptr<Frame> &frame) {};

    virtual void process(std::shared_ptr<Frame> &frame) = 0;

    virtual void print() = 0;

    virtual void readInfo() = 0;

    static const string title() { return "Base Class"; }


    virtual ~BasicAnalysis() = default;
};


class GmxTrj : public BasicAnalysis {
    enum class PBCType {
        None,
        OneAtom,
        OneMol
    } pbc_type;
    int step = 0;
    string grofilename;
    string xtcfilename;
    string trrfilename;
    string mdcrdfilename;
    XTCWriter xtc;
    TRRWriter trr;
    NetCDFWriter mdcrd;

    int num;

    bool enable_xtc = true;
    bool enable_trr = true;
    bool enable_gro = true;
    bool enable_mdcrd = true;
public:
    void process(std::shared_ptr<Frame> &frame) override;

    void print() override;

    void readInfo() override;

    static const string title() {
        return "Gromacs XTC & TRR & GRO & NetCDF Output";
    }


};

void GmxTrj::process(std::shared_ptr<Frame> &frame) {
    if (pbc_type == PBCType::OneAtom) {
        auto center_x = frame->atom_map[this->num]->x;
        auto center_y = frame->atom_map[this->num]->y;
        auto center_z = frame->atom_map[this->num]->z;

        for (auto &mol : frame->molecule_list) {
            for (auto &atom : mol->atom_list) {
                atom->x -= center_x;
                atom->y -= center_y;
                atom->z -= center_z;
            }
        }
    } else if (pbc_type == PBCType::OneMol) {
        auto center_mol = frame->atom_map[this->num]->molecule.lock();
        center_mol->calc_center(frame);
        auto center_x = center_mol->center_x;
        auto center_y = center_mol->center_y;
        auto center_z = center_mol->center_z;
        for (auto &mol : frame->molecule_list) {
            for (auto &atom : mol->atom_list) {
                atom->x -= center_x;
                atom->y -= center_y;
                atom->z -= center_z;
            }
        }
    }
    if (pbc_type != PBCType::None) {
        for (auto &mol : frame->molecule_list) {
            mol->calc_center(frame);
            double x_move = 0.0;
            double y_move = 0.0;
            double z_move = 0.0;
            while (mol->center_x > frame->a_axis_half) {
                mol->center_x -= frame->a_axis;
                x_move -= frame->a_axis;
            }
            while (mol->center_x < -frame->a_axis_half) {
                mol->center_x += frame->a_axis;
                x_move += frame->a_axis;
            }
            while (mol->center_y > frame->b_axis_half) {
                mol->center_y -= frame->b_axis;
                y_move -= frame->b_axis;
            }
            while (mol->center_y < -frame->b_axis_half) {
                mol->center_y += frame->b_axis;
                y_move += frame->b_axis;
            }
            while (mol->center_z > frame->c_axis_half) {
                mol->center_z -= frame->c_axis;
                z_move -= frame->c_axis;
            }
            while (mol->center_z < -frame->c_axis_half) {
                mol->center_z += frame->c_axis;
                z_move += frame->c_axis;
            }
            for (auto &atom : mol->atom_list) {
                atom->x += x_move + frame->a_axis_half;
                atom->y += y_move + frame->b_axis_half;
                atom->z += z_move + frame->c_axis_half;
            }
        }
    }

    if (step == 0) {
        GROWriter w;
        w.write(grofilename, frame);
        if (enable_xtc) xtc.open(xtcfilename);
        if (enable_trr) trr.open(trrfilename);
        if (enable_mdcrd) mdcrd.open(mdcrdfilename, frame->atom_list.size());
    }
    if (enable_xtc) xtc.write(frame);
    if (enable_trr) trr.write(frame);
    if (enable_mdcrd) mdcrd.write(frame);
    step++;
}

void GmxTrj::print() {
    if (enable_xtc) xtc.close();
    if (enable_trr) trr.close();
    if (enable_mdcrd) mdcrd.close();
}

void GmxTrj::readInfo() {
    string line = input("output gro file [Y]: ");
    boost::trim(line);
    if (!line.empty()) {
        if (line[0] == 'N' or line[0] == 'n') {
            enable_gro = false;
        }
    }
    if (enable_gro) {
        grofilename = input("GRO filename : ");
        boost::trim(grofilename);
    }
    line = input("output xtc file [Y]: ");
    boost::trim(line);
    if (!line.empty()) {
        if (line[0] == 'N' or line[0] == 'n') {
            enable_xtc = false;
        }
    }
    if (enable_xtc) {
        xtcfilename = input("XTC filename : ");
        boost::trim(xtcfilename);
    }

    line = input("output trr file [Y]: ");
    boost::trim(line);
    if (!line.empty()) {
        if (line[0] == 'N' or line[0] == 'n') {
            enable_trr = false;
        }
    }
    if (enable_trr) {
        trrfilename = input("TRR filename : ");
        boost::trim(trrfilename);
    }

    line = input("output amber nc file [Y]: ");
    boost::trim(line);
    if (!line.empty()) {
        if (line[0] == 'N' or line[0] == 'n') {
            enable_mdcrd = false;
        }
    }
    if (enable_mdcrd) {
        mdcrdfilename = input("nc filename : ");
        boost::trim(mdcrdfilename);
    }
    cout << "PBC transform option" << endl;
    cout << "(0) Do nothing" << endl;
    cout << "(1) Make atom i as center" << endl;
    cout << "(2) Make molecule i as center" << endl;

    int option = 0;
    while (true) {
        line = input("Choose:");
        boost::trim(line);

        if (!line.empty()) {

            try {
                option = std::stoi(line);
                string line2;
                switch (option) {
                    case 0:
                        pbc_type = PBCType::None;
                        break;
                    case 1:
                        pbc_type = PBCType::OneAtom;
                        while (true) {
                            line2 = input("Plese enter the atom NO. : ");
                            try {
                                num = std::stoi(line2);
                            } catch (std::invalid_argument &e) {
                                cerr << "must be a number !" << e.what() << endl;
                                continue;
                            }


                            if (num <= 0) {
                                cerr << "the atom number must be positive" << endl;
                                continue;
                            }
                            break;
                        }
                        break;
                    case 2:
                        pbc_type = PBCType::OneMol;
                        while (true) {
                            line2 = input("Plese enter one atom NO. that the molecule include:");
                            try {
                                num = std::stoi(line2);
                            } catch (std::invalid_argument &e) {
                                cerr << "must be a number !" << e.what() << endl;
                                continue;
                            }
                            if (num <= 0) {
                                cerr << "the atom number must be positive" << endl;
                                continue;
                            }
                            break;
                        }
                        break;
                    default:
                        cerr << "option not found !" << endl;
                        continue;
                }
                break;

            } catch (std::invalid_argument &e) {
                cerr << "must be a number !" << e.what() << endl;
            }
        }
    }


}


// Distance
class Distance : public BasicAnalysis {
    std::list<double> group_dist_list;

    Atom::AtomIndenter ids1;
    Atom::AtomIndenter ids2;

    std::unordered_set<shared_ptr<Atom>> group1;
    std::unordered_set<shared_ptr<Atom>> group2;


public:
    Distance() {
        enable_outfile = true;
    }

    void process(std::shared_ptr<Frame> &frame) override;

    void print() override;

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void readInfo() override;

    static const string title() {
        return "Distance";
    }

};

void Distance::process(std::shared_ptr<Frame> &frame) {

    double x1, y1, z1, x2, y2, z2;
    x1 = y1 = z1 = x2 = y2 = z2 = 0.0;
    double weigh1, weigh2;
    weigh1 = weigh2 = 0.0;
    for (auto &atom1 : group1) {
        double mass = forcefield.find_mass(atom1);
        x1 += atom1->x * mass;
        y1 += atom1->y * mass;
        z1 += atom1->z * mass;
        weigh1 += mass;
    }
    x1 /= weigh1;
    y1 /= weigh1;
    z1 /= weigh1;

    for (auto &atom2 : group2) {
        double mass = forcefield.find_mass(atom2);
        x2 += atom2->x * mass;
        y2 += atom2->y * mass;
        z2 += atom2->z * mass;
        weigh2 += mass;
    }
    x2 /= weigh2;
    y2 /= weigh2;
    z2 /= weigh2;

    double xr = x2 - x1;
    double yr = y2 - y1;
    double zr = z2 - z1;

    frame->image(xr, yr, zr);

    double dist = sqrt(xr * xr + yr * yr + zr * zr);
    group_dist_list.push_back(dist);
}


void Distance::print() {
    outfile << "****************************************" << endl;
    outfile << "Distance between " << endl;
    outfile << "group 1 :" << ids1 << std::endl;
    outfile << "group 2 :" << ids2 << std::endl;
    outfile << "****************************************" << endl;

    int cyc = 1;
    for (auto &dist : group_dist_list) {
        outfile << cyc << "    " << dist << std::endl;
        cyc++;
    }
}

void Distance::readInfo() {
    enable_forcefield = true;
    Atom::select2group(ids1, ids2);
}

void Distance::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->ids1)) this->group1.insert(atom);
                      if (Atom::is_match(atom, this->ids2)) this->group2.insert(atom);
                  });
}
// End of Distance


class CoordinateNumPerFrame : public BasicAnalysis {

    Atom::AtomIndenter ids1;
    Atom::AtomIndenter ids2;

    std::unordered_set<shared_ptr<Atom>> group1;
    std::unordered_set<shared_ptr<Atom>> group2;


    double dist_cutoff;
    list<int> cn_list;
public:
    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

public:
    CoordinateNumPerFrame() {
        enable_outfile = true;
    }

    void process(std::shared_ptr<Frame> &frame) override;

    void print() override;

    void readInfo() override;

    static const string title() {
        return "Coordinate Number per Frame";
    }
};

void CoordinateNumPerFrame::process(std::shared_ptr<Frame> &frame) {

    int cn_sum = 0;
    for (auto &atom1 : group1) {
        for (auto &atom2 : group2) {
            if (atom_distance(atom1, atom2, frame) <= this->dist_cutoff) {
                cn_sum++;
            }
        }
    }
    cn_list.push_back(cn_sum);

}

void CoordinateNumPerFrame::print() {
    outfile << "***************************" << endl;
    outfile << "******* CN per Frame ******" << endl;
    outfile << "type 1 :" << ids1 << endl;
    outfile << "type 2 :" << ids2 << endl;
    outfile << "cutoff :" << dist_cutoff << endl;
    outfile << "***************************" << endl;
    int cyc = 1;
    for (auto &cn : cn_list) {
        outfile << cyc << "   " << cn << endl;
        cyc++;
    }
    outfile << "***************************" << endl;
}

void CoordinateNumPerFrame::readInfo() {
    Atom::select2group(ids1, ids2);
    dist_cutoff = choose(0.0, GMX_DOUBLE_MAX, "Please enter distance cutoff:");
}

void CoordinateNumPerFrame::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->ids1)) this->group1.insert(atom);
                      if (Atom::is_match(atom, this->ids2)) this->group2.insert(atom);
                  });
}


class FirstCoordExchangeSearch : public BasicAnalysis {

public:
    FirstCoordExchangeSearch() {
        enable_outfile = true;
    }

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print() override;

    void readInfo() override;

    static const string title() {
        return "Water exchange analysis";
    }


private:

    Atom::AtomIndenter ids1;
    Atom::AtomIndenter ids2;

    std::unordered_set<shared_ptr<Atom>> group1;
    std::unordered_set<shared_ptr<Atom>> group2;

    double dist_cutoff, tol_dist, time_cutoff;
    enum class Direction {
        IN, OUT
    };
    typedef struct {
        Direction direction;
        int seq;
        int exchange_frame;
    } ExchangeItem;

    list<ExchangeItem> exchange_list;

    int step = 0;

    typedef struct {
        bool inner;
    } State;
    map<int, State> state_machine;

    std::set<int> init_seq_in_shell;

};

void FirstCoordExchangeSearch::process(std::shared_ptr<Frame> &frame) {

    step++;

    for (auto &atom1 : group1) {
        for (auto &atom2 : group2) {
            if (step == 1) {
                State state;
                state.inner = atom_distance(atom1, atom2, frame) <= this->dist_cutoff;
                state_machine[atom2->seq] = state;
                if (state.inner) init_seq_in_shell.insert(atom2->seq);
            } else {
                auto &state = state_machine[atom2->seq];
                if (state.inner) {
                    if (atom_distance(atom1, atom2, frame) >= this->dist_cutoff + tol_dist) {
                        state.inner = false;
                        ExchangeItem item;
                        item.seq = atom2->seq;
                        item.direction = Direction::OUT;
                        item.exchange_frame = step;
                        exchange_list.push_back(item);
                    }
                } else {
                    if (atom_distance(atom1, atom2, frame) <= this->dist_cutoff - tol_dist) {
                        state.inner = true;
                        ExchangeItem item;
                        item.seq = atom2->seq;
                        item.direction = Direction::IN;
                        item.exchange_frame = step;
                        exchange_list.push_back(item);
                    }
                }
            }
        }

    }
}

void FirstCoordExchangeSearch::print() {
    outfile << "***************************" << endl;
    outfile << "***** Exchange Search *****" << endl;
    outfile << "type 1 :" << ids1 << endl;
    outfile << "type 2 :" << ids2 << endl;

    outfile << "cutoff :" << dist_cutoff << endl;
    outfile << "tol dist :" << tol_dist << endl;
    outfile << "time_cutoff :" << time_cutoff << endl;
    outfile << "***************************" << endl;
    for (int seq : init_seq_in_shell) {
        outfile << "  " << seq;
    }
    outfile << "\n***************************" << endl;

    outfile << "** seq **  direction ***** exchange frame *****" << endl;
    for (auto it = exchange_list.begin(); it != exchange_list.end(); it++) {
        outfile << boost::format("%10d%6s%15d   !   ")
                   % it->seq
                   % (it->direction == Direction::IN ? "IN" : "OUT")
                   % it->exchange_frame;
        if (it->direction == Direction::IN) {
            init_seq_in_shell.insert(it->seq);
        } else {
            init_seq_in_shell.erase(it->seq);
        }
        for (int seq : init_seq_in_shell) {
            outfile << "  " << seq;
        }
        outfile << std::endl;
    }
    outfile << "***************************" << endl;


}

void FirstCoordExchangeSearch::readInfo() {
    Atom::select2group(ids1, ids2);
    dist_cutoff = choose(0.0, GMX_DOUBLE_MAX, "Please enter distance cutoff:");
    tol_dist = choose(0.0, GMX_DOUBLE_MAX, "Please enter tol dist:");
    time_cutoff = choose(0.0, GMX_DOUBLE_MAX, "Please enter timecutoff:");
}

void FirstCoordExchangeSearch::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->ids1)) this->group1.insert(atom);
                      if (Atom::is_match(atom, this->ids2)) this->group2.insert(atom);
                  });
}


class RadicalDistribtuionFunction : public BasicAnalysis {

    double rmax;
    double width;

    bool intramol = false;
    int nframe = 0;
    int numj = 0, numk = 0;

    int nbin;
    map<int, int> hist;
    map<int, double> gr, gs, integral;

    double xbox, ybox, zbox;

    Atom::AtomIndenter ids1;
    Atom::AtomIndenter ids2;

    std::unordered_set<shared_ptr<Atom>> group1;
    std::unordered_set<shared_ptr<Atom>> group2;

public:
    RadicalDistribtuionFunction() {
        enable_outfile = true;
    }

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print() override;

    void readInfo() override;

    static const string title() {
        return "Radical Distribution Function";
    }
};


void RadicalDistribtuionFunction::process(std::shared_ptr<Frame> &frame) {


    nframe++;
    xbox = frame->a_axis;
    ybox = frame->b_axis;
    zbox = frame->c_axis;

    for (auto &atom_j : group1) {
        double xj = atom_j->x;
        double yj = atom_j->y;
        double zj = atom_j->z;
        auto molj = atom_j->molecule.lock();

        for (auto &atom_k : group2) {
            if (atom_j != atom_k) {
                auto mol_k = atom_k->molecule.lock();
                if (intramol or (molj not_eq mol_k)) {
                    double dx = atom_k->x - xj;
                    double dy = atom_k->y - yj;
                    double dz = atom_k->z - zj;
                    frame->image(dx, dy, dz);
                    double rjk = sqrt(dx * dx + dy * dy + dz * dz);
                    int ibin = int(rjk / width) + 1;
                    if (ibin <= nbin) {
                        hist[ibin] += 1;
                    }
                }
            }
        }
    }
}

void RadicalDistribtuionFunction::readInfo() {

    Atom::select2group(ids1, ids2);
    rmax = choose(0.0, GMX_DOUBLE_MAX, "Enter Maximum Distance to Accumulate[10.0 Ang]:", true, 10.0);
    width = choose(0.0, GMX_DOUBLE_MAX, "Enter Width of Distance Bins [0.01 Ang]:", true, 0.01);
    string inputline = input("Include Intramolecular Pairs in Distribution[N]:");
    boost::trim(inputline);
    if (!inputline.empty()) {
        if (boost::to_lower_copy(inputline) == "y") {
            intramol = true;
        }
    }
    nbin = int(rmax / width);
    for (int i = 0; i <= nbin; i++) {
        hist[i] = 0;
        gr[i] = gs[i] = 0.0;
    }

}

void RadicalDistribtuionFunction::print() {
    if (numj != 0 and numk != 0) {
        double factor = (4.0 / 3.0) * M_PI * nframe;
        int pairs = numj * numk;
        double volume = xbox * ybox * zbox;
        factor *= pairs / volume;
        for (int i = 1; i <= nbin; i++) {
            double rupper = i * width;
            double rlower = rupper - width;
            double expect = factor * (pow(rupper, 3) - pow(rlower, 3));
            gr[i] = hist[i] / expect;
            if (i == 1) {
                integral[i] = hist[i] / double(nframe);
            } else {
                integral[i] = hist[i] / double(nframe) + integral[i - 1];
            }
        }

    }

//     find the 5th degree polynomial smoothed distribution function

    if (nbin >= 5) {
        gs[1] = (69.0 * gr[1] + 4.0 * gr[2] - 6.0 * gr[3] + 4.0 * gr[4] - gr[5]) / 70.0;
        gs[2] = (2.0 * gr[1] + 27.0 * gr[2] + 12.0 * gr[3] - 8.0 * gr[4] + 2.0 * gr[5]) / 35.0;
        for (int i = 3; i <= nbin - 2; i++) {
            gs[i] = (-3.0 * gr[i - 2] + 12.0 * gr[i - 1] +
                     17.0 * gr[i] + 12.0 * gr[i + 1] - 3.0 * gr[i + 2]) / 35.0;
        }
        gs[nbin - 1] =
                (2.0 * gr[nbin - 4] - 8.0 * gr[nbin - 3] +
                 12.0 * gr[nbin - 2] + 27.0 * gr[nbin - 1] + 2.0 * gr[nbin]) / 35.0;
        gs[nbin] =
                (-gr[nbin - 4] + 4.0 * gr[nbin - 3] - 6.0 * gr[nbin - 2]
                 + 4.0 * gr[nbin - 1] + 69.0 * gr[nbin]) / 70.0;
        for (int i = 1; i <= nbin; i++) {
            gs[i] = max(0.0, gs[i]);
        }
    }

    outfile << "************************************************" << endl;
    outfile << "***** Pairwise Radial Distribution Function ****" << endl;

    outfile << "First Type : " << ids1 << " Second Type : " << ids2 << endl;

    outfile << "************************************************" << endl;
    outfile << "Bin    Counts    Distance    Raw g(r)  Smooth g(r)   Integral" << endl;

    for (int i = 1; i <= nbin; i++) {
        outfile << fmt::sprintf("%d        %d      %.4f      %.4f     %.4f      %.4f\n",
                                i, hist[i], (i - 0.5) * width, gr[i], gs[i], integral[i]);
    }

    outfile << "************************************************" << endl;
}

void RadicalDistribtuionFunction::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->ids1)) this->group1.insert(atom);
                      if (Atom::is_match(atom, this->ids2)) this->group2.insert(atom);
                  });
    numj = group1.size();
    numk = group2.size();
}

class RMSDCal : public BasicAnalysis {

    int steps = 0; // current frame number
    std::map<int, double> rmsd_map;
    bool first_frame = true;

    double x1[ATOM_MAX], y1[ATOM_MAX], z1[ATOM_MAX];
    double x2[ATOM_MAX], y2[ATOM_MAX], z2[ATOM_MAX];

    double rmsvalue(shared_ptr<Frame> &frame);

    Atom::AtomIndenter ids;

    std::unordered_set<shared_ptr<Atom>> group;


    void find_matched_atoms(shared_ptr<Frame> &frame) {
        if (first_frame) {
            std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                          [this](shared_ptr<Atom> &atom) {
                              if (Atom::is_match(atom, ids)) group.insert(atom);
                          });
        }
    }

public:
    RMSDCal() {
        enable_outfile = true;
    }

    void process(std::shared_ptr<Frame> &frame) override;

    void print() override;

    void readInfo() override;

    static const string title() {
        return "RMSD Calculator";
    }

    static double rmsfit(double x1[], double y1[], double z1[],
                         double x2[], double y2[], double z2[], int nfit);

    static void jacobi(int n, double a[4][4], double d[], double v[4][4]);

    static void quatfit(int n1, double x1[], double y1[], double z1[],
                        int n2, double x2[], double y2[], double z2[], int nfit);

    static void center(int n1, double x1[], double y1[], double z1[],
                       int n2, double x2[], double y2[], double z2[],
                       double mid[], int nfit);

};

void RMSDCal::process(std::shared_ptr<Frame> &frame) {
    find_matched_atoms(frame);
    steps++;
    rmsd_map[steps] = rmsvalue(frame);
}

void RMSDCal::print() {
    outfile << "***************************" << endl;
    outfile << "****** RMSD Calculator ****" << endl;
    outfile << "GROUP:" << ids << endl;
    outfile << "***************************" << endl;
    for (int cyc = 1; cyc <= steps; cyc++) {
        outfile << cyc << "     " << rmsd_map[cyc] << endl;
    }
    outfile << "***************************" << endl;
}

void RMSDCal::readInfo() {
    Atom::select1group(ids, "Please enter atom group:");
}


double RMSDCal::rmsvalue(shared_ptr<Frame> &frame) {

    int nfit, n;
    nfit = n = static_cast<int>(this->group.size());
    BOOST_ASSERT_MSG(n < ATOM_MAX, "need to increase ATOM_MAX");

    if (first_frame) {
        first_frame = false;
        int index = 0;
        bool first_atom = true;
        double first_x, first_y, first_z;
        for (auto &atom : this->group) {
            if (first_atom) {
                first_atom = false;
                first_x = x1[index] = atom->x;
                first_y = y1[index] = atom->y;
                first_z = z1[index] = atom->z;
            } else {
                double xr = atom->x - first_x;
                double yr = atom->y - first_y;
                double zr = atom->z - first_z;
                frame->image(xr, yr, zr);
                x1[index] = first_x + xr;
                y1[index] = first_y + yr;
                z1[index] = first_z + zr;
            }
            index++;
        }
        return 0.0;
    } else {
        int index = 0;
        bool first_atom = true;
        double first_x, first_y, first_z;
        for (auto &atom : this->group) {
            if (first_atom) {
                first_atom = false;
                first_x = x2[index] = atom->x;
                first_y = y2[index] = atom->y;
                first_z = z2[index] = atom->z;
            } else {
                double xr = atom->x - first_x;
                double yr = atom->y - first_y;
                double zr = atom->z - first_z;
                frame->image(xr, yr, zr);
                x2[index] = first_x + xr;
                y2[index] = first_y + yr;
                z2[index] = first_z + zr;
            }
            index++;
        }
    }

    double mid[3];
    center(n, x1, y1, z1, n, x2, y2, z2, mid, nfit);
    quatfit(n, x1, y1, z1, n, x2, y2, z2, nfit);
    double rms = rmsfit(x1, y1, z1, x2, y2, z2, nfit);
    return rms;
}

void RMSDCal::center(int n1, double x1[], double y1[], double z1[],
                     int n2, double x2[], double y2[], double z2[],
                     double mid[], int nfit) {
    mid[0] = mid[1] = mid[2] = 0.0;
    double norm = 0.0;
    for (int i = 0; i < nfit; ++i) {
        mid[0] += x2[i];
        mid[1] += y2[i];
        mid[2] += z2[i];
        norm += 1.0;
    }
    mid[0] /= norm;
    mid[1] /= norm;
    mid[2] /= norm;
    for (int i = 0; i < n2; ++i) {
        x2[i] -= mid[0];
        y2[i] -= mid[1];
        z2[i] -= mid[2];
    }

    mid[0] = mid[1] = mid[2] = 0.0;
    norm = 0.0;
    for (int i = 0; i < nfit; ++i) {
        mid[0] += x1[i];
        mid[1] += y1[i];
        mid[2] += z1[i];
        norm += 1.0;
    }
    mid[0] /= norm;
    mid[1] /= norm;
    mid[2] /= norm;
    for (int i = 0; i < n1; ++i) {
        x1[i] -= mid[0];
        y1[i] -= mid[1];
        z1[i] -= mid[2];
    }


}

void RMSDCal::quatfit(int /* n1 */, double x1[], double y1[], double z1[],
                      int n2, double x2[], double y2[], double z2[], int nfit) {
    int i;
    //    int i1, i2;
    //    double weigh;
    double xrot, yrot, zrot;
    double xxyx, xxyy, xxyz;
    double xyyx, xyyy, xyyz;
    double xzyx, xzyy, xzyz;
    double q[4], d[4];
    double rot[3][3];
    double c[4][4], v[4][4];

    xxyx = 0.0;
    xxyy = 0.0;
    xxyz = 0.0;
    xyyx = 0.0;
    xyyy = 0.0;
    xyyz = 0.0;
    xzyx = 0.0;
    xzyy = 0.0;
    xzyz = 0.0;

    for (i = 0; i < nfit; ++i) {
        xxyx += x1[i] * x2[i];
        xxyy += y1[i] * x2[i];
        xxyz += z1[i] * x2[i];
        xyyx += x1[i] * y2[i];
        xyyy += y1[i] * y2[i];
        xyyz += z1[i] * y2[i];
        xzyx += x1[i] * z2[i];
        xzyy += y1[i] * z2[i];
        xzyz += z1[i] * z2[i];
    }

    c[0][0] = xxyx + xyyy + xzyz;
    c[0][1] = xzyy - xyyz;
    c[1][1] = xxyx - xyyy - xzyz;
    c[0][2] = xxyz - xzyx;
    c[1][2] = xxyy + xyyx;
    c[2][2] = xyyy - xzyz - xxyx;
    c[0][3] = xyyx - xxyy;
    c[1][3] = xzyx + xxyz;
    c[2][3] = xyyz + xzyy;
    c[3][3] = xzyz - xxyx - xyyy;

    jacobi(4, c, d, v);

    q[0] = v[0][3];
    q[1] = v[1][3];
    q[2] = v[2][3];
    q[3] = v[3][3];

    rot[0][0] = q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3];
    rot[1][0] = 2.0 * (q[1] * q[2] - q[0] * q[3]);
    rot[2][0] = 2.0 * (q[1] * q[3] + q[0] * q[2]);
    rot[0][1] = 2.0 * (q[2] * q[1] + q[0] * q[3]);
    rot[1][1] = q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3];
    rot[2][1] = 2.0 * (q[2] * q[3] - q[0] * q[1]);
    rot[0][2] = 2.0 * (q[3] * q[1] - q[0] * q[2]);
    rot[1][2] = 2.0 * (q[3] * q[2] + q[0] * q[1]);
    rot[2][2] = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3];

    for (i = 0; i < n2; ++i) {
        xrot = x2[i] * rot[0][0] + y2[i] * rot[0][1] + z2[i] * rot[0][2];
        yrot = x2[i] * rot[1][0] + y2[i] * rot[1][1] + z2[i] * rot[1][2];
        zrot = x2[i] * rot[2][0] + y2[i] * rot[2][1] + z2[i] * rot[2][2];
        x2[i] = xrot;
        y2[i] = yrot;
        z2[i] = zrot;
    }

}

double RMSDCal::rmsfit(double x1[], double y1[], double z1[],
                       double x2[], double y2[], double z2[], int nfit) {

    double fit;
    double xr, yr, zr, dist2;
    double norm;

    fit = 0.0;
    norm = 0.0;
    for (int i = 0; i < nfit; ++i) {
        xr = x1[i] - x2[i];
        yr = y1[i] - y2[i];
        zr = z1[i] - z2[i];
        dist2 = xr * xr + yr * yr + zr * zr;
        norm += 1.0;
        fit += dist2;
    }
    return sqrt(fit / norm);
}

void RMSDCal::jacobi(int n, double a[4][4], double d[], double v[4][4]) {
    // taken from tinker

    int i, j, k;
    int ip, iq;
    int nrot, maxrot;
    double sm, tresh, s, c, t;
    double theta, tau, h, g, p;
    //    double *b; //traditional point
    //    double *z;
    double b[4];
    double z[4];
    //    b = new double[n];
    //    z = new double[n];

    maxrot = 100;
    nrot = 0;
    for (ip = 0; ip < n; ++ip) {
        for (iq = 0; iq < n; ++iq) {
            v[ip][iq] = 0.0;
        }
        v[ip][ip] = 1.0;
    }
    for (ip = 0; ip < n; ++ip) {
        b[ip] = a[ip][ip];
        d[ip] = b[ip];
        z[ip] = 0.0;
    }

    //  perform the jacobi rotations

    for (i = 0; i < maxrot; ++i) {
        sm = 0.0;
        for (ip = 0; ip < n - 1; ++ip) {
            for (iq = ip + 1; iq < n; ++iq) {
                sm += abs(a[ip][iq]);
            }
        }
        if (sm == 0.0) goto label_10;
        if (i < 3) {
            tresh = 0.2 * sm / (n * n);
        } else {
            tresh = 0.0;
        }
        for (ip = 0; ip < n - 1; ++ip) {
            for (iq = ip + 1; iq < n; ++iq) {
                g = 100.0 * abs(a[ip][iq]);
                if (i > 3 && abs(d[ip]) + g == abs(d[ip]) && abs(d[iq]) + g == abs(d[iq]))
                    a[ip][iq] = 0.0;
                else if (abs(a[ip][iq]) > tresh) {
                    h = d[iq] - d[ip];
                    if (abs(h) + g == abs(h))
                        t = a[ip][iq] / h;
                    else {
                        theta = 0.5 * h / a[ip][iq];
                        t = 1.0 / (abs(theta) + sqrt(1.0 + theta * theta));
                        if (theta < 0.0) t = -t;
                    }
                    c = 1.0 / sqrt(1.0 + t * t);
                    s = t * c;
                    tau = s / (1.0 + c);
                    h = t * a[ip][iq];
                    z[ip] -= h;
                    z[iq] += h;
                    d[ip] -= h;
                    d[iq] += h;
                    a[ip][iq] = 0.0;
                    for (j = 0; j <= ip - 1; ++j) {
                        g = a[j][ip];
                        h = a[j][iq];
                        a[j][ip] = g - s * (h + g * tau);
                        a[j][iq] = h + s * (g - h * tau);
                    }
                    for (j = ip + 1; j <= iq - 1; ++j) {
                        g = a[ip][j];
                        h = a[j][iq];
                        a[ip][j] = g - s * (h + g * tau);
                        a[j][iq] = h + s * (g - h * tau);
                    }
                    for (j = iq + 1; j < n; ++j) {
                        g = a[ip][j];
                        h = a[iq][j];
                        a[ip][j] = g - s * (h + g * tau);
                        a[iq][j] = h + s * (g - h * tau);
                    }
                    for (j = 0; j < n; ++j) {
                        g = v[j][ip];
                        h = v[j][iq];
                        v[j][ip] = g - s * (h + g * tau);
                        v[j][iq] = h + s * (g - h * tau);
                    }
                    ++nrot;
                }
            }
        }
        for (ip = 0; ip < n; ++ip) {
            b[ip] += z[ip];
            d[ip] = b[ip];
            z[ip] = 0.0;
        }
    }

    label_10:
    //    delete [] b; b = nullptr;
    //    delete [] z; z = nullptr;

    if (nrot == maxrot)
        std::cerr << " JACOBI  --  Matrix Diagonalization not Converged" << std::endl;

    for (i = 0; i < n - 1; ++i) {
        k = i;
        p = d[i];
        for (j = i + 1; j < n; ++j) {
            if (d[j] < p) {
                k = j;
                p = d[j];
            }
        }
        if (k != i) {
            d[k] = d[i];
            d[i] = p;
            for (j = 0; j < n; ++j) {
                p = v[j][i];
                v[j][i] = v[j][k];
                v[j][k] = p;
            }
        }
    }


}

class RMSFCal : public BasicAnalysis {


    Atom::AtomIndenter ids;

    std::unordered_set<shared_ptr<Atom>> group;


    void find_matched_atoms(shared_ptr<Frame> &frame) {
        if (first_frame) {
            std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                          [this](shared_ptr<Atom> &atom) {
                              if (Atom::is_match(atom, ids)) group.insert(atom);
                          });
        }
    }

    int steps = 0; // current frame number
    bool first_frame = true;

    static double rmsfit(double x1[], double y1[], double z1[],
                         double x2[], double y2[], double z2[], int nfit) {
        return RMSDCal::rmsfit(x1, y1, z1, x2, y2, z2, nfit);
    }

    static void jacobi(int n, double a[4][4], double d[], double v[4][4]) {
        RMSDCal::jacobi(n, a, d, v);
    }

    static void quatfit(int n1, double x1[], double y1[], double z1[],
                        int n2, double x2[], double y2[], double z2[], int nfit) {
        RMSDCal::quatfit(n1, x1, y1, z1, n2, x2, y2, z2, nfit);
    }

    static void center(int n, double x[], double y[], double z[], double mid[], int nfit) {
        mid[0] = mid[1] = mid[2] = 0.0;
        double norm = 0.0;
        for (int i = 0; i < nfit; ++i) {
            mid[0] += x[i];
            mid[1] += y[i];
            mid[2] += z[i];
            norm += 1.0;
        }
        mid[0] /= norm;
        mid[1] /= norm;
        mid[2] /= norm;

        for (int i = 0; i < n; ++i) {
            x[i] -= mid[0];
            y[i] -= mid[1];
            z[i] -= mid[2];
        }
    }

    double x_avg[ATOM_MAX], y_avg[ATOM_MAX], z_avg[ATOM_MAX];
    double x1[ATOM_MAX], y1[ATOM_MAX], z1[ATOM_MAX];
    double x2[ATOM_MAX], y2[ATOM_MAX], z2[ATOM_MAX];

    double rmsvalue(int index);

    map<int, map<int, double>> x, y, z;

public:
    RMSFCal() {
        enable_outfile = true;
    }

    void process(std::shared_ptr<Frame> &frame) override;

    void print() override;

    void readInfo() override;

    static const string title() {
        return "RMSF Calculator";
    }
};

void RMSFCal::process(std::shared_ptr<Frame> &frame) {
    steps++;
    int nfit, n;
    nfit = n = static_cast<int>(this->group.size());
    BOOST_ASSERT_MSG(group.size() < ATOM_MAX, "need to increase ATOM_MAX");

    if (!first_frame) {
        first_frame = false;
        int index = 0;
        bool first_atom = true;
        double first_x, first_y, first_z;
        for (auto &atom : this->group) {
            if (first_atom) {
                first_atom = false;
                first_x = x1[index] = atom->x;
                first_y = y1[index] = atom->y;
                first_z = z1[index] = atom->z;
            } else {
                double xr = atom->x - first_x;
                double yr = atom->y - first_y;
                double zr = atom->z - first_z;
                frame->image(xr, yr, zr);
                x1[index] = first_x + xr;
                y1[index] = first_y + yr;
                z1[index] = first_z + zr;
            }
            index++;
        }

        double mid[3];
        center(n, x1, y1, z1, mid, nfit);
        map<int, double> f1x, f1y, f1z;
        for (unsigned index = 0; index < this->group.size(); index++) {
            x_avg[index] = x1[index];
            y_avg[index] = y1[index];
            z_avg[index] = z1[index];
            f1x[index] = x1[index];
            f1y[index] = y1[index];
            f1z[index] = z1[index];
        }
        x[steps] = f1x;
        y[steps] = f1y;
        z[steps] = f1z;
    } else {
        int index = 0;
        bool first_atom = true;
        double first_x, first_y, first_z;
        for (auto &atom : this->group) {
            if (first_atom) {
                first_atom = false;
                first_x = x2[index] = atom->x;
                first_y = y2[index] = atom->y;
                first_z = z2[index] = atom->z;
            } else {
                double xr = atom->x - first_x;
                double yr = atom->y - first_y;
                double zr = atom->z - first_z;
                frame->image(xr, yr, zr);
                x2[index] = first_x + xr;
                y2[index] = first_y + yr;
                z2[index] = first_z + zr;
            }
            index++;
        }
        double mid[3];
        center(n, x2, y2, z2, mid, nfit);
        quatfit(n, x1, y1, z1, n, x2, y2, z2, nfit);
        map<int, double> f2x, f2y, f2z;
        for (unsigned index = 0; index < this->group.size(); index++) {
            x_avg[index] += x2[index];
            y_avg[index] += y2[index];
            z_avg[index] += z2[index];
            f2x[index] = x2[index];
            f2y[index] = y2[index];
            f2z[index] = z2[index];
        }
        x[steps] = f2x;
        y[steps] = f2y;
        z[steps] = f2z;
    }

}

void RMSFCal::print() {

    for (unsigned int index = 0; index < this->group.size(); index++) {
        x_avg[index] /= this->steps;
        y_avg[index] /= this->steps;
        z_avg[index] /= this->steps;
    }

    outfile << "***************************" << endl;
    outfile << "****** RMSF Calculator ****" << endl;
    outfile << "SET:" << ids << endl;
    outfile << "***************************" << endl;
    int index = 0;
    for (auto &at : group) {
        outfile << at->seq << "     " << rmsvalue(index) << endl;
        index++;
    }
    outfile << "***************************" << endl;
}

void RMSFCal::readInfo() {
    Atom::select1group(ids, "Please enter group:");
}

double RMSFCal::rmsvalue(int index) {

    double dx2_y2_z2 = 0.0;
    double dx, dy, dz;
    for (int frame = 1; frame <= steps; frame++) {
        dx = x[frame][index] - x_avg[index];
        dy = y[frame][index] - y_avg[index];
        dz = z[frame][index] - z_avg[index];
        dx2_y2_z2 += dx * dx + dy * dy + dz * dz;
    }

    return sqrt(dx2_y2_z2 / steps);

}

enum class Symbol {
    Hydrogen,
    Carbon,
    Nitrogen,
    Oxygen,
    X,
    Unknown
};


Symbol which(const shared_ptr<Atom> &atom) {
    double mass = forcefield.find_mass(atom);
    if (mass < 2.0)
        return Symbol::Hydrogen;
    else if (mass > 11.5 and mass < 13.0)
        return Symbol::Carbon;
    else if (mass > 13.0 and mass < 15.0)
        return Symbol::Nitrogen;
    else if (mass >= 15.0 and mass < 17.0)
        return Symbol::Oxygen;
    else return Symbol::Unknown;
}

enum class Selector {
    Acceptor,
    Donor,
    Both
};

enum class HBondType {
    VMDVerion,
    GMXVersion
};

class HBond : public BasicAnalysis {

public:
    HBond() {
        enable_outfile = true;
        enable_forcefield = true;
    }

    void process(std::shared_ptr<Frame> &frame) override;

    void print() override;

    void readInfo() override;

    static const string title() {
        return "Hydrogen Bond";
    }

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

private:
    std::map<int, int> hbonds;

    HBondType hbond_type = HBondType::VMDVerion;


    void Selector_Both(std::shared_ptr<Frame> &frame);

    void Selector_Donor(std::shared_ptr<Frame> &frame);

    void Selector_Acceptor(std::shared_ptr<Frame> &frame);

    Atom::AtomIndenter ids1;
    Atom::AtomIndenter ids2;

    std::unordered_set<shared_ptr<Atom>> group1;
    std::unordered_set<shared_ptr<Atom>> group2;

    Selector mode;
    double donor_acceptor_dist_cutoff;
    double angle_cutoff;
    int steps = 0;

};


void HBond::print() {
    outfile << "***************************" << endl;
    outfile << "****** Hydrogen Bond ******" << endl;
    outfile << "TYPE:";
    switch (mode) {
        case Selector::Both:
            outfile << "Both" << endl;
            outfile << "SET1:" << ids1 << endl;
            break;
        case Selector::Donor:
            outfile << "Donor" << endl;
            outfile << "SET1:" << ids1 << endl;
            outfile << "SET2:" << ids2 << endl;
            break;
        case Selector::Acceptor:
            outfile << "Acceptor" << endl;
            outfile << "SET1:" << ids1 << endl;
            outfile << "SET2:" << ids2 << endl;
            break;
    }
    switch (hbond_type) {
        case HBondType::VMDVerion:
            outfile << "HBond criteria : VMD version" << endl;
            break;
        case HBondType::GMXVersion:
            outfile << "HBond criteria : GMX version" << endl;
            break;
    }
    outfile << "distance cutoff : " << this->donor_acceptor_dist_cutoff << endl;
    outfile << "angle cutoff : " << this->angle_cutoff << endl;
    outfile << "***************************" << endl;
    for (int cyc = 1; cyc <= steps; cyc++) {
        outfile << cyc << "   " << hbonds[cyc] << endl;
    }
    outfile << "***************************" << endl;
}

void HBond::process(std::shared_ptr<Frame> &frame) {


    steps++;
    hbonds[steps] = 0;
    switch (mode) {
        case Selector::Both:
            Selector_Both(frame);
            break;
        case Selector::Donor:
            Selector_Donor(frame);
            break;
        case Selector::Acceptor:
            Selector_Acceptor(frame);
            break;
    }

}


void HBond::readInfo() {

    Atom::select1group(ids1, "Please enter group:");

    auto input_line = input("Which selector: [(1)Acceptor | (2)Donor | (3)Both]:");
    auto mode_no = stoi(input_line);
    switch (mode_no) {
        case 1:
            this->mode = Selector::Acceptor;
            break;
        case 2:
            this->mode = Selector::Donor;
            break;
        case 3:
            this->mode = Selector::Both;
            break;
        default:
            cerr << "Error type" << endl;
            exit(1);
    }
    if (this->mode not_eq Selector::Both) {
        Atom::select1group(ids2, "Please enter group2:");
    }
    cout << "HBond criteria:\n(1)VMD version\n(2)GMX version\n";
    switch (choose(1, 2, "choose:")) {
        case 1:
            hbond_type = HBondType::VMDVerion;
            break;
        case 2:
            hbond_type = HBondType::GMXVersion;
            break;
        default:
            cerr << "wrong type! \n";
            exit(2);
    }

    this->donor_acceptor_dist_cutoff =
            choose(0.0, static_cast<double>(std::numeric_limits<int>::max()), "Donor-Acceptor Distance:");
    this->angle_cutoff =
            choose(0.0, static_cast<double>(std::numeric_limits<int>::max()), "Angle cutoff:");

}

void HBond::Selector_Acceptor(std::shared_ptr<Frame> &frame) {
    for (auto &atom2 : group2) {
        // Atom2
        if (which(atom2) == Symbol::Hydrogen) {
            // Atom1
            auto atom1 = frame->atom_map[atom2->con_list.front()];
            auto atom1_symbol = which(atom1);
            if (atom1_symbol == Symbol::Nitrogen or atom1_symbol == Symbol::Oxygen) {
                for (auto &atom3 : group1) {
                    // Atom3
                    auto atom3_symbol = which(atom3);
                    if (atom3_symbol == Symbol::Nitrogen or atom3_symbol == Symbol::Oxygen) {
                        auto distance = atom_distance(atom1, atom3, frame);
                        if (distance <= this->donor_acceptor_dist_cutoff) {
                            auto vx1 = atom2->x - atom1->x;
                            auto vy1 = atom2->y - atom1->y;
                            auto vz1 = atom2->z - atom1->z;
                            frame->image(vx1, vy1, vz1);
                            auto len1 = std::sqrt(vx1 * vx1 + vy1 * vy1 + vz1 * vz1);

                            double vx2, vy2, vz2;
                            switch (hbond_type) {
                                case HBondType::VMDVerion:
                                    // atom2 -> atom3
                                    vx2 = atom3->x - atom2->x;
                                    vy2 = atom3->y - atom2->y;
                                    vz2 = atom3->z - atom2->z;
                                    break;
                                case HBondType::GMXVersion:
                                    // atom1 -> atom3
                                    vx2 = atom3->x - atom1->x;
                                    vy2 = atom3->y - atom1->y;
                                    vz2 = atom3->z - atom1->z;
                                    break;
                            }

                            frame->image(vx2, vy2, vz2);
                            auto len2 = std::sqrt(vx2 * vx2 + vy2 * vy2 + vz2 * vz2);
                            auto cosine = (vx1 * vx2 + vy1 * vy2 + vz1 * vz2) / (len1 * len2);
                            auto ang = radian * std::acos(cosine);
                            if (std::abs(ang) <= this->angle_cutoff) {
                                hbonds[steps]++;
                            }
                        }
                    }
                }
            }
        }
    }
}

void HBond::Selector_Donor(std::shared_ptr<Frame> &frame) {
    for (auto &atom2 : group1) {
        // Atom2
        if (which(atom2) == Symbol::Hydrogen) {
            // Atom1
            auto atom1 = frame->atom_map[atom2->con_list.front()];
            auto atom1_symbol = which(atom1);
            if (atom1_symbol == Symbol::Nitrogen or atom1_symbol == Symbol::Oxygen) {
                for (auto &atom3: group2) {
                    // Atom3
                    auto atom3_symbol = which(atom3);
                    if (atom3_symbol == Symbol::Nitrogen or atom3_symbol == Symbol::Oxygen) {
                        auto distance = atom_distance(atom1, atom3, frame);
                        if (distance <= this->donor_acceptor_dist_cutoff) {
                            auto vx1 = atom2->x - atom1->x;
                            auto vy1 = atom2->y - atom1->y;
                            auto vz1 = atom2->z - atom1->z;
                            frame->image(vx1, vy1, vz1);
                            auto len1 = std::sqrt(vx1 * vx1 + vy1 * vy1 + vz1 * vz1);

                            double vx2, vy2, vz2;
                            switch (hbond_type) {
                                case HBondType::VMDVerion:
                                    // atom2 -> atom3
                                    vx2 = atom3->x - atom2->x;
                                    vy2 = atom3->y - atom2->y;
                                    vz2 = atom3->z - atom2->z;
                                    break;
                                case HBondType::GMXVersion:
                                    // atom1 -> atom3
                                    vx2 = atom3->x - atom1->x;
                                    vy2 = atom3->y - atom1->y;
                                    vz2 = atom3->z - atom1->z;
                                    break;
                            }

                            frame->image(vx2, vy2, vz2);
                            auto len2 = std::sqrt(vx2 * vx2 + vy2 * vy2 + vz2 * vz2);
                            auto cosine = (vx1 * vx2 + vy1 * vy2 + vz1 * vz2) / (len1 * len2);
                            auto ang = radian * std::acos(cosine);
                            if (std::abs(ang) <= this->angle_cutoff) {
                                hbonds[steps]++;
                            }
                        }
                    }
                }
            }
        }
    }
}

void HBond::Selector_Both(std::shared_ptr<Frame> &frame) {
    for (auto &atom2 : group1) {
        // Atom2
        if (which(atom2) == Symbol::Hydrogen) {
            // Atom1
            auto &atom1 = frame->atom_map[atom2->con_list.front()];
            auto atom1_symbol = which(atom1);
            if (atom1_symbol == Symbol::Nitrogen or atom1_symbol == Symbol::Oxygen) {
                for (auto &atom3 : group1) {
                    // atom1 and atom3 can not the same atom
                    if (atom1 == atom3) continue;
                    auto atom3_symbol = which(atom3);
                    if (atom3_symbol == Symbol::Nitrogen or atom3_symbol == Symbol::Oxygen) {
                        auto distance = atom_distance(atom1, atom3, frame);
                        if (distance <= this->donor_acceptor_dist_cutoff) {
                            if (atom1 != atom3 && !atom1->adj(atom3)) {
                                auto vx1 = atom2->x - atom1->x;
                                auto vy1 = atom2->y - atom1->y;
                                auto vz1 = atom2->z - atom1->z;
                                // atom1 -> atom2
                                frame->image(vx1, vy1, vz1);
                                auto len1 = std::sqrt(vx1 * vx1 + vy1 * vy1 + vz1 * vz1);

                                double vx2, vy2, vz2;
                                switch (hbond_type) {
                                    case HBondType::VMDVerion:
                                        // atom2 -> atom3
                                        vx2 = atom3->x - atom2->x;
                                        vy2 = atom3->y - atom2->y;
                                        vz2 = atom3->z - atom2->z;
                                        break;
                                    case HBondType::GMXVersion:
                                        // atom1 -> atom3
                                        vx2 = atom3->x - atom1->x;
                                        vy2 = atom3->y - atom1->y;
                                        vz2 = atom3->z - atom1->z;
                                        break;
                                }

                                frame->image(vx2, vy2, vz2);
                                auto len2 = std::sqrt(vx2 * vx2 + vy2 * vy2 + vz2 * vz2);

                                auto cosine = (vx1 * vx2 + vy1 * vy2 + vz1 * vz2) / (len1 * len2);
                                auto ang = radian * std::acos(cosine);
                                if (std::abs(ang) <= this->angle_cutoff) {
                                    hbonds[steps]++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void HBond::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->ids1)) this->group1.insert(atom);
                      if (this->mode not_eq Selector::Both && Atom::is_match(atom, this->ids2))
                          this->group2.insert(atom);
                  });
}

class ResidenceTime : public BasicAnalysis {

public:

    ResidenceTime() {
        enable_tbb = true;
        enable_outfile = true;
    }

    ~ResidenceTime() {
        boost::checked_array_delete(time_array);
        boost::checked_array_delete(Rt_array);
        if (mark) {
            for (int i = 0; i < atom_num; ++i)
                boost::checked_array_delete(mark[i]);
            boost::checked_array_delete(mark);
        }
    };

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print() override;

    void readInfo() override;

    static const string title() {
        return "ResidenceTime";
    }

private:

    std::map<int, std::list<bool>> mark_map;

    void calculate();

    void setSteps(size_t steps, int atom_num);

    double dis_cutoff;
    size_t steps = 0;
    int atom_num = 0;
    int *time_array = nullptr;
    double *Rt_array = nullptr;
    int **mark = nullptr;
    double time_star = 0;

    Atom::AtomIndenter ids1;
    Atom::AtomIndenter ids2;

    std::unordered_set<shared_ptr<Atom>> group1;
    std::unordered_set<shared_ptr<Atom>> group2;
};

void ResidenceTime::setSteps(size_t steps, int atom_num) {
    this->steps = steps;
    delete[] time_array;
    time_array = new int[steps - 1];
    bzero(time_array, sizeof(int) * (steps - 1));

    delete[] Rt_array;
    Rt_array = new double[steps - 1];
    bzero(Rt_array, sizeof(double) * (steps - 1));
    this->atom_num = atom_num;
    mark = new int *[atom_num];
    for (int i = 0; i < atom_num; ++i)
        mark[i] = new int[steps];

}

void ResidenceTime::calculate() {

    auto atom_star_map = new std::list<pair<unsigned int, unsigned int>>[atom_num];
    for (int atom = 0; atom < atom_num; atom++) {
        unsigned int count1 = 0;
        bool swi = true;
        for (unsigned int k = 0; k < steps; k++) {
            if ((!mark[atom][k]) && swi) {
                swi = false;
                count1 = k;
            } else if (mark[atom][k] && (!swi)) {
                swi = true;
                atom_star_map[atom].emplace_back(count1, k - 1);
            }
        }

    }
    class Body {
    public:
        size_t steps;
        int atom_num;
        int *time_array = nullptr;
        double *Rt_array = nullptr;
        int **mark = nullptr;
        double time_star;
        std::list<std::pair<unsigned int, unsigned int>> *atom_star_map;

        Body(size_t &steps, int &atom_num, int * /* time_array */,
             double * /* Rt_array */, int **mark, double time_star,
             std::list<std::pair<unsigned int, unsigned int>> *atom_star_map) :
                steps(steps), atom_num(atom_num), mark(mark),
                time_star(time_star), atom_star_map(atom_star_map) {
            this->time_array = new int[steps - 1];
            bzero(this->time_array, sizeof(int) * (steps - 1));
            this->Rt_array = new double[steps - 1];
            bzero(this->Rt_array, sizeof(double) * (steps - 1));
        }

        Body(Body &body, tbb::split) :
                steps(body.steps), atom_num(body.atom_num), mark(body.mark),
                time_star(body.time_star), atom_star_map(body.atom_star_map) {
            this->time_array = new int[body.steps - 1];
            bzero(this->time_array, sizeof(int) * (body.steps - 1));
            this->Rt_array = new double[body.steps - 1];
            bzero(this->Rt_array, sizeof(double) * (body.steps - 1));
        }

        ~Body() {
            delete[] time_array;
            delete[] Rt_array;
        }

        void join(const Body &y) {
            for (unsigned int step = 0; step < steps - 1; step++) {
                time_array[step] += y.time_array[step];
                Rt_array[step] += y.Rt_array[step];
            }
        }

        void operator()(const tbb::blocked_range<size_t> &r) const {
            for (size_t i = r.begin(); i != r.end(); ++i) {
                int CN = 0;
                for (int atom = 0; atom < atom_num; atom++)
                    CN += mark[atom][i];
                if (CN == 0) {
                    std::cerr << "Warning !!! Coordination Number is zero,  skip frame (start from 0) = " << i << "\n";
                    continue;
                }
                for (size_t j = i + 1; j < steps; j++) {
                    double value = 0.0;
                    for (int atom = 0; atom < atom_num; atom++) {
                        if (mark[atom][i] && mark[atom][j]) {
                            auto &li = atom_star_map[atom];
                            unsigned int maxcount = 0;
                            for (auto &pi : li) {
                                unsigned int a = pi.first;
                                unsigned int b = pi.second;
                                if (i < a and j < b) continue;
                                else if (i < a and b < j) {
                                    if (maxcount < b - a) maxcount = b - a;
                                } else if (a < i and b < j) continue;
                                else {
                                    cerr << fmt::sprintf(" i = %d, j = %d, a = %d, d = %d\n", i, j, a, b);
                                    cerr << "error " << __FILE__ << " : " << __LINE__ << endl;
                                    exit(1);
                                }
                            }
                            if (maxcount <= time_star) value++;
                        }
                    }
                    time_array[j - i - 1]++;
                    Rt_array[j - i - 1] += value / CN;
                }
            }

        }
    } body(steps, atom_num, time_array, Rt_array, mark, time_star, atom_star_map);
    tbb::parallel_reduce(tbb::blocked_range<size_t>(0, steps - 1), body, tbb::auto_partitioner());
    for (unsigned int step = 0; step < steps - 1; step++) {
        Rt_array[step] = body.Rt_array[step];
        time_array[step] = body.time_array[step];
    }

    for (unsigned int i = 0; i < steps - 1; i++)
        Rt_array[i] /= time_array[i];
    delete[] atom_star_map;
}


void ResidenceTime::process(std::shared_ptr<Frame> &frame) {

    steps++;
    int atom_no = 0;
    for (auto &atom1 : group1) {
        for (auto &atom2 : group2) {
            double xr = atom1->x - atom2->x;
            double yr = atom1->y - atom2->y;
            double zr = atom1->z - atom2->z;
            frame->image(xr, yr, zr);
            double dist = std::sqrt(xr * xr + yr * yr + zr * zr);
            mark_map[atom_no].push_back(dist <= dis_cutoff);
            atom_no += 1;
        }
    }
}

void ResidenceTime::print() {
    if (steps < 2) {
        cerr << "Too few frame number :" << steps << endl;
        exit(1);
    }
    setSteps(steps, static_cast<int>(mark_map.size()));
    for (const auto &it : mark_map) {
        auto atom_no = it.first;
        int cyc = 0;
        for (auto value : it.second) {
            mark[atom_no][cyc] = int(value);
            cyc++;
        }
    }
    calculate();

    outfile << "# typ1: " << ids1 << ",  typ2: " << ids2 << "  dist_cutoff = " << dis_cutoff << " t* = "
            << time_star
            << std::endl;

    outfile << "# Frame        R " << std::endl;
    for (unsigned int i = 0; i < steps - 1; i++)
        outfile << i + 1 << "    " << Rt_array[i] << std::endl;
}

void ResidenceTime::readInfo() {

    Atom::select2group(ids1, ids2);

    std::cout << "Please enter distance cutoff1:";
    std::cin >> dis_cutoff;
    std::cout << "Please enter t*: ( unit: frame)";
    std::cin >> time_star;

}

void ResidenceTime::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->ids1)) this->group1.insert(atom);
                      if (Atom::is_match(atom, this->ids2)) this->group2.insert(atom);
                  });
}


// Use Green-Kubo equation to calculate self-diffuse coefficients
class GreenKubo : public BasicAnalysis {

public:
    GreenKubo() {
        enable_read_velocity = true;
        enable_tbb = true;
        enable_outfile = true;
    }

    void process(shared_ptr<Frame> &) override;

    void print() override;

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void readInfo() override;

    static const string title() {
        return "Green-Kubo";
    }

private:
    Atom::AtomIndenter ids;
    std::unordered_set<shared_ptr<Atom>> group;

    double timestep;
    size_t steps = 0;

    map<int, double> vecx_map;
    map<int, double> vecy_map;
    map<int, double> vecz_map;
};


void GreenKubo::print() {

    if (steps < 2) {
        cerr << "Too few frame number :" << steps << endl;
        exit(1);
    }
    auto vecx = new double[steps];
    auto vecy = new double[steps];
    auto vecz = new double[steps];
    for (unsigned int i = 0; i < steps; i++) {
        vecx[i] = vecx_map[i + 1];
        vecy[i] = vecy_map[i + 1];
        vecz[i] = vecz_map[i + 1];
    }

    class Body {
    public:
        double *vecx;
        double *vecy;
        double *vecz;

        double *cxx = nullptr;
        double *cyy = nullptr;
        double *czz = nullptr;
        int *numbers = nullptr;
        size_t steps;

        Body(double *vecx, double *vecy, double *vecz, size_t steps) :
                vecx(vecx), vecy(vecy), vecz(vecz), steps(steps) {
            allocate();
        }

        Body(const Body &c, tbb::split) :
                vecx(c.vecx), vecy(c.vecy), vecz(c.vecz), steps(c.steps) {
            allocate();
        }

        void allocate() {
            cxx = new double[steps];
            cyy = new double[steps];
            czz = new double[steps];
            numbers = new int[steps];
            bzero(cxx, sizeof(double) * steps);
            bzero(cyy, sizeof(double) * steps);
            bzero(czz, sizeof(double) * steps);
            bzero(numbers, sizeof(int) * steps);
        }

        ~Body() {
            delete[] cxx;
            delete[] cyy;
            delete[] czz;
            delete[] numbers;
        }

        void join(const Body &y) {
            for (unsigned int step = 0; step < steps - 1; step++) {
                cxx[step] += y.cxx[step];
                cyy[step] += y.cyy[step];
                czz[step] += y.czz[step];
                numbers[step] += y.numbers[step];
            }
        }

        void operator()(const tbb::blocked_range<size_t> &r) const {
            for (size_t i = r.begin(); i != r.end(); i++) {
                for (size_t j = i; j < steps; j++) {
                    cxx[j - i] += vecx[i] * vecx[j];
                    cyy[j - i] += vecy[i] * vecy[j];
                    czz[j - i] += vecz[i] * vecz[j];
                    numbers[j - i]++;
                }
            }
        }


    } body(vecx, vecy, vecz, steps);

    tbb::parallel_reduce(tbb::blocked_range<size_t>(0, steps - 1), body, tbb::auto_partitioner());


    for (unsigned int i = 0; i < steps; i++) {
        int no = body.numbers[i];
        body.cxx[i] /= no;
        body.cyy[i] /= no;
        body.czz[i] /= no;
    }
    outfile << "Green-Kubo self-diffuse " << endl;
    outfile << "selected group : " << ids << endl;
    outfile << "Time (ps)             DA(10^-9 m2/s)" << endl;
    double pre_c = body.cxx[0] + body.cyy[0] + body.czz[0];
    double integral = 0.0;
    for (unsigned int i = 1; i < steps; i++) {
        double c = body.cxx[i] + body.cyy[i] + body.czz[i];
        integral += 5 * (pre_c + c) * timestep / 3;
        pre_c = c;
        outfile << i * timestep << "    " << integral << endl;
    }
    delete[] vecx;
    delete[] vecy;
    delete[] vecz;
}


void GreenKubo::readInfo() {
    Atom::select1group(ids, "Please enter atom group");
    this->timestep = choose(0.0, static_cast<double>(numeric_limits<int>::max()),
                            "Please enter time step for each frame(ps):");

}

void GreenKubo::process(shared_ptr<Frame> &) {


    steps++;
    double xv, yv, zv;
    xv = yv = zv = 0.0;
    for (auto &atom_ptr : group) {

        xv += atom_ptr->vx;
        yv += atom_ptr->vy;
        zv += atom_ptr->vz;
        break;

    }

    auto atom_nums = group.size();
    xv /= atom_nums;
    yv /= atom_nums;
    zv /= atom_nums;

    vecx_map[steps] = xv;
    vecy_map[steps] = yv;
    vecz_map[steps] = zv;

}

void GreenKubo::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->ids)) this->group.insert(atom);
                  });
}

class Cluster : public BasicAnalysis {

public:
    Cluster() {
        enable_tbb = true;
        enable_outfile = true;
    }

    void process(std::shared_ptr<Frame> &frame) override;

    void print() override;

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void readInfo() override;

    static const string title() {
        return "Cluster Analysis(linkage)";
    }

private:

    Atom::AtomIndenter ids;
    std::unordered_set<shared_ptr<Atom>> group;


    double cutoff = 0.0;
    int steps = 0; // current frame number
    map<pair<int, int>, double> rmsd_map;

    double rmsfit(double x1[], double y1[], double z1[],
                  double x2[], double y2[], double z2[], int nfit) {
        return RMSDCal::rmsfit(x1, y1, z1, x2, y2, z2, nfit);
    }

    void jacobi(int n, double a[4][4], double d[], double v[4][4]) {
        RMSDCal::jacobi(n, a, d, v);
    }

    void quatfit(int n1, double x1[], double y1[], double z1[],
                 int n2, double x2[], double y2[], double z2[], int nfit) {
        RMSDCal::quatfit(n1, x1, y1, z1, n2, x2, y2, z2, nfit);
    }

    void center(int n1, double x1[], double y1[], double z1[],
                int n2, double x2[], double y2[], double z2[],
                double mid[], int nfit) {
        RMSDCal::center(n1, x1, y1, z1, n2, x2, y2, z2, mid, nfit);
    }

    map<int, map<int, double>> x, y, z;

    double rmsvalue(int index1, int index2);

};


void Cluster::process(std::shared_ptr<Frame> &frame) {


    map<int, double> xx, yy, zz;
    int index = 0;
    bool first_atom = true;
    double first_x, first_y, first_z;
    BOOST_ASSERT_MSG(group.size() < ATOM_MAX, "need to increase ATOM_MAX");
    for (auto &atom : group) {

        if (first_atom) {
            first_atom = false;
            first_x = xx[index] = atom->x;
            first_y = yy[index] = atom->y;
            first_z = zz[index] = atom->z;
        } else {
            double xr = atom->x - first_x;
            double yr = atom->y - first_y;
            double zr = atom->z - first_z;
            frame->image(xr, yr, zr);
            xx[index] = first_x + xr;
            yy[index] = first_y + yr;
            zz[index] = first_z + zr;
        }
        index++;
    }
    x[steps] = xx;
    y[steps] = yy;
    z[steps] = zz;
    steps++;
}

double Cluster::rmsvalue(int index1, int index2) {
    int nfit, n;
    nfit = n = group.size();

    double x1[ATOM_MAX], y1[ATOM_MAX], z1[ATOM_MAX];
    double x2[ATOM_MAX], y2[ATOM_MAX], z2[ATOM_MAX];

    auto &f1x = x[index1];
    auto &f1y = y[index1];
    auto &f1z = z[index1];
    for (int i = 0; i < n; i++) {
        x1[i] = f1x[i];
        y1[i] = f1y[i];
        z1[i] = f1z[i];

    }

    auto &f2x = x[index2];
    auto &f2y = y[index2];
    auto &f2z = z[index2];
    for (int i = 0; i < n; i++) {
        x2[i] = f2x[i];
        y2[i] = f2y[i];
        z2[i] = f2z[i];

    }

    double mid[3];
    center(n, x1, y1, z1, n, x2, y2, z2, mid, nfit);
    quatfit(n, x1, y1, z1, n, x2, y2, z2, nfit);
    double rms = rmsfit(x1, y1, z1, x2, y2, z2, nfit);
    return rms;
}

class rms_m {
public:
    rms_m(int i, int j, double rms) : i(i), j(j), rms(rms) {}

    int i;
    int j;
    double rms;
};

class rms_m2 {
public:
    rms_m2(int conf, int clust) : conf(conf), clust(clust) {};
    int conf;
    int clust;
};

void Cluster::print() {

    //  list<rms_m> rms_list;


    class CalCore {
    public:
        list<rms_m> local_rms_list;
        int steps;
        double cutoff;

        Cluster *parent;

        CalCore(int steps, double cutoff, Cluster *parent) : steps(steps), cutoff(cutoff), parent(parent) {};

        CalCore(CalCore &c, tbb::split) {
            steps = c.steps;
            cutoff = c.cutoff;
            parent = c.parent;
        }

        void join(CalCore &c) {
            local_rms_list.merge(c.local_rms_list,
                                 [](const rms_m &m1, const rms_m &m2) { return (m1.rms < m2.rms); });
        }

        void operator()(const tbb::blocked_range<int> &range) {
            double rms;
            for (int index1 = range.begin(); index1 != range.end(); ++index1) {
                for (int index2 = index1 + 1; index2 < steps; ++index2) {
                    rms = parent->rmsvalue(index1, index2);
                    if (rms < cutoff) local_rms_list.emplace_back(index1, index2, rms);
                }
            }
            local_rms_list.sort([](const rms_m &m1, const rms_m &m2) { return (m1.rms < m2.rms); });
        }
    } core(steps, cutoff, this);

    tbb::parallel_reduce(tbb::blocked_range<int>(0, steps - 1), core, tbb::auto_partitioner());

    //    cout << "Sorting rms values... ("<< core.local_rms_list.size()<<" numbers)";
    //    rms_m** vect = new rms_m*[core.local_rms_list.size()];

    //    int num = 0;
    //    for (auto &k : core.local_rms_list){
    //        vect[num] = &k;
    //        if (k.i == 0 and k.j < 100) cout << k.j <<"  " << k.rms << endl;
    //        num++;
    //    }
    //
    //    tbb::parallel_sort(vect,vect+num, [](const rms_m *m1, const rms_m *m2){return (m1->rms < m2->rms);});
    //    for (int n =0 ; n< 100; n++){
    //        cout <<vect[n]->i << "  " <<  vect[n]->j <<"  " << vect[n]->rms << endl;
    //    }
    //    cout << vect[num-1]->i << "  " << vect[num-1]->j <<"  " <<vect[num-1]->rms << endl;
    //  core.local_rms_list.sort([](const rms_m &m1, const rms_m &m2){return (m1.rms < m2.rms);});


    //    for (int index1 = 0; index1 < 100; ++index1){
    //        cout << rmsvalue(0,index1) << endl;
    //    }

    //    for (int index1 = 0; index1 < steps - 1; ++index1) {
    //        for (int index2 = index1 + 1; index2 < steps; ++index2) {
    //            rms_list.emplace_back(index1, index2, rmsvalue(index1, index2));
    //        }
    //    }

    //    rms_list.sort([](const rms_m &m1, const rms_m &m2){return (m1.rms < m2.rms);});
    //   cout << "    Done" << endl;

    int n = 0;
    for (auto &v : core.local_rms_list) {
        cout << v.i << "  " << v.j << "  " << v.rms << endl;
        n++;
        if (n > 99) break;
    }
    vector<rms_m2> c;
    for (int i = 0; i < this->steps; i++) {
        c.emplace_back(i, i);
    }
    bool bChange;
    do {
        bChange = false;
        for (auto &k : core.local_rms_list) {
            //        for(int n = 0; n < num; n++){
            //            rms_m k = *(vect[n]);
            if (k.rms >= this->cutoff) break;
            int diff = c[k.j].clust - c[k.i].clust;
            if (diff) {
                bChange = true;
                if (diff > 0) {
                    c[k.j].clust = c[k.i].clust;
                } else {
                    c[k.i].clust = c[k.j].clust;
                }
            }
        }
    } while (bChange);
    //   delete [] vect;
    cout << "Sorting and renumbering clusters..." << endl;
    tbb::parallel_sort(c.begin(), c.end(), [](const rms_m2 &i, const rms_m2 &j) { return (i.clust < j.clust); });
    int cid = 1;
    unsigned int k;
    for (k = 1; k < c.size(); k++) {
        if (c[k].clust != c[k - 1].clust) {
            c[k - 1].clust = cid;
            cid++;
        } else {
            c[k - 1].clust = cid;
        }
    }
    c[k - 1].clust = cid;
    tbb::parallel_sort(c.begin(), c.end(), [](const rms_m2 &i, const rms_m2 &j) { return (i.conf < j.conf); });
    outfile << "***************************" << endl;
    outfile << "*Cluster Analysis(Linkage)*" << endl;
    outfile << "SET:" << ids << endl;
    outfile << "cutoff : " << this->cutoff << endl;
    outfile << "***************************" << endl;
    for (unsigned int k = 0; k < c.size(); k++) {
        outfile << k + 1 << "   " << c[k].clust << endl;
    }
    outfile << "***************************" << endl;

}

void Cluster::readInfo() {

    Atom::select1group(ids, "Please enter group:");
    this->cutoff = choose(0.0, static_cast<double>(numeric_limits<int>::max()), "Cutoff : ");

}

void Cluster::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->ids)) this->group.insert(atom);
                  });
}


enum class AminoAcidType {
    H3N_Ala, H3N_Gly, H3N_Ile, H3N_Leu, H3N_Pro, H3N_Val, H3N_Phe, H3N_Trp, H3N_Tyr, H3N_Asp,
    H3N_Glu, H3N_Arg, H3N_His, H3N_Lys, H3N_Ser, H3N_Thr, H3N_Cys, H3N_Met, H3N_Asn, H3N_Gln,
    Ala, Gly, Ile, Leu, Pro, Val, Phe, Trp, Tyr, Asp,
    Glu, Arg, His, Lys, Ser, Thr, Cys, Met, Asn, Gln
};


boost::bimap<boost::bimaps::set_of<AminoAcidType>, boost::bimaps::set_of<std::string>> aminotype_str_bimap =
        boost::assign::list_of<boost::bimap<boost::bimaps::set_of<AminoAcidType>, boost::bimaps::set_of<std::string>>::relation>
                (AminoAcidType::H3N_Ala, "H3N_Ala")
                (AminoAcidType::H3N_Gly, "H3N_Gly")
                (AminoAcidType::H3N_Pro, "H3N_Pro")
                (AminoAcidType::H3N_Trp, "H3N_Trp")
                (AminoAcidType::H3N_Asp, "H3N_Asp")
                (AminoAcidType::H3N_Glu, "H3N_Glu")
                (AminoAcidType::H3N_Arg, "H3N_Arg")

                (AminoAcidType::Ala, "Ala")
                (AminoAcidType::Gly, "Gly")
                (AminoAcidType::Pro, "Pro")
                (AminoAcidType::Trp, "Trp")
                (AminoAcidType::Asp, "Asp")
                (AminoAcidType::Glu, "Glu")
                (AminoAcidType::Arg, "Arg");

class AminoAcid {
public:
    AminoAcidType type;
    map<int, string> atom_no_map;
    int sequence_no;

};

class AminoTop {
public:
    class AminoItem {
    public:
        int no;
        Symbol symbol;
        std::list<int> linked_atom_nos;
        shared_ptr<Atom> atom;
        string H_;

    };

    std::map<int, std::shared_ptr<AminoItem>> topmap;

    void atom_null() {
        for (auto item : topmap) {
            item.second->atom.reset();
        }
    }

    AminoAcidType type;

    void readTop(const string &filename) {
        fstream f;
        f.open(filename);
        string line;
        while (true) {
            getline(f, line);
            auto field = split(line);
            if (field.empty()) break;
            auto item = make_shared<AminoItem>();
            item->no = stoi(field[0]);
            auto sym = field[1];
            if (sym == "H") item->symbol = Symbol::Hydrogen;
            else if (sym == "C") item->symbol = Symbol::Carbon;
            else if (sym == "N") item->symbol = Symbol::Nitrogen;
            else if (sym == "O") item->symbol = Symbol::Oxygen;
            else if (sym == "X") item->symbol = Symbol::X;


            item->H_ = field[2];

            for (unsigned int i = 3; i < field.size(); i++) {
                item->linked_atom_nos.push_back(stoi(field[i]));
            }
            topmap[item->no] = item;

        }

    }
};

class NMRRange : public BasicAnalysis {


    void recognize_amino_acid(std::shared_ptr<Frame> &frame);

    bool recognize_walk(shared_ptr<Atom> atom, shared_ptr<AminoTop::AminoItem> item,
                        AminoTop &top, shared_ptr<Frame> &frame,
                        list<shared_ptr<Atom>> &atom_list, list<shared_ptr<AminoTop::AminoItem>> &item_list);

    bool first_frame = true;


    void loadTop() {
        string path = std::getenv("ANALYSIS_TOP_PATH");
        for (auto &t : aminotype_str_bimap) {
            AminoTop top;
            top.readTop(path + "/" + t.right + ".top");
            top.type = t.left;
            amino_top_list.push_back(top);
        }

    }

    list<AminoTop> amino_top_list;

    list<shared_ptr<AminoAcid>> amino_acid_list;

    map<pair<int, int>, list<double>> dist_range_map;

    map<int, string> name_map;

public:
    NMRRange() {
        enable_outfile = true;
    }

    void process(std::shared_ptr<Frame> &frame) override;

    void print() override;

    void readInfo() override;

    static const string title() {
        return "NMRRange Analysis";
    }
};

void NMRRange::process(std::shared_ptr<Frame> &frame) {
    if (first_frame) {
        loadTop();
        recognize_amino_acid(frame);
        first_frame = false;
        vector<int> need_to_calc_list;
        for (auto &aminoacid : amino_acid_list) {
            for (auto &item : aminoacid->atom_no_map)
                need_to_calc_list.push_back(item.first);
        }
        for (size_t i = 0; i < need_to_calc_list.size() - 1; i++) {
            for (size_t j = i + 1; j < need_to_calc_list.size(); j++) {
                dist_range_map[make_pair(need_to_calc_list[i], need_to_calc_list[j])] = list<double>();
                dist_range_map[make_pair(need_to_calc_list[i], need_to_calc_list[j])].push_back(
                        atom_distance(frame->atom_map[need_to_calc_list[i]], frame->atom_map[need_to_calc_list[j]],
                                      frame));
            }
        }
    } else {
        for (auto &item : dist_range_map) {
            item.second.push_back(
                    atom_distance(frame->atom_map[item.first.first], frame->atom_map[item.first.second], frame));
        }
    }

}


void NMRRange::print() {

    for (auto &item : dist_range_map) {

        double sum = 0.0;
        int count = 0;
        for (auto v : item.second) {
            sum += pow(v, -6);
            count++;
        }
        double avg = sum / count;
        double e = -1.0 / 6;
        double value = pow(avg, e);
        outfile << name_map[item.first.first] << " <->\t"
                << name_map[item.first.second] << "\t" << value << endl;
    }
}

void NMRRange::readInfo() {
    enable_forcefield = true;

}


void NMRRange::recognize_amino_acid(std::shared_ptr<Frame> &frame) {
    cout << "First Frame, Analyze..." << endl;
    int amino_seq_no = 0;
    for (auto atom : frame->atom_list) {
        if (which(atom) == Symbol::Nitrogen) {
            for (auto &amino_top : amino_top_list) {
                list<shared_ptr<Atom>> child_atom_list;
                list<shared_ptr<AminoTop::AminoItem>> child_item_list;
                bool match = recognize_walk(atom, amino_top.topmap[1], amino_top, frame,
                                            child_atom_list, child_item_list);


                if (match) {
                    amino_seq_no++;
                    cout << aminotype_str_bimap.left.find(amino_top.type)->second << "   " << endl;

                    auto new_amino = make_shared<AminoAcid>();
                    new_amino->type = amino_top.type;
                    new_amino->sequence_no = amino_seq_no;

                    for (auto &item :amino_top.topmap) {
                        if (item.second->symbol == Symbol::Hydrogen) {
                            new_amino->atom_no_map[item.second->atom->seq] = item.second->H_;

                            name_map[item.second->atom->seq] = to_string(item.second->atom->seq) + ":"
                                                               + aminotype_str_bimap.left.find(amino_top.type)->second
                                                               + ":" + to_string(amino_seq_no) + " " + item.second->H_;
                        }
                        cout << item.first << " " << item.second->H_ <<
                             "," << item.second->atom->seq << " " << item.second->atom->atom_name << endl;
                        //if (item.second->symbol == Symbol::X)
                        item.second->atom->mark = false;
                    }
                    amino_acid_list.push_back(new_amino);
                    amino_top.atom_null();
                    break;
                }


            }
        }
    }
    cout << endl;

}

bool NMRRange::recognize_walk(shared_ptr<Atom> atom, shared_ptr<AminoTop::AminoItem> item,
                              AminoTop &top, shared_ptr<Frame> &frame,
                              list<shared_ptr<Atom>> &atom_list,
                              list<shared_ptr<AminoTop::AminoItem>> &item_list) {
    if (item->symbol == Symbol::X) {
        atom->mark = true;
        item->atom = atom;
        atom_list.push_back(atom);
        item_list.push_back(item);
        return true;
    }
    if (which(atom) != item->symbol) {
        return false;
    }
    if (atom->con_list.size() != item->linked_atom_nos.size()) {

        return false;
    }
    atom->mark = true;
    item->atom = atom;

    bool match = true;
    list<shared_ptr<Atom>> child_atom_list;
    list<shared_ptr<AminoTop::AminoItem>> child_item_list;
    for (int i : atom->con_list) {
        auto next_atom = frame->atom_map[i];
        if (next_atom->mark) continue;
        bool ok = false;

        for (int j : item->linked_atom_nos) {
            auto next_item = top.topmap[j];
            if (next_item->atom) continue;
            ok = recognize_walk(next_atom, next_item, top, frame, child_atom_list, child_item_list);
            if (ok) break;
        }

        match = ok && match;
        if (!match) break;
    }
    if (match) {
        atom_list.push_back(atom);
        item_list.push_back(item);
        for (auto &child_atom: child_atom_list) {
            atom_list.push_back(child_atom);
        }
        for (auto &child_item : child_item_list) {
            item_list.push_back(child_item);
        }

        return true;
    }

    atom->mark = false;
    item->atom.reset();

    for (auto &child_atom: child_atom_list) {
        child_atom->mark = false;
    }
    for (auto &child_item : child_item_list) {
        child_item->atom.reset();
    }

    return false;
}

class RotAcfCutoff : public BasicAnalysis {

public:

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    explicit RotAcfCutoff() {
        enable_outfile = true;
        enable_forcefield = true;
    }

    void process(std::shared_ptr<Frame> &frame) override;

    void print() override;

    void readInfo() override;

    static const string title() {
        return "Rotational autocorrelation function (Cutoff)";
    }

    struct InnerAtom {
        int index;
        std::list<std::tuple<double, double, double>> *list_ptr = nullptr;

        InnerAtom(int index, std::list<std::tuple<double, double, double>> *list_ptr)
                : index(index), list_ptr(list_ptr) {}
    };

    struct InnerAtomHasher {
        typedef InnerAtom argument_type;
        typedef std::size_t result_type;

        result_type operator()(argument_type const &s) const noexcept {
            result_type seed = 0;
            boost::hash_combine(seed, s.index);
            boost::hash_combine(seed, s.list_ptr);
            return seed;
        }
    };


private:

    double time_increment_ps = 0.1;
    double cutoff2;

    Atom::AtomIndenter ids1;
    Atom::AtomIndenter ids2;

    std::unordered_set<shared_ptr<Atom>> group1;
    std::unordered_set<shared_ptr<Atom>> group2;

    std::unordered_set<InnerAtom, InnerAtomHasher> inner_atoms;

    std::list<std::list<std::tuple<double, double, double>> *> rots;

    auto find_in(int seq);

    tuple<double, double, double> calVector(shared_ptr<Molecule> &mol, shared_ptr<Frame> &frame);


};

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

    double cutoff = choose(0.0, GMX_DOUBLE_MAX, "Please enter distance cutoff:");
    this->cutoff2 = cutoff * cutoff;

    this->time_increment_ps = choose(0.0, GMX_DOUBLE_MAX, "Enter the Time Increment in Picoseconds [0.1]:", true, 0.1);

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

void RotAcfCutoff::print() {
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

    outfile << "*********************************************************" << endl;
    outfile << "cutoff : " << std::sqrt(cutoff2) << std::endl;
    outfile << "First Type : " << ids1 << " Second Type : " << ids2 << endl;
    outfile << " rotational autocorrelation function" << endl;

    outfile << "    Time Gap      ACF       intergrate" << endl;
    outfile << "      (ps)                    (ps)" << endl;

    for (std::size_t t = 0; t < acf.size(); t++) {
        outfile << boost::format("%12.2f%18.14f%15.5f") % (t * time_increment_ps) % acf[t].second % integrate[t]
                << endl;
    }
    outfile << "*********************************************************" << endl;

}

void RotAcfCutoff::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->ids1)) this->group1.insert(atom);
                      if (Atom::is_match(atom, this->ids2)) this->group2.insert(atom);
                  });
}


class DiffuseCutoff : public BasicAnalysis {

public:

    DiffuseCutoff() {
        enable_outfile = true;
        enable_forcefield = true;
    }

    void process(std::shared_ptr<Frame> &frame) override;

    void print() override;

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void readInfo() override;

    static const string title() {
        return "Diffusion Cutoff Coefficient by Einstein equation";
    }

    struct InnerAtom {
        int index;
        std::list<std::tuple<double, double, double>> *list_ptr = nullptr;

        InnerAtom(int index, std::list<std::tuple<double, double, double>> *list_ptr)
                : index(index), list_ptr(list_ptr) {}
    };

    struct InnerAtomHasher {
        typedef InnerAtom argument_type;
        typedef std::size_t result_type;

        result_type operator()(argument_type const &s) const noexcept {
            result_type
            const h1(std::hash<int>()
            (s.index));
            result_type
            const h2(std::hash<std::list<std::tuple<double, double, double>> *>()
            (s.list_ptr));
            return h1 ^ (h2 << 1); // or use boost::hash_combine (see Discussion)
        }
    };

private:

    double time_increment_ps = 0.1;
    double cutoff2;

    Atom::AtomIndenter ids1;
    Atom::AtomIndenter ids2;

    std::unordered_set<shared_ptr<Atom>> group1;
    std::unordered_set<shared_ptr<Atom>> group2;

    std::unordered_set<InnerAtom, InnerAtomHasher> inner_atoms;

    std::list<std::list<std::tuple<double, double, double>> *> rcm;

    auto find_in(int seq);

};

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
        outfile << fmt::sprintf("%12.2f%12.2f%12.2f%12.2f%12.2f%12.4f",
                                delta, xvalue, yvalue, zvalue, rvalue, dvalue) << endl;
    }
    outfile << "*********************************************************" << endl;
}

void DiffuseCutoff::readInfo() {

    Atom::select2group(ids1, ids2);

    double cutoff = choose(0.0, GMX_DOUBLE_MAX, "Please enter distance cutoff:");
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


class Diffuse : public BasicAnalysis {

public:
    Diffuse() {
        enable_forcefield = true;
        enable_outfile = true;
    }

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print() override;

    void readInfo() override;

    static const string title() {
        return "Diffusion Coefficient by Einstein equation";
    }

private:

    Atom::AtomIndenter ids;

    std::unordered_set<shared_ptr<Atom>> group;
    bool first_round = true;

    bool first_frame = true;
    int total_frame_number;
    int steps = 0;
    double time_increment_ps = 0.1;
    int total_mol = 0;


    bool bSerial = true;
    bool bTradition = true;

    Eigen::MatrixXd xcm;
    Eigen::MatrixXd ycm;
    Eigen::MatrixXd zcm;

};

void Diffuse::process(std::shared_ptr<Frame> &frame) {

    if (first_frame) {
        first_frame = false;
        for (auto &mol : frame->molecule_list) {
            mol->bExculde = true;
        }

        for (auto &atom : group) {
            auto mol = atom->molecule.lock();
            if (mol) {
                mol->bExculde = false;
            }
        }

        for (auto &mol : frame->molecule_list) {
            if (mol->bExculde) continue;
            mol->calc_mass();
            total_mol++;
        }
        cout << "Total molecule number : " << total_mol << endl;
        xcm = Eigen::MatrixXd::Zero(total_frame_number, total_mol);
        ycm = Eigen::MatrixXd::Zero(total_frame_number, total_mol);
        zcm = Eigen::MatrixXd::Zero(total_frame_number, total_mol);

        int mol_index = 0;
        for (auto &mol: frame->molecule_list) {
            if (!mol->bExculde) {
                auto coord = mol->calc_weigh_center(frame);
                xcm(0, mol_index) = get<0>(coord);
                ycm(0, mol_index) = get<1>(coord);
                zcm(0, mol_index) = get<2>(coord);
                mol_index++;
            }
        }
    } else {
        int mol_index = 0;
        for (auto &mol: frame->molecule_list) {
            if (!mol->bExculde) {
                auto coord = mol->calc_weigh_center(frame);
                double xold = xcm(steps - 1, mol_index);
                double yold = ycm(steps - 1, mol_index);
                double zold = zcm(steps - 1, mol_index);
                double xr = get<0>(coord) - xold;
                double yr = get<1>(coord) - yold;
                double zr = get<2>(coord) - zold;
                frame->image(xr, yr, zr);
                xcm(steps, mol_index) = xr + xold;
                ycm(steps, mol_index) = yr + yold;
                zcm(steps, mol_index) = zr + zold;
                mol_index++;
            }
        }
    }
    steps++;

}

void Diffuse::print() {

    vector<int> ntime(total_frame_number, 0);
    vector<double> xmsd(total_frame_number, 0.0);
    vector<double> ymsd(total_frame_number, 0.0);
    vector<double> zmsd(total_frame_number, 0.0);

    if (bSerial) {
        if (bTradition) {
            for (int i = 0; i < total_frame_number - 1; i++) {
                for (int j = i + 1; j < total_frame_number; j++) {
                    int m = j - i - 1;
                    ntime[m]++;
                    for (int k = 0; k < total_mol; k++) {
                        double xdiff = xcm(j, k) - xcm(i, k);
                        double ydiff = ycm(j, k) - ycm(i, k);
                        double zdiff = zcm(j, k) - zcm(i, k);
                        xmsd[m] += xdiff * xdiff;
                        ymsd[m] += ydiff * ydiff;
                        zmsd[m] += zdiff * zdiff;

                    }
                }
            }
        } else {
            for (int m = 0; m < total_frame_number - 1; m++) {
                for (int i = 0; i < total_frame_number - 1; i += m + 1) {
                    int j = i + m + 1;
                    if (j < total_frame_number) {
                        ntime[m]++;
                        for (int k = 0; k < total_mol; k++) {
                            double xdiff = xcm(j, k) - xcm(i, k);
                            double ydiff = ycm(j, k) - ycm(i, k);
                            double zdiff = zcm(j, k) - zcm(i, k);
                            xmsd[m] += xdiff * xdiff;
                            ymsd[m] += ydiff * ydiff;
                            zmsd[m] += zdiff * zdiff;
                        }
                    }
                }
            }
        }

    } else {

        class Body {
        public:
            int total_frame_number;
            vector<int> ntime;
            vector<double> xmsd;
            vector<double> ymsd;
            vector<double> zmsd;

            int total_mol;
            Eigen::MatrixXd &xcm;
            Eigen::MatrixXd &ycm;
            Eigen::MatrixXd &zcm;

            Body(int total_frame_number, int total_mol, Eigen::MatrixXd &xcm, Eigen::MatrixXd &ycm,
                 Eigen::MatrixXd &zcm) :
                    total_frame_number(total_frame_number),
                    ntime(total_frame_number - 1), xmsd(total_frame_number - 1),
                    ymsd(total_frame_number - 1), zmsd(total_frame_number - 1),
                    total_mol(total_mol), xcm(xcm), ycm(ycm), zcm(zcm) {}

            Body(const Body &body, tbb::split) :
                    total_frame_number(body.total_frame_number),
                    ntime(body.total_frame_number - 1), xmsd(body.total_frame_number - 1),
                    ymsd(body.total_frame_number - 1), zmsd(body.total_frame_number - 1),
                    total_mol(body.total_mol), xcm(body.xcm), ycm(body.ycm), zcm(body.zcm) {}

            void join(const Body &body) {
                for (int i = 0; i < total_frame_number - 1; i++) {
                    ntime[i] += body.ntime[i];
                    xmsd[i] += body.xmsd[i];
                    ymsd[i] += body.ymsd[i];
                    zmsd[i] += body.zmsd[i];
                }

            }

            void operator()(const tbb::blocked_range<int> &range) {
                for (int i = range.begin(); i != range.end(); i++) {
                    for (int j = i + 1; j < total_frame_number; j++) {
                        int m = j - i - 1;
                        ntime[m]++;
                        for (int k = 0; k < total_mol; k++) {
                            double xdiff = xcm(j, k) - xcm(i, k);
                            double ydiff = ycm(j, k) - ycm(i, k);
                            double zdiff = zcm(j, k) - zcm(i, k);
                            xmsd[m] += xdiff * xdiff;
                            ymsd[m] += ydiff * ydiff;
                            zmsd[m] += zdiff * zdiff;

                        }
                    }
                }
            }
        } body(total_frame_number, total_mol, xcm, ycm, zcm);

        tbb::parallel_reduce(tbb::blocked_range<int>(0, total_frame_number - 1), body, tbb::auto_partitioner());


        ntime = body.ntime;
        xmsd = body.xmsd;
        ymsd = body.ymsd;
        zmsd = body.zmsd;
    }

    const double dunits = 10.0;


    for (int i = 0; i < total_frame_number - 1; i++) {
        double counts = total_mol * ntime[i];
        xmsd[i] /= counts;
        ymsd[i] /= counts;
        zmsd[i] /= counts;
    }

    outfile << "*********************************************************" << endl;
    outfile << "Group :" << ids << endl;
    outfile << "Mean Squared Displacements and Self-Diffusion Constant" << endl;
    outfile << "    Time Gap      X MSD       Y MSD       Z MSD       R MSD        Diff Const" << endl;
    outfile << "      (ps)       (Ang^2)     (Ang^2)     (Ang^2)     (Ang^2)     (x 10^-5 cm**2/sec)" << endl;

    for (int i = 0; i < total_frame_number - 1; i++) {
        double delta = time_increment_ps * (i + 1);
        double xvalue = xmsd[i];
        double yvalue = ymsd[i];
        double zvalue = zmsd[i];
        double rvalue = xmsd[i] + ymsd[i] + zmsd[i];
        double dvalue = dunits * rvalue / delta / 6.0;
        outfile << fmt::sprintf("%12.2f%12.2f%12.2f%12.2f%12.2f%12.4f",
                                delta, xvalue, yvalue, zvalue, rvalue, dvalue) << endl;

    }
    outfile << "*********************************************************" << endl;

}

void Diffuse::readInfo() {
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
    while (true) {
        string input_line = input(" Enter the Total Frame Number: ");
        boost::trim(input_line);
        if (!input_line.empty()) {
            total_frame_number = stoi(input_line);
            if (total_frame_number <= 0) {
                cerr << "error total frame number " << total_frame_number << endl;
                continue;
            }
            break;
        }
    }
    Atom::select1group(ids, "select group :");

    while (true) {
        string input_line = input(" serial ? [T]:");
        boost::trim(input_line);
        if (!input_line.empty()) {
            if (input_line.compare("T") == 0) {
                bSerial = true;
                while (true) {
                    string input_line = input(" trandition ? [T]:");
                    boost::trim(input_line);
                    if (!input_line.empty()) {
                        if (input_line.compare("T") == 0) {
                            bTradition = true;
                            break;
                        } else if (input_line.compare("F") == 0) {
                            bTradition = false;
                            break;
                        }
                    }
                }
                break;
            } else if (input_line.compare("F") == 0) {
                bSerial = false;
                enable_tbb = true;
                break;
            }
        }
    }


}

void Diffuse::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->ids)) this->group.insert(atom);
                  });
}


class DipoleAngle : public BasicAnalysis {
public:
    DipoleAngle() {
        enable_forcefield = true;
        enable_outfile = true;
    }

    virtual ~DipoleAngle() = default;

    void process(std::shared_ptr<Frame> &frame) override;

    void print() override;

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void readInfo() override;

    static const string title() {
        return "Dipole Angle";
    }

protected:

    Atom::AtomIndenter ids1;
    Atom::AtomIndenter ids2;

    std::unordered_set<shared_ptr<Atom>> group1;
    std::unordered_set<shared_ptr<Atom>> group2;


    double distance_width;
    double angle_width;

    int distance_bins;
    int angle_bins;

    map<pair<int, int>, size_t> hist;
};

void DipoleAngle::process(std::shared_ptr<Frame> &frame) {

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
            mol->ow = atom;
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

            double x1 = mol->ow->x;
            double y1 = mol->ow->y;
            double z1 = mol->ow->z;

            double xr = x1 - ref_x;
            double yr = y1 - ref_y;
            double zr = z1 - ref_z;

            frame->image(xr, yr, zr);

            double distance = sqrt(xr * xr + yr * yr + zr * zr);

            auto dipole = mol->calc_dipole(frame);

            double dipole_scalar = sqrt(pow(get<0>(dipole), 2) + pow(get<1>(dipole), 2) + pow(get<2>(dipole), 2));

            double _cos =
                    (xr * get<0>(dipole) + yr * get<1>(dipole) + zr * get<2>(dipole)) / (distance * dipole_scalar);

            double angle = acos(_cos) * 180.0 / 3.1415926;

            int i_distance_bin = int(distance / distance_width) + 1;
            int i_angle_bin = int(angle / angle_width) + 1;

            if (i_distance_bin <= distance_bins and i_angle_bin <= angle_bins) {
                hist[make_pair(i_distance_bin, i_angle_bin)] += 1;
            }
        }
    }
}

void DipoleAngle::print() {

//    outfile << boost::format("%5s") % "Dist";
//    for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
//        outfile << boost::format("%7.1f") % (i_angle  * angle_width);
//    }
//    outfile << endl;
//    for (int i_distance = 1; i_distance < distance_bins; i_distance++) {
//        outfile << boost::format("%5.3f") % (i_distance * distance_width);
//        size_t total = 0;
//        for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
//            total += hist[make_pair(i_distance,i_angle)];
//        }
//        for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
//            outfile << boost::format("%7.4f") % ( total == 0 ? 0.0 : double(hist[make_pair(i_distance,i_angle)]) / total);
//        }
//        outfile << endl;
//    }
    for (int i_distance = 1; i_distance < distance_bins; i_distance++) {
        size_t total = 0;
        for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
            total += hist[make_pair(i_distance, i_angle)];
        }
        for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
            outfile << boost::format("%5.3f%10.3f%10.4f\n")
                       % ((i_distance - 0.5) * distance_width)
                       % ((i_angle - 0.5) * angle_width)
                       % (total == 0 ? 0.0 : double(hist[make_pair(i_distance, i_angle)]) / total);
        }
    }
}

void DipoleAngle::readInfo() {
    Atom::select2group(ids1, ids2);
    double rmax = choose(0.0, GMX_DOUBLE_MAX, "Enter Maximum Distance to Accumulate[10.0 Ang]:", true, 10.0);
    distance_width = choose(0.0, GMX_DOUBLE_MAX, "Enter Width of Distance Bins [0.01 Ang]:", true, 0.01);
    double angle_max = choose(0.0, 180.0, "Enter Maximum Angle to Accumulate[180.0 degree]:", true, 180.0);
    angle_width = choose(0.0, 180.0, "Enter Width of Angle Bins [0.5 degree]:", true, 0.5);

    distance_bins = int(rmax / distance_width);
    angle_bins = int(angle_max / angle_width);

    for (int i_distance = 1; i_distance <= distance_bins; i_distance++) {
        for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
            hist[make_pair(i_distance, i_angle)] = 0;
        }
    }
}

void DipoleAngle::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->ids1)) this->group1.insert(atom);
                      if (Atom::is_match(atom, this->ids2)) this->group2.insert(atom);
                  });
}

class DipoleAngleSingleDistanceNormal : public DipoleAngle {
public:
    void print() override;

    static const string title() {
        return "Dipole Angle of single distance normal";
    };

};

void DipoleAngleSingleDistanceNormal::print() {
    outfile << "DipoleAngleSingleDistanceNormal : " << std::endl;
    double factor = (4.0 / 3.0) * M_PI;
    for (int i_distance = 1; i_distance < distance_bins; i_distance++) {
        size_t total = 0;
        for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
            total += hist[make_pair(i_distance, i_angle)];
        }
        double dv = factor * (pow(i_distance * distance_width, 3) - pow((i_distance - 1) * distance_width, 3));
        for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
            outfile << boost::format("%5.3f   %10.2f   %g\n")
                       % ((i_distance - 0.5) * distance_width)
                       % ((i_angle - 0.5) * angle_width)
                       % (total == 0 ? 0.0 : double(hist[make_pair(i_distance, i_angle)]) / (total * dv * angle_width));
        }
    }
}

class DipoleAngleVolumeNormal : public DipoleAngle {
public:
    void print() override;

    static const string title() {
        return "Dipole Angle of volume normal";
    }
};

void DipoleAngleVolumeNormal::print() {
    outfile << "DipoleAngleVolumeNormal : " << std::endl;
    size_t total = 0;
    for (int i_distance = 1; i_distance < distance_bins; i_distance++) {
        for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
            total += hist[make_pair(i_distance, i_angle)];
        }

    }
    double factor = (4.0 / 3.0) * M_PI;
    for (int i_distance = 1; i_distance < distance_bins; i_distance++) {
        double dv = factor * (pow(i_distance * distance_width, 3) - pow((i_distance - 1) * distance_width, 3));
        for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
            outfile << boost::format("%5.3f   %10.3f   %g\n")
                       % ((i_distance - 0.5) * distance_width)
                       % ((i_angle - 0.5) * angle_width)
                       % (total == 0 ? 0.0 : double(hist[make_pair(i_distance, i_angle)]) / (total * dv * angle_width));
        }
    }
}


class DipoleAngle2Gibbs : public DipoleAngle {
public:
    void print() override;

    void readInfo() override;

    static const string title() {
        return "Dipole Angle to Gibbs Free Energy";
    };

protected:

    const double kb = 1.380649e-23; // unit: J/K
    double temperature;  // unit: K
    const double avogadro_constant = 6.022140857e23;

};

void DipoleAngle2Gibbs::print() {

    double factor = -kb * temperature * avogadro_constant / 4184.0;
    double max_value = 0.0;

    for (int i_distance = 1; i_distance < distance_bins; i_distance++) {
        double dv = pow(i_distance * distance_width, 3) - pow((i_distance - 1) * distance_width, 3);
        for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
            max_value = max(max_value, hist[make_pair(i_distance, i_angle)] / (dv));
        }
    }
    for (int i_distance = 1; i_distance < distance_bins; i_distance++) {
        double dv = pow(i_distance * distance_width, 3) - pow((i_distance - 1) * distance_width, 3);
        for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
            double pop = double(hist[make_pair(i_distance, i_angle)]) / (max_value * dv);
            outfile << boost::format("%5.3f   %10.3f   %15.6f\n")
                       % ((i_distance - 0.5) * distance_width)
                       % ((i_angle - 0.5) * angle_width)
                       % (pop == 0.0 ? 100.0 : factor * log(pop));
            //               %  pop;
        }
    }

    /*for (int i_distance = 1; i_distance < distance_bins; i_distance++) {
        for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
            max_value = max(max_value,hist[make_pair(i_distance, i_angle)] / (angle_width));
        }
    }
    for (int i_distance = 1; i_distance < distance_bins; i_distance++) {
        double dv = pow(i_distance * distance_width, 3) - pow((i_distance - 1) * distance_width, 3);
        for (int i_angle = 1; i_angle <= angle_bins; i_angle++) {
            double pop = double(hist[make_pair(i_distance, i_angle)]) / (max_value * angle_width);
            outfile << boost::format("%5.3f%7.1f%10.6f\n")
                       % (i_distance * distance_width)
                       % (i_angle * angle_width)
                     //  % (pop == 0.0 ? 10000.0 : factor * log(pop));
                       %  pop;
        }
    }
    */
}

void DipoleAngle2Gibbs::readInfo() {
    DipoleAngle::readInfo();
    temperature = choose(0.0, 10000.0, "Temperature [298] (K):", true, 298.0);
}


class ShellDensity : public BasicAnalysis {
public:
    ShellDensity() { enable_outfile = true; }

    void process(std::shared_ptr<Frame> &frame) override;

    void print() override;

    void readInfo() override;

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    static const string title() {
        return "Shell Density function";
    }

protected:
    Atom::AtomIndenter ids1;
    Atom::AtomIndenter ids2;

    std::unordered_set<shared_ptr<Atom>> group1;
    std::unordered_set<shared_ptr<Atom>> group2;

    double distance_width;

    int distance_bins;

    map<int, size_t> hist;

    size_t nframe = 0;

};

void ShellDensity::process(std::shared_ptr<Frame> &frame) {
    nframe++;
    for (auto &ref :group1) {
        for (auto &atom : group2) {
            int ibin = int(atom_distance(ref, atom, frame) / distance_width) + 1;
            if (ibin <= distance_bins) {
                hist[ibin]++;
            }
        }

    }
}

void ShellDensity::print() {
    outfile << "************************************************" << endl;
    outfile << "***** Shell Density Function ****" << endl;

    outfile << "First Type : " << ids1 << " Second Type : " << ids2 << endl;

    outfile << "************************************************" << endl;
    outfile << "Bin    Distance    Densitry (count / Ang3 / frame)" << endl;

    for (int i = 1; i <= distance_bins; i++) {
        double dv = (4.0 / 3.0) * M_PI * (pow(i * distance_width, 3) - pow((i - 1) * distance_width, 3));
        outfile << fmt::sprintf("%d      %.4f      %g \n",
                                i, (i - 0.5) * distance_width, hist[i] / (nframe * dv));
    }

    outfile << "************************************************" << endl;
}


void ShellDensity::readInfo() {
    Atom::select2group(ids1, ids2);
    double rmax = choose(0.0, GMX_DOUBLE_MAX, "Enter Maximum Distance to Accumulate[10.0 Ang]:", true, 10.0);
    distance_width = choose(0.0, GMX_DOUBLE_MAX, "Enter Width of Distance Bins [0.01 Ang]:", true, 0.01);
    distance_bins = int(rmax / distance_width);
    for (int i = 1; i <= distance_bins; i++) {
        hist[i] = 0;
    }
}

void ShellDensity::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->ids1)) this->group1.insert(atom);
                      if (Atom::is_match(atom, this->ids2)) this->group2.insert(atom);
                  });
}

class SearchInteractionResidue : public BasicAnalysis {
public:
    SearchInteractionResidue() { enable_outfile = true; }

    void process(std::shared_ptr<Frame> &frame) override;

    void print() override;

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void readInfo() override;

    static const string title() { return "Search Interaction Residue between two groups"; }

private:
    Atom::AtomIndenter ids1;
    Atom::AtomIndenter ids2;

    std::unordered_set<shared_ptr<Atom>> group1;
    std::unordered_set<shared_ptr<Atom>> group2;

    double cutoff;

    std::list<std::unordered_set<std::string>> interaction_residues;
    int total_frames = 0;

    enum class OutputStyle {
        BOOL = 0, NUMBER = 1
    };
    OutputStyle style;

};

void SearchInteractionResidue::process(std::shared_ptr<Frame> &frame) {

    std::unordered_set<std::string> residue_set; // resname:no

    for (auto &atomA : group1) {
        for (auto &atomB : group2) {
            if (atom_distance(atomA, atomB, frame) <= cutoff) {
                residue_set.insert(
                        atomB->residue_name.get() + "-" + boost::lexical_cast<std::string>(atomB->residue_num.get()));
            }
        }
    }
    interaction_residues.push_back(residue_set);
    total_frames++;
}


void SearchInteractionResidue::print() {
    outfile << "************************************************\n";
    outfile << "*****" << SearchInteractionResidue::title() << " ****\n";

    outfile << "First Group  " << ids1 << " Second Group  " << ids2 << '\n';
    outfile << "cutoff :" << cutoff << " Ang\n";

    outfile << "************************************************\n";


    using ResItem = struct {
        std::string name;
        int count;
    };

    std::unordered_map<std::string, ResItem *> map;
    for (auto &set : interaction_residues) {
        for (auto &item : set) {
            auto it = map.find(item);
            if (it == map.end()) {
                map[item] = new ResItem{item, 1};
            } else {
                it->second->count++;
            }
        }
    }

    std::vector<ResItem *> itemVec;
    for (auto &item : map) {
        itemVec.push_back(item.second);
    }

    std::sort(itemVec.begin(), itemVec.end(), [](ResItem *i1, ResItem *i2) { return i1->count > i2->count; });

    std::size_t nframe = 1;
    outfile << boost::format("%10s") % "name";
    for (auto &item : itemVec) {
        outfile << boost::format("%10s") % item->name;
    }
    outfile << boost::format("\n%10s") % "Freq";;
    for (auto &item : itemVec) {
        outfile << boost::format("%9.1f%%") % (item->count * 100.0 / total_frames);
    }
    outfile << '\n';
    for (auto &set : interaction_residues) {
        outfile << boost::format("%10d") % nframe;
        int index = 1;
        for (auto &item : itemVec) {
            outfile << boost::format("%10d") % (style == OutputStyle::NUMBER ?
                                                (set.count(item->name) ? index : 0) : set.count(item->name));
            index++;
        }
        outfile << '\n';
        nframe++;
    }

    outfile << "************************************************" << endl;


}

void SearchInteractionResidue::readInfo() {
    std::cout << "The output residues is in the second group\n";
    Atom::select2group(ids1, ids2);
    cutoff = choose(0.0, GMX_DOUBLE_MAX, "cutoff [Ang]:", false);
    std::cout << "(0) bool style\n(1) number style\n";
    style = static_cast<OutputStyle>(choose(0, 1, "which style ? [ 0 ] ", true, 0));
}

void SearchInteractionResidue::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->ids1)) this->group1.insert(atom);
                      if (Atom::is_match(atom, this->ids2)) this->group2.insert(atom);
                  });
}


class FindMinBetweenTwoGroups : public BasicAnalysis {
public:
    FindMinBetweenTwoGroups() { enable_outfile = true; }

    void process(std::shared_ptr<Frame> &frame) override;

    void print() override;

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void readInfo() override;

    static const string title() { return "Find Min distance between two groups"; }

private:
    Atom::AtomIndenter ids;

    int total_frames = 0;

    std::vector<std::shared_ptr<Molecule>> mol_list;

    std::list<std::vector<double>> results;
};


void FindMinBetweenTwoGroups::process(std::shared_ptr<Frame> &frame) {

    std::size_t length = mol_list.size();

    std::vector<double> line_rest;
    for (std::size_t i = 0; i < length - 1; i++) {
        for (std::size_t j = i + 1; j < length; j++) {
            auto &mol1 = mol_list[i];
            auto &mol2 = mol_list[j];
            line_rest.push_back(min_distance(mol1, mol2, frame));
        }
    }
    results.push_back(line_rest);
}

void FindMinBetweenTwoGroups::print() {
    outfile << "************************************************\n";
    outfile << "*****" << FindMinBetweenTwoGroups::title() << " ****\n";
    outfile << "Group  " << ids << '\n';
    outfile << "************************************************\n";


    outfile << boost::format("%6s") % "Frame";
    std::size_t length = mol_list.size();
    for (std::size_t i = 0; i < length - 1; i++) {
        for (std::size_t j = i + 1; j < length; j++) {
            outfile << boost::format("  %4d-%-4d  ") % i % j;
        }
    }

    outfile << '\n';

    std::size_t nframe = 0;

    for (auto &v : results) {
        nframe++;

        outfile << boost::format("%6d") % nframe;
        for (auto value: v) outfile << boost::format("  %9.2f  ") % value;
        outfile << '\n';
    }


    outfile << "************************************************\n";
}


void FindMinBetweenTwoGroups::readInfo() {
    Atom::select1group(ids, "Input Residue Name Mask: ");
}

void FindMinBetweenTwoGroups::processFirstFrame(std::shared_ptr<Frame> &frame) {
    for (auto &mol : frame->molecule_list) {
        for (auto &atom : mol->atom_list) {
            if (Atom::is_match(atom, this->ids)) {
                mol_list.push_back(mol);
                break;
            }
        }
    }
}


class DemixIndexOfTwoGroup : public BasicAnalysis {
public:
    DemixIndexOfTwoGroup() { enable_outfile = true; }

    void process(std::shared_ptr<Frame> &frame) override;

    void print() override;

    void readInfo() override;

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    static const string title() { return "Calculate demix index of two groups"; }

private:

    auto calculate_grid_index(const std::shared_ptr<Atom> &atom, const std::shared_ptr<Frame> &frame);

    Atom::AtomIndenter ids1;
    Atom::AtomIndenter ids2;

    std::unordered_set<shared_ptr<Atom>> group1;
    std::unordered_set<shared_ptr<Atom>> group2;

    int grid_x;
    int grid_y;
    int grid_z;

    std::list<std::tuple<double, double>> demix_index_list;
};


auto
DemixIndexOfTwoGroup::calculate_grid_index(const std::shared_ptr<Atom> &atom, const std::shared_ptr<Frame> &frame) {

    auto box_index_x = int(atom->x / (frame->a_axis / grid_x)) % grid_x;
    auto box_index_y = int(atom->y / (frame->b_axis / grid_y)) % grid_y;
    auto box_index_z = int(atom->z / (frame->c_axis / grid_z)) % grid_z;


    if (!atom->mass) {
        std::cerr << "ERROR !!  Atom mass not available !\n";
        exit(EXIT_FAILURE);
    }

    assert(box_index_x >= 0 && box_index_x < grid_x);
    assert(box_index_y >= 0 && box_index_y < grid_y);
    assert(box_index_z >= 0 && box_index_z < grid_z);
    return make_tuple(box_index_x, box_index_y, box_index_z);
}

void DemixIndexOfTwoGroup::process(std::shared_ptr<Frame> &frame) {


    double group1_dens[grid_x][grid_y][grid_z];
    double group2_dens[grid_x][grid_y][grid_z];

    bzero(group1_dens, grid_x * grid_y * grid_z * sizeof(double));
    bzero(group2_dens, grid_x * grid_y * grid_z * sizeof(double));

    for (auto &atom : group1) {
        auto[box_index_x, box_index_y, box_index_z] =  calculate_grid_index(atom, frame);
        group1_dens[box_index_x][box_index_y][box_index_z] += atom->mass.get();
    }
    for (auto &atom : group2) {
        auto[box_index_x, box_index_y, box_index_z] =  calculate_grid_index(atom, frame);
        group2_dens[box_index_x][box_index_y][box_index_z] += atom->mass.get();
    }

    double d_sum = 0.0;
    double d1_sum = 0.0;
    double d2_sum = 0.0;
    for (int i = 0; i < grid_x; i++) {
        for (int j = 0; j < grid_y; j++) {
            for (int k = 0; k < grid_z; k++) {
                auto d1 = group1_dens[i][j][k];
                d1_sum += d1;
                auto d2 = group2_dens[i][j][k];
                d2_sum += d2;
                if (d1 == 0.0 or d2 == 0.0) continue;
                d_sum += 1 / ((1 / d1) + (1 / d2));
            }
        }
    }


    d_sum /= frame->volume();

    auto d_ideal = 1 / (1 / d1_sum + 1 / d2_sum);
    d_ideal /= frame->volume();

    demix_index_list.emplace_back(d_sum, d_ideal);
}

void DemixIndexOfTwoGroup::readInfo() {
    Atom::select2group(ids1, ids2, "Please select group1 > ", "Please select group2 > ");
    grid_x = choose<int>(0, 1000, "Grid in X dememsion  :  ");
    grid_y = choose<int>(0, 1000, "Grid in Y dememsion  :  ");
    grid_z = choose<int>(0, 1000, "Grid in Z dememsion  :  ");
}

void DemixIndexOfTwoGroup::print() {

//    constexpr double AvogadroConstant =  6.02214076E23;
//    constexpr double AngstromToCentimeter = 1E-8;
//
//    const double Unit = 1 / (pow(AngstromToCentimeter, 3) * AvogadroConstant);

    outfile << "#####################################\n";
    outfile << "#    Demix Index (normalization)\n";
    outfile << "#    Group1  " << ids1 << '\n';
    outfile << "#    Group2  " << ids2 << '\n';
    outfile << "#    Grid    X = " << grid_x << "  Y = " << grid_y << "  Z = " << grid_z << '\n';
    outfile << "#####################################\n";

    outfile << "@   title \"Demix Index\"\n";
    outfile << "@    xaxis  label \"Frame Number\"\n";
    outfile << "@    yaxis  label \"Demix Index\"\n";
    outfile << "@TYPE xy\n";
    outfile << "@ legend on\n";
    outfile << "@ legend length 1\n";
    outfile << "@ s0 legend \"Demix Index\"\n";
//    outfile << "@ s1 legend \"Demix:Ideal\"\n";
    outfile << "# Frame      Demix Index \n";
    int cyc = 1;
    for (auto[d, d_ideal] : demix_index_list) {
        outfile << cyc << "        " << d / d_ideal << '\n';
        cyc++;
    }


}

void DemixIndexOfTwoGroup::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->ids1)) this->group1.insert(atom);
                      if (Atom::is_match(atom, this->ids2)) this->group2.insert(atom);
                  });
}


class PrintTopolgy {
public:
    PrintTopolgy() {}

    void action(const std::string &topology_filename);

    static const string title() { return "Print Selected Atoms in Topolgoy File"; }
};

void PrintTopolgy::action(const std::string &topology_filename) {
    TrajectoryReader reader;
    reader.add_topology(topology_filename);
    auto frame = reader.readTopology();

    int sequence = 1;
    for (auto &mol : frame->molecule_list) {
        mol->sequence = sequence++;
    }

    enum class Mode {
        Mass,
        Geom,
        Noop,
    } mode;

    for (;;) {
        beg:
        CenterRuleNode r;
        selectCentergroup(r, "> ");
        Atom::AtomIndenter id;

        if (auto ret = boost::get<std::shared_ptr<MassCenterRuleNode>>(&r)) {
            id.ast = (*ret)->SelectionMask;
            mode = Mode::Mass;
        } else if (auto ret = boost::get<std::shared_ptr<GeomCenterRuleNode>>(&r)) {
            id.ast = (*ret)->SelectionMask;
            mode = Mode::Geom;
        } else if (auto ret = boost::get<std::shared_ptr<NoopRuleNode>>(&r)) {
            id.ast = (*ret)->SelectionMask;
            mode = Mode::Noop;
        }

        std::cout << boost::format("%-6s %-7s %4s %-7s %4s %-6s  %6s %8s  %8s%8s%8s\n")
                     % "#Atom" % "Name" % "#Res" % "Name" % "#Mol" % "Type" % "Charge" % "Mass"
                     % "X(Ang)" % "Y(Ang)" % "Z(Ang)";
        double weight = 0;
        double sum_x = 0.0;
        double sum_y = 0.0;
        double sum_z = 0.0;
        for (auto &atom : frame->atom_list) {
            if (Atom::is_match(atom, id)) {
                std::cout << boost::format("%6d %-7s %4s %-7s %4s %-6s % 6.4f %8s %8.3f%8.3f%8.3f\n")
                             % atom->seq % atom->atom_name
                             % (atom->residue_num ? boost::lexical_cast<std::string>(atom->residue_num.get()) : "-")
                             % (atom->residue_name ? atom->residue_name.get() : "-")
                             % atom->molecule.lock()->sequence.get()
                             % atom->type_name
                             % atom->charge
                             % (atom->mass ? (boost::format("%8.4f") % atom->mass.get()).str() : "-")
                             % atom->x
                             % atom->y
                             % atom->z;
                switch (mode) {
                    case Mode::Mass:
                        if (!atom->mass) {
                            std::cerr << "atom mass not available !\n";
                            goto beg;
                        }
                        sum_x += atom->x * atom->mass.get();
                        sum_y += atom->y * atom->mass.get();
                        sum_z += atom->z * atom->mass.get();
                        weight += atom->mass.get();
                        break;
                    case Mode::Geom:
                        sum_x += atom->x;
                        sum_y += atom->y;
                        sum_z += atom->z;
                        weight++;
                        break;
                }
            }
        }
        if (weight != 0.0) {
            switch (mode) {
                case Mode::Mass:

                    std::cout << boost::format("Mass Center %8s%8s%8s\n") % "X(Ang)" % "Y(Ang)" % "Z(Ang)";
                    std::cout << boost::format("            %8.3f%8.3f%8.3f\n")
                                 % (sum_x / weight)
                                 % (sum_y / weight)
                                 % (sum_z / weight);

                    break;
                case Mode::Geom:
                    std::cout << boost::format("Geom Center %8s%8s%8s\n") % "X(Ang)" % "Y(Ang)" % "Z(Ang)";
                    std::cout << boost::format("            %8.3f%8.3f%8.3f\n")
                                 % (sum_x / weight)
                                 % (sum_y / weight)
                                 % (sum_z / weight);
                    break;
            }
        }

    }
}


void processOneFrame(shared_ptr<Frame> &frame,
                     shared_ptr<list<shared_ptr<BasicAnalysis>>> &task_list) {
    for (auto &task : *task_list) {
        task->process(frame);
    }
}

void processFirstFrame(shared_ptr<Frame> &frame,
                       shared_ptr<list<shared_ptr<BasicAnalysis>>> &task_list) {
    for (auto &task : *task_list) {
        task->processFirstFrame(frame);
    }
}


template<typename T1, typename T2>
struct add_item {
    add_item(T1 &v1, T2 &v2) : v1(v1), v2(v2) {};

    template<typename T>
    void operator()(boost::type<T>) {
        v1.emplace_back(bind(make_shared<T>));
        v2.emplace_back((boost::format("(%d) %s") % v1.size() % T::title()).str());
    }

    T1 &v1;
    T2 &v2;
};

auto getTasks() {
    auto task_list = make_shared<list<shared_ptr<BasicAnalysis>>>();

    using components = mpl::vector<
            GmxTrj,
            Distance,
            CoordinateNumPerFrame,
            RadicalDistribtuionFunction,
            ResidenceTime,
            GreenKubo,
            HBond,
            RMSDCal,
            RMSFCal,
            Cluster,
            NMRRange,
            Diffuse,
            DiffuseCutoff,
            FirstCoordExchangeSearch,
            RotAcfCutoff,
            DipoleAngle,
            DipoleAngle2Gibbs,
            DipoleAngleSingleDistanceNormal,
            DipoleAngleVolumeNormal,
            ShellDensity,
            SearchInteractionResidue,
            FindMinBetweenTwoGroups,
            DemixIndexOfTwoGroup
    >;

    BOOST_MPL_ASSERT((mpl::equal<mpl::unique<components, is_same<mpl::_1, mpl::_2> >::type, components>));

    std::vector<std::function<shared_ptr<BasicAnalysis>()>> task_vec;
    std::vector<string> item_menu;

    mpl::for_each<components, boost::type<mpl::_>>(add_item<
            std::vector<std::function<shared_ptr<BasicAnalysis>()>>,
            std::vector<string>
    >(task_vec, item_menu));

    auto menu1 = [&item_menu]() {
        std::cout << "Please select the desired operation (Trajectory Analysis)" << std::endl;
        std::cout << "(0) Start\n";
        for_each(item_menu.cbegin(), item_menu.cend(), std::cout << phoenix::placeholders::_1 << '\n');
        return choose<int>(0, mpl::size<components>::value, "select :");
    };
    while (true) {
        int num = menu1();
        if (num == 0) return task_list;
        shared_ptr<BasicAnalysis> task = task_vec[num - 1]();

        string line(item_menu[num - 1].size() + 6, '-');

        std::cout << line << "\n";
        std::cout << "<- " << item_menu[num - 1] << " ->\n";
        std::cout << line << "\n";
        task->readInfo();
        task_list->push_back(task);
    }
}

int main(int argc, char *argv[]) {

    std::cout << "Build DateTime : " << __DATE__ << " " << __TIME__ << endl;

    std::cout << "current work dir : " << boost::filesystem::current_path() << std::endl;

    auto desc = make_program_options();
    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    } catch (std::exception &e) {
        std::cerr << e.what() << '\n';
        std::cout << desc;
        exit(EXIT_FAILURE);
    }

    if (vm.count("help")) {
        cout << desc;
        exit(EXIT_SUCCESS);
    }

    std::vector<std::string> xyzfiles;
    if (vm.count("file")) {
        xyzfiles = vm["file"].as<std::vector<string>>();
        for (auto &xyzfile : xyzfiles) {
            if (!boost::filesystem::exists(xyzfile)) {
                cerr << "The file " << xyzfile << " is bad !" << endl;
                exit(EXIT_FAILURE);
            }
        }
    }


    auto mainMenu = []() {
        std::cout << "Main Menu\n";
        std::cout << "(0) Trajectory Analysis\n";
        std::cout << "(1) Print Topology\n";
        return choose<int>(0, 1, "select :");
    };

    if (mainMenu() == 1) {
        PrintTopolgy printer;
        if (vm.count("topology")) {
            string topol = vm["topology"].as<string>();
            if (file_exist(topol)) {
                printer.action(topol);
                return EXIT_SUCCESS;
            }
            cout << "topology file " << topol << " is bad ! please retype !" << endl;
        }
        printer.action(choose_file("input topology file : ", true));
        return EXIT_SUCCESS;
    }

    if (!vm.count("file")) {
        cerr << "input trajectory file is not set !" << endl;
        cerr << desc;
        exit(EXIT_FAILURE);
    }

    auto task_list = getTasks();

    while (enable_forcefield) {
        if (vm.count("prm")) {
            auto ff = vm["prm"].as<string>();
            if (boost::filesystem::exists(ff)) {
                forcefield.read(ff);
                break;
            }
            cout << "force field file " << ff << " is bad ! please retype !" << endl;

        }
        forcefield.read(choose_file("force field filename:", true));
        break;
    }
    int start = choose(1, INT32_MAX, "Enter the start frame[1]:", true, 1);
    int step_size = choose(1, INT32_MAX, "Enter the step size[1]:", true, 1);
    int total_frames = choose(0, INT32_MAX, "How many frames to read [all]:", true);

    tbb::task_scheduler_init tbb_init(tbb::task_scheduler_init::deferred);
    if (enable_tbb) {
        int threads = choose<int>(0, sysconf(_SC_NPROCESSORS_ONLN),
                                  "How many cores to used in parallel[automatic]:", true);
        tbb_init.initialize(threads == 0 ? tbb::task_scheduler_init::automatic : threads);
    }
    int current_frame_num = 0;

    auto reader = make_shared<TrajectoryReader>();
    bool b_added_topology = false;
    for (auto &xyzfile : xyzfiles) {
        reader->add_filename(xyzfile);
        string ext = ext_filename(xyzfile);
        if (ext == "traj" && xyzfiles.size() != 1) {
            cout << "traj file can not use multiple files" << endl;
            exit(EXIT_FAILURE);
        }
        if (!b_added_topology) {
            if (vm.count("topology")) {
                string topol = vm["topology"].as<string>();
                if (file_exist(topol)) {
                    reader->add_topology(topol);
                    b_added_topology = true;
                    continue;
                }
                cout << "topology file " << topol << " is bad ! please retype !" << endl;
            }
            reader->add_topology(choose_file("input topology file : ", true));
            b_added_topology = true;
        }
    }


    auto input_line = input("Do you want to use multiple files [No]:");
    boost::trim(input_line);
    if (!input_line.empty()) {
        if (input_line[0] == 'Y' or input_line[0] == 'y') {
            while (true) {
                input_line = choose_file("next file [Enter for End]:", true, "", true);
                if (input_line.empty()) break;
                if (ext_filename(input_line) == "traj") {
                    cout << "traj file can not use multiple files [retype]" << endl;
                    continue;
                }
                reader->add_filename(input_line);
            }
        }
    }


    if (enable_outfile) {
        outfile.open(vm.count("output") ? vm["output"].as<string>() : choose_file("Output file: ", false),
                     std::fstream::out);
    }
    shared_ptr<Frame> frame;
    int Clear = 0;
    while ((frame = reader->readOneFrame())) {
        current_frame_num++;
        if (total_frames != 0 and current_frame_num > total_frames)
            break;
        if (current_frame_num % 10 == 0) {
            if (Clear) {
                std::cout << "\r";
            }
            cout << "Processing Coordinate Frame  " << current_frame_num << "   " << std::flush;
            Clear = 1;
        }
        if (current_frame_num >= start && (current_frame_num - start) % step_size == 0) {
            if (current_frame_num == start) {
                processFirstFrame(frame, task_list);
            }
            processOneFrame(frame, task_list);
        }
    }
    std::cout << std::endl;


    for (auto &task : *task_list)
        task->print();
    if (outfile.is_open()) outfile.close();
    std::cout << "Mission Complete" << std::endl;

    return EXIT_SUCCESS;

}



