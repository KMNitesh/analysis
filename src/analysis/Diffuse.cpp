//
// Created by xiamr on 6/14/19.
//
#include <tbb/tbb.h>

#include "Diffuse.hpp"
#include "frame.hpp"
#include "molecule.hpp"
#include <boost/range/algorithm.hpp>

using namespace std;

Diffuse::Diffuse() {
    enable_forcefield = true;
    enable_outfile = true;
    enable_tbb = true;
}

namespace {
    class wrap {
        const std::tuple<double, double, double> &rhs;
    public:
        explicit wrap(const std::tuple<double, double, double> &rhs) : rhs(rhs) {}

        operator Eigen::Array3d() { return {std::get<0>(rhs), std::get<1>(rhs), std::get<2>(rhs)}; }
    };
}

void Diffuse::process(std::shared_ptr<Frame> &frame) {
    int mol_index = 0;
    if (steps == 0) {
        xyzcm.resize(total_frame_number, mols.size());
        for (auto &mol : mols) {
            xyzcm(steps, mol_index) = wrap(mol->calc_weigh_center(frame));
            mol_index++;
        }
    } else {
        for (auto &mol : mols) {
            Eigen::Array3d coord = wrap(mol->calc_weigh_center(frame));
            auto &xyzold = xyzcm(steps - 1, mol_index);
            Eigen::Array3d r = coord - xyzold;
            frame->image(r);
            xyzcm(steps, mol_index) = r + xyzold;
            mol_index++;
        }
    }
    steps++;
}

void Diffuse::print(std::ostream &os) {

    class Body {
    public:
        int total_frame_number;
        std::vector<Eigen::Array3d, tbb::tbb_allocator<Eigen::Array3d>> msd;

        int total_mol;
        Eigen::Matrix<Eigen::Array3d, Eigen::Dynamic, Eigen::Dynamic> &xyzcm;

        Body(int total_frame_number, int total_mol,
             Eigen::Matrix<Eigen::Array3d, Eigen::Dynamic, Eigen::Dynamic> &xyzcm) :
                total_frame_number(total_frame_number),
                msd(total_frame_number - 1),
                total_mol(total_mol), xyzcm(xyzcm) {}

        Body(const Body &body, tbb::split) :
                total_frame_number(body.total_frame_number),
                msd(body.total_frame_number - 1),
                total_mol(body.total_mol), xyzcm(body.xyzcm) {}

        void join(const Body &body) {
            for (int i = 0; i < total_frame_number - 1; i++) {
                msd[i] += body.msd[i];
            }
        }

        void operator()(const tbb::blocked_range<int> &range) {
            for (int i = range.begin(); i != range.end(); i++) {
                for (int j = i + 1; j < total_frame_number; j++) {
                    int m = j - i - 1;
                    for (int k = 0; k < total_mol; k++) {
                        Eigen::Array3d diff = xyzcm(j, k) - xyzcm(i, k);
                        msd[m] += diff.square();
                    }
                }
            }
        }
    } body(total_frame_number, mols.size(), xyzcm);

    tbb::parallel_reduce(tbb::blocked_range<int>(0, total_frame_number - 1), body, tbb::auto_partitioner());


    constexpr double dunits = 10.0;

    for (int i = 0; i < total_frame_number - 1; i++) {
        auto counts = mols.size() * (total_frame_number - (i + 1));
        body.msd[i] /= counts;
    }

    os << "*********************************************************\n";
    os << description();
    os << "Mean Squared Displacements and Self-Diffusion Constant\n";
    os << "    Time Gap      X MSD       Y MSD       Z MSD       R MSD        Diff Const\n";
    os << "      (ps)       (Ang^2)     (Ang^2)     (Ang^2)     (Ang^2)     (x 10^-5 cm**2/sec)\n";

    for (int i = 0; i < total_frame_number - 1; i++) {
        double delta = time_increment_ps * (i + 1);
        double xvalue = body.msd[i][0];
        double yvalue = body.msd[i][1];
        double zvalue = body.msd[i][2];
        double rvalue = xvalue + yvalue + zvalue;
        double dvalue = dunits * rvalue / delta / 6.0;
        os << boost::format("%12.2f%12.2f%12.2f%12.2f%12.2f%12.4f\n") %
              delta % xvalue % yvalue % zvalue % rvalue % dvalue;
    }
    os << "*********************************************************\n";

}

void Diffuse::setParameters(const Atom::Node &mask, double time_increment_ps, int total_frames,
                            const std::string &outfilename) {

    this->ids = mask;

    if (time_increment_ps <= 0) {
        throw runtime_error("`time_increment_ps' must great than zero");
    }
    this->time_increment_ps = time_increment_ps;


    if (total_frames <= 0) {
        throw runtime_error("`total_frames' must great than zero");
    }
    this->total_frame_number = total_frames;

    this->outfilename = outfilename;
    boost::trim(this->outfilename);
    if (this->outfilename.empty()) {
        throw runtime_error("outfilename cannot empty");
    }
}

void Diffuse::readInfo() {
    time_increment_ps = choose(0.0, "Enter the Time Increment in Picoseconds [0.1]: ", Default(0.1));
    total_frame_number = choose(1, "Enter the Total Frame Number: ");
    Atom::select1group(ids, "select group :");
}

void Diffuse::processFirstFrame(std::shared_ptr<Frame> &frame) {
    boost::for_each(frame->atom_list,
                    [this](shared_ptr<Atom> &atom) {
                        if (Atom::is_match(atom, ids)) mols.insert(atom->molecule.lock());
                    });
}

string Diffuse::description() {
    stringstream ss;
    string title_line = "------ " + title() + " ------";
    ss << title_line << "\n";
    ss << " mask              = [ " << ids << " ]\n";
    ss << " time_increment_ps = " << time_increment_ps << " (ps)\n";
    ss << " total_frames      = " << total_frame_number << " (frames)\n";
    ss << " outfilename       = " << outfilename << "\n";
    ss << string(title_line.size(), '-') << '\n';
    return ss.str();
}

