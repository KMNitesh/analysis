//
// Created by xiamr on 6/14/19.
//
#include "Diffuse.hpp"

#include <tbb/tbb.h>

#include <boost/range/algorithm.hpp>

#include "data_structure/frame.hpp"
#include "data_structure/molecule.hpp"

Diffuse::Diffuse() {
    enable_forcefield = true;
    enable_outfile = true;
    enable_tbb = true;
}

void Diffuse::process(std::shared_ptr<Frame> &frame) {
    int mol_index = 0;
    if (steps == 0) {
        xyzcm.resize(total_frame_number, mols.size());
        for (auto &mol : mols) {
            xyzcm(steps, mol_index) = mol->calc_weigh_center(frame);
            mol_index++;
        }
    } else {
        for (auto &mol : mols) {
            auto coord = mol->calc_weigh_center(frame);
            auto &xyzold = xyzcm(steps - 1, mol_index);
            auto r = coord - xyzold;
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
        std::vector<std::tuple<double, double, double>, tbb::tbb_allocator<std::tuple<double, double, double>>> msd;

        int total_mol;
        Eigen::Matrix<std::tuple<double, double, double>, Eigen::Dynamic, Eigen::Dynamic> &xyzcm;

        Body(int total_frame_number, int total_mol,
             Eigen::Matrix<std::tuple<double, double, double>, Eigen::Dynamic, Eigen::Dynamic> &xyzcm)
            : total_frame_number(total_frame_number), msd(total_frame_number - 1), total_mol(total_mol), xyzcm(xyzcm) {}

        Body(const Body &body, tbb::split)
            : total_frame_number(body.total_frame_number),
              msd(body.total_frame_number - 1),
              total_mol(body.total_mol),
              xyzcm(body.xyzcm) {}

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
                        auto [xr, yr, zr] = xyzcm(j, k) - xyzcm(i, k);
                        msd[m] += std::make_tuple(xr * xr, yr * yr, zr * zr);
                    }
                }
            }
        }
    } body(total_frame_number, mols.size(), xyzcm);

    tbb::parallel_reduce(tbb::blocked_range<int>(0, total_frame_number - 1), body);

    constexpr double dunits = 10.0;

    for (int i = 0; i < total_frame_number - 1; i++) {
        auto counts = mols.size() * (total_frame_number - (i + 1));
        body.msd[i] /= counts;
    }

    os << description();
    os << "    Time Gap      X MSD       Y MSD       Z MSD       R MSD      Diff Const\n";
    os << "      (ps)       (Ang^2)     (Ang^2)     (Ang^2)     (Ang^2)   (x 10^-5 cm**2/sec)\n";

    for (int i = 0; i < total_frame_number - 1; i++) {
        double delta = time_increment_ps * (i + 1);
        auto [xvalue, yvalue, zvalue] = body.msd[i];
        double rvalue = xvalue + yvalue + zvalue;
        double dvalue = dunits * rvalue / delta / 6.0;
        os << boost::format("%12.2f%12.2f%12.2f%12.2f%12.2f%14.4f\n") % delta % xvalue % yvalue % zvalue % rvalue %
                  dvalue;
    }
    os << std::string(50, '#') << '\n';
}

void Diffuse::setParameters(const AmberMask &mask, double time_increment_ps, int total_frames,
                            const std::string &outfilename) {
    this->mask = mask;

    if (time_increment_ps <= 0) {
        throw std::runtime_error("`time_increment_ps' must great than zero");
    }
    this->time_increment_ps = time_increment_ps;

    if (total_frames <= 0) {
        throw std::runtime_error("`total_frames' must great than zero");
    }
    this->total_frame_number = total_frames;

    this->outfilename = outfilename;
    boost::trim(this->outfilename);
    if (this->outfilename.empty()) {
        throw std::runtime_error("outfilename cannot empty");
    }
}

void Diffuse::readInfo() {
    time_increment_ps = choose(0.0, "Enter the Time Increment in Picoseconds [0.1]: ", Default(0.1));
    total_frame_number = choose(1, "Enter the Total Frame Number: ");
    select1group(mask, "select group :");
}

void Diffuse::processFirstFrame(std::shared_ptr<Frame> &frame) {
    boost::for_each(frame->atom_list, [this](std::shared_ptr<Atom> &atom) {
        if (is_match(atom, mask)) mols.insert(atom->molecule.lock());
    });
}

std::string Diffuse::description() {
    std::stringstream ss;
    std::string title_line = "------ " + std::string(title()) + " ------";
    ss << title_line << "\n";
    ss << " mask              = [ " << mask << " ]\n";
    ss << " time_increment_ps = " << time_increment_ps << " (ps)\n";
    ss << " total_frames      = " << total_frame_number << " (frames)\n";
    if (!outfilename.empty())    ss << " outfilename       = " << outfilename << "\n";
    ss << std::string(title_line.size(), '-') << '\n';
    return ss.str();
}
