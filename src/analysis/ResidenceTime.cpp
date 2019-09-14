//
// Created by xiamr on 6/14/19.
//

#include <tbb/tbb.h>
#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm.hpp>
#include "ResidenceTime.hpp"
#include "frame.hpp"

using namespace std;

void ResidenceTime::calculate() {

    std::vector<std::vector<int>> atom_star_map(atom_num);
    for (int atom = 0; atom < atom_num; atom++) {
        int count1 = 0;
        bool swi = true;
        for (int k = 0; k < steps; k++) {
            if ((!mark(k, atom)) && swi) {
                swi = false;
                count1 = k;
            } else if (mark(k, atom) && (!swi)) {
                swi = true;
                if (timeStarMode == TimeStarMode::Loose) {
                    if (k - 1 - count1 > time_star) atom_star_map[atom].emplace_back(count1);
                } else {
                    if (k - count1 > time_star) atom_star_map[atom].emplace_back(count1);
                }
            }
        }
    }

    std::vector<std::vector<std::pair<int, int>>> hydrationed_atoms(steps);

    for (int step = 0; step < steps; ++step) {
        for (int atom = 0; atom < atom_num; atom++) {
            if (mark(step, atom)) {
                int out_frame = steps + 1;
                for (auto out_frame_index : atom_star_map[atom]) {
                    if (step < out_frame_index) {
                        out_frame = out_frame_index;
                        break;
                    }
                }
                hydrationed_atoms[step].emplace_back(atom, out_frame);
            }
        }
    }
    atom_star_map.clear();

    class Body {
    public:
        int steps;

        Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> &mark;
        std::vector<int, tbb::tbb_allocator<int>> time_array;
        std::vector<double, tbb::tbb_allocator<double>> Rt_array;
        std::vector<std::vector<std::pair<int, int>>> &hydrationed_atoms;

        Body(int steps, Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> &mark,
             std::vector<std::vector<std::pair<int, int>>> &hydrationed_atoms) :
                steps(steps), mark(mark),
                time_array(steps - 1), Rt_array(steps - 1), hydrationed_atoms(hydrationed_atoms) {}

        Body(Body &body, tbb::split) :
                steps(body.steps), mark(body.mark),
                time_array(body.steps - 1), Rt_array(body.steps - 1), hydrationed_atoms(body.hydrationed_atoms) {}

        void join(const Body &rhs) {
            for (int step = 0; step < steps - 1; step++) {
                time_array[step] += rhs.time_array[step];
                Rt_array[step] += rhs.Rt_array[step];
            }
        }

        void operator()(const tbb::blocked_range<int> &r) {
            for (auto i = r.begin(); i != r.end(); ++i) {
                auto CN = hydrationed_atoms[i].size();
                if (CN == 0) {
                    std::cerr << "Warning !!! Coordination Number is zero,  skip frame (start from 0) = " << i << "\n";
                    continue;
                }
                auto factor = 1.0 / CN;
                for (int j = i + 1; j < steps; j++) {
                    int value = 0.0;
                    auto n = j - i - 1;
                    for (auto[atom, outframe] : hydrationed_atoms[i]) {
                        if (mark(j, atom) and j < outframe) value++;
                    }
                    ++time_array[n];
                    Rt_array[n] += value * factor;
                }
            }
        }
    } body(steps, mark, hydrationed_atoms);
    tbb::parallel_reduce(tbb::blocked_range<int>(0, steps - 1), body);
    for (int i = 0; i < steps - 1; i++)
        body.Rt_array[i] /= body.time_array[i];

    Rt_array = std::move(body.Rt_array);
}


void ResidenceTime::process(std::shared_ptr<Frame> &frame) {

    steps++;
    int atom_no = 0;
    for (auto &atom1 : center_atom_group) {
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

void ResidenceTime::print(std::ostream &os) {
    if (steps < 2) {
        cerr << "Too few frame number :" << steps << endl;
        exit(1);
    }
    atom_num = static_cast<int>(mark_map.size());
    mark = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>::Zero(steps, atom_num);
    for (const auto &it : mark_map) {
        auto atom_no = it.first;
        for (const auto &ele : it.second | boost::adaptors::indexed(0)) {
            mark(ele.index(), atom_no) = ele.value();
        }
    }
    calculate();

    os << "# " << title() << '\n';
    os << "# t* mode : " << (timeStarMode == TimeStarMode::Loose ? "Loose" : "Strict") << '\n';
    os << "# center_atom_mask: " << center_atom_mask << ",  Ow_atom_mask: " << Ow_atom_mask << '\n';
    os << "# dist_cutoff (Ang) = " << dis_cutoff << '\n';
    os << "# t* (frame) = " << time_star << '\n';

    os << format("#%15s %15s\n", "Frame", "R(t)");
    for (const auto &ele : Rt_array | boost::adaptors::indexed(1)) {
        os << format(" %15d %15.6f\n", ele.index(), ele.value());
    }
}

void ResidenceTime::readInfo() {

    Atom::select2group(center_atom_mask, Ow_atom_mask, "Enter mask for center atom > ", "Enter mask for Ow atoms > ");

    dis_cutoff = choose(0.0, "Please enter distance cutoff1: ");
    time_star = choose(0, "Please enter t*: ( unit: frame) : ");

    timeStarMode = static_cast<TimeStarMode >(choose(0, 1, "t* mode\n(0) Loose\n(1) Strict\nchoose : "));
}

void ResidenceTime::setParameters(const Atom::Node &id1, const Atom::Node &id2,
                                  double cutoff, int t_star, const std::string &outfilename) {
    this->center_atom_mask = id1;
    this->Ow_atom_mask = id2;


    if (cutoff <= 0) {
        throw runtime_error("cutoff must large than zero");
    }

    this->dis_cutoff = cutoff;

    if (t_star == 0) {
        throw runtime_error("t_star(t*) must large than zero");
    }
    this->time_star = t_star;

    this->outfilename = outfilename;
    boost::trim(this->outfilename);
    if (this->outfilename.empty()) {
        throw runtime_error("outfilename cannot empty");
    }

}

void ResidenceTime::processFirstFrame(std::shared_ptr<Frame> &frame) {
    boost::for_each(frame->atom_list,
                    [this](shared_ptr<Atom> &atom) {
                        if (Atom::is_match(atom, this->center_atom_mask)) this->center_atom_group.insert(atom);
                        if (Atom::is_match(atom, this->Ow_atom_mask)) this->group2.insert(atom);
                    });
}

ResidenceTime::~ResidenceTime() = default;

