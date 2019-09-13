//
// Created by xiamr on 6/14/19.
//

#include <tbb/tbb.h>
#include <boost/range/adaptors.hpp>
#include "ResidenceTime.hpp"
#include "frame.hpp"

using namespace std;

void ResidenceTime::setSteps(size_t steps, int atom_num) {
    this->steps = steps;
    delete[] time_array;
    time_array = new int[steps - 1];
    bzero(time_array, sizeof(int) * (steps - 1));

    delete[] Rt_array;
    Rt_array = new double[steps - 1];
    bzero(Rt_array, sizeof(double) * (steps - 1));
    this->atom_num = atom_num;
    mark = Eigen::MatrixXi::Zero(steps, atom_num);
}

void ResidenceTime::calculate() {

    std::vector<std::vector<std::pair<unsigned int, unsigned int>>> atom_star_map(atom_num);
    for (int atom = 0; atom < atom_num; atom++) {
        unsigned int count1 = 0;
        bool swi = true;
        for (unsigned int k = 0; k < steps; k++) {
            if ((!mark(k, atom)) && swi) {
                swi = false;
                count1 = k;
            } else if (mark(k, atom) && (!swi)) {
                swi = true;
                if (k - 1 - count1 > time_star) atom_star_map[atom].emplace_back(count1, k - 1);
            }
        }
    }

    std::vector<std::vector<std::pair<int, int>>> hydrationed_atoms(steps);

    for (std::size_t step = 0; step < steps; ++step) {
        for (int atom = 0; atom < atom_num; atom++) {
            if (mark(step, atom)) {
                int out_frame = steps;
                for (auto[a, b] : atom_star_map[atom]) {
                    if (step < a) {
                        out_frame = b;
                        break;
                    }
                }
                hydrationed_atoms[step].emplace_back(atom, out_frame);
            }
        }
    }

    class Body {
    public:
        size_t steps;
        int atom_num;

        std::vector<int, tbb::tbb_allocator<int>> time_array;
        std::vector<double, tbb::tbb_allocator<double>> Rt_array;
        std::vector<std::vector<std::pair<int, int>>> &hydrationed_atoms;

        Body(std::size_t steps, int atom_num,
             std::vector<std::vector<std::pair<int, int>>> &hydrationed_atoms) :
                steps(steps), atom_num(atom_num),
                time_array(steps - 1), Rt_array(steps - 1), hydrationed_atoms(hydrationed_atoms) {}

        Body(Body &body, tbb::split) :
                steps(body.steps), atom_num(body.atom_num),
                time_array(body.steps - 1), Rt_array(body.steps - 1), hydrationed_atoms(body.hydrationed_atoms) {}

        void join(const Body &rhs) {
            for (std::size_t step = 0; step < steps - 1; step++) {
                time_array[step] += rhs.time_array[step];
                Rt_array[step] += rhs.Rt_array[step];
            }
        }

        void operator()(const tbb::blocked_range<std::size_t> &r) {
            for (size_t i = r.begin(); i != r.end(); ++i) {
                auto CN = hydrationed_atoms[i].size();
                if (CN == 0) {
                    std::cerr << "Warning !!! Coordination Number is zero,  skip frame (start from 0) = " << i << "\n";
                    continue;
                }
                auto factor = 1.0 / CN;
                for (size_t j = i + 1; j < steps; j++) {
                    int value = 0.0;
                    auto n = j - i - 1;
                    for (auto[atom, outframe] : hydrationed_atoms[i]) {
                        if (j < outframe) value++;
                    }
                    ++time_array[n];
                    Rt_array[n] += value * factor;
                }
            }

        }
    } body(steps, atom_num, hydrationed_atoms);
    tbb::parallel_reduce(tbb::blocked_range<size_t>(0, steps - 1), body);
    for (std::size_t step = 0; step < steps - 1; step++) {
        Rt_array[step] = body.Rt_array[step];
        time_array[step] = body.time_array[step];
    }

    for (std::size_t i = 0; i < steps - 1; i++)
        Rt_array[i] /= time_array[i];
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

void ResidenceTime::print(std::ostream &os) {
    if (steps < 2) {
        cerr << "Too few frame number :" << steps << endl;
        exit(1);
    }
    setSteps(steps, static_cast<int>(mark_map.size()));
    for (const auto &it : mark_map) {
        auto atom_no = it.first;
        for (const auto &ele : it.second | boost::adaptors::indexed(0)) {
            mark(ele.index(), atom_no) = int(ele.value());
        }
    }
    calculate();

    os << "# typ1: " << ids1 << ",  typ2: " << ids2 << "  dist_cutoff = " << dis_cutoff << " t* = "
       << time_star
       << std::endl;

    os << "# Frame        R " << std::endl;
    for (unsigned int i = 0; i < steps - 1; i++)
        os << i + 1 << "    " << Rt_array[i] << std::endl;
}

void ResidenceTime::readInfo() {

    Atom::select2group(ids1, ids2);

    std::cout << "Please enter distance cutoff1:";
    std::cin >> dis_cutoff;
    std::cout << "Please enter t*: ( unit: frame)";
    std::cin >> time_star;

}

void ResidenceTime::setParameters(const Atom::Node &id1, const Atom::Node &id2,
                                  double cutoff, int t_star, const std::string &outfilename) {
    this->ids1 = id1;
    this->ids2 = id2;


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
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->ids1)) this->group1.insert(atom);
                      if (Atom::is_match(atom, this->ids2)) this->group2.insert(atom);
                  });
}

ResidenceTime::~ResidenceTime() {
    boost::checked_array_delete(time_array);
    boost::checked_array_delete(Rt_array);
}

