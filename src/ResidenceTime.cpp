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
                                    cerr << boost::format(" i = %d, j = %d, a = %d, d = %d\n") % i % j % a % b;
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

void ResidenceTime::print(std::ostream &os) {
    if (steps < 2) {
        cerr << "Too few frame number :" << steps << endl;
        exit(1);
    }
    setSteps(steps, static_cast<int>(mark_map.size()));
    for (const auto &it : mark_map) {
        auto atom_no = it.first;
        for (const auto &ele : it.second | boost::adaptors::indexed(0)) {
            mark[atom_no][ele.index()] = int(ele.value());
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
    if (mark) {
        for (int i = 0; i < atom_num; ++i)
            boost::checked_array_delete(mark[i]);
        boost::checked_array_delete(mark);
    }
}

