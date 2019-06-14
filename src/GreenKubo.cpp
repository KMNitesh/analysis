//
// Created by xiamr on 6/14/19.
//
#include <tbb/tbb.h>

#include "GreenKubo.hpp"
#include "frame.hpp"


using namespace std;

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

void GreenKubo::process(std::shared_ptr<Frame> &) {


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