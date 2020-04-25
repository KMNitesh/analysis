//
// Created by xiamr on 6/14/19.
//
#include "GreenKubo.hpp"

#include <tbb/tbb.h>
#include <tbb/tbb_allocator.h>

#include "data_structure/frame.hpp"
#include "utils/common.hpp"

using namespace std;

GreenKubo::GreenKubo() {
    enable_read_velocity = true;
    enable_tbb = true;
    enable_outfile = true;
}

void GreenKubo::print(std::ostream &os) {
    if (steps < 2) {
        cerr << "Too few frame number :" << steps << endl;
        exit(1);
    }
    std::vector<double> vecx(steps);
    std::vector<double> vecy(steps);
    std::vector<double> vecz(steps);
    for (unsigned int i = 0; i < steps; i++) {
        vecx[i] = vecx_map[i + 1];
        vecy[i] = vecy_map[i + 1];
        vecz[i] = vecz_map[i + 1];
    }

    class Body {
    public:
        const std::vector<double> &vecx;
        const std::vector<double> &vecy;
        const std::vector<double> &vecz;

        size_t steps;

        std::vector<double, tbb::tbb_allocator<double>> cxx, cyy, czz;
        std::vector<size_t, tbb::tbb_allocator<size_t>> numbers;

        Body(const std::vector<double> &vecx, const std::vector<double> &vecy, const std::vector<double> &vecz,
             size_t steps)
            : vecx(vecx), vecy(vecy), vecz(vecz), steps(steps), cxx(steps), cyy(steps), czz(steps), numbers(steps) {}

        Body(const Body &c, tbb::split)
            : vecx(c.vecx), vecy(c.vecy), vecz(c.vecz), steps(c.steps), cxx(c.steps), cyy(c.steps), czz(c.steps),
              numbers(c.steps) {}

        void join(const Body &y) {
            for (size_t step = 0; step < steps - 1; step++) {
                cxx[step] += y.cxx[step];
                cyy[step] += y.cyy[step];
                czz[step] += y.czz[step];
                numbers[step] += y.numbers[step];
            }
        }

        void operator()(const tbb::blocked_range<size_t> &r) {
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
    os << "Green-Kubo self-diffuse " << endl;
    os << "selected group : " << ids << endl;
    os << "Time (ps)             D(10^-5 cm**2/sec)" << endl;
    double pre_c = body.cxx[0] + body.cyy[0] + body.czz[0];
    double integral = 0.0;
    for (unsigned int i = 1; i < steps; i++) {
        double c = body.cxx[i] + body.cyy[i] + body.czz[i];
        integral += 5 * (pre_c + c) * timestep / 3;
        pre_c = c;
        os << i * timestep << "    " << integral << endl;
    }
}

void GreenKubo::readInfo() {
    select1group(ids, "Please enter atom group");
    this->timestep =
        choose(0.0, static_cast<double>(numeric_limits<int>::max()), "Please enter time step for each frame(ps):");
}

void GreenKubo::process(std::shared_ptr<Frame> &frame) {
    if (!frame->has_velocity) {
        std::cerr << "Trajectory does not have velocity\n";
        std::exit(EXIT_FAILURE);
    }

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
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(), [this](shared_ptr<Atom> &atom) {
        if (is_match(atom, this->ids))
            this->group.insert(atom);
    });
}
