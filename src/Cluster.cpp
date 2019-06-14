//
// Created by xiamr on 6/14/19.
//
#include <tbb/tbb.h>

#include "Cluster.hpp"
#include "frame.hpp"

using namespace std;

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