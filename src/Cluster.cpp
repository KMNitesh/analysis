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


void Cluster::print() {

    class CalCore {
    public:
        list<rmsd_matrix> local_rms_list;
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
                                 [](const rmsd_matrix &m1, const rmsd_matrix &m2) { return (m1.rms < m2.rms); });
        }

        void operator()(const tbb::blocked_range<int> &range) {
            double rms;
            for (int index1 = range.begin(); index1 != range.end(); ++index1) {
                for (int index2 = index1 + 1; index2 < steps; ++index2) {
                    rms = parent->rmsvalue(index1, index2);
                    if (rms < cutoff) local_rms_list.emplace_back(index1, index2, rms);
                }
            }
            local_rms_list.sort([](const rmsd_matrix &m1, const rmsd_matrix &m2) { return (m1.rms < m2.rms); });
        }
    } core(steps, cutoff, this);

    tbb::parallel_reduce(tbb::blocked_range<int>(0, steps - 1), core, tbb::auto_partitioner());


    auto &rmsd_list = core.local_rms_list;

    int n = 0;
    for (auto &v : rmsd_list) {
        cout << v.i << "  " << v.j << "  " << v.rms << endl;
        n++;
        if (n > 99) break;
    }


    vector<conf_clust> c = do_cluster(rmsd_list, steps);


    cout << "Sorting and renumbering clusters..." << endl;
    int cid = do_sort_and_renumber_parallel(c);


    outfile << "***************************" << endl;
    outfile << "*Cluster Analysis(Linkage)*" << endl;
    outfile << "SET:" << ids << endl;
    outfile << "cutoff : " << this->cutoff << endl;
    outfile << "Total cluster number : " << cid << '\n';
    outfile << "***************************" << endl;
    for (unsigned int k = 0; k < c.size(); k++) {
        outfile << k + 1 << "   " << c[k].clust << endl;
    }
    outfile << "***************************" << endl;

}

/*
 * Sort and Renumber the cluster vector, make cluster ID start with 1 and return the total cluster amount
 *
 */

int Cluster::do_sort_and_renumber_parallel(vector<conf_clust> &conf_clust_vector) const {

    do_sort_clust_parallel(conf_clust_vector);
    int total_clust = do_renumber_clust(conf_clust_vector);
    do_sort_conf_parallel(conf_clust_vector);
    return total_clust;
}

int Cluster::do_renumber_clust(vector<Cluster::conf_clust> &conf_clust_vector) const {
    int cid = 1;
    unsigned int k;
    for (k = 1; k < conf_clust_vector.size(); k++) {
        if (conf_clust_vector[k].clust != conf_clust_vector[k - 1].clust) {
            conf_clust_vector[k - 1].clust = cid;
            cid++;
        } else {
            conf_clust_vector[k - 1].clust = cid;
        }
    }
    conf_clust_vector[k - 1].clust = cid;

    return cid;
}

void Cluster::do_sort_conf_parallel(vector<Cluster::conf_clust> &conf_clust_vector) const {
    tbb::parallel_sort(conf_clust_vector.begin(), conf_clust_vector.end(),
                       [](const conf_clust &i, const conf_clust &j) { return (i.conf < j.conf); });
}

void Cluster::do_sort_clust_parallel(vector<Cluster::conf_clust> &conf_clust_vector) const {
    tbb::parallel_sort(conf_clust_vector.begin(), conf_clust_vector.end(),
                       [](const conf_clust &i, const conf_clust &j) { return (i.clust < j.clust); });
}

/*
 *  give the rmsd matrix between frames, return cluster number assgined to frames
 *  the cluster numbers are not have order or continuous
 *
 */

vector<Cluster::conf_clust> Cluster::do_cluster(const list<rmsd_matrix> &rmsd_list, int conf_size) const {
    vector<Cluster::conf_clust> c = initialize_conf_clust_vector(conf_size);

    // The algorithm of blow code block comes from gromacs
    bool bChange;
    do {
        bChange = false;
        for (auto &k : rmsd_list) {
            //        for(int n = 0; n < num; n++){
            //            rmsd_matrix k = *(vect[n]);
            if (k.rms >= cutoff) break;
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
    /////////////////////////////////////////


    return c;
}

vector<Cluster::conf_clust> Cluster::initialize_conf_clust_vector(int conf_size) const {
    vector<conf_clust> c;
    for (int i = 0; i < conf_size; i++) {
        c.emplace_back(i, i);
    }
    return c;
}

void Cluster::readInfo() {

    Atom::select1group(ids, "Please enter group > ");
    this->cutoff = choose(0.0, static_cast<double>(numeric_limits<int>::max()), "Cutoff > ");
}

void Cluster::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->ids)) this->group.insert(atom);
                  });
}