//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_CLUSTER_HPP
#define TINKER_CLUSTER_HPP

#include <memory>
#include <unordered_set>
#include <unordered_map>
#include <string>
#include <map>
#include <vector>
#include <list>
#include <utility>

#include "common.hpp"
#include "BasicAnalysis.hpp"
#include "atom.hpp"
#include "RMSDCal.hpp"

class Frame;

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

    static const std::string title() {
        return "Cluster Analysis(linkage)";
    }

private:

    Atom::AtomIndenter ids;
    std::unordered_set<std::shared_ptr<Atom>> group;


    double cutoff = 0.0;
    int steps = 0; // current frame number
    std::map<std::pair<int, int>, double> rmsd_map;

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

    std::map<int, std::map<int, double>> x, y, z;

    double rmsvalue(int index1, int index2);

public:

    class rmsd_matrix {
    public:
        rmsd_matrix(int i, int j, double rms) : i(i), j(j), rms(rms) {}

        int i;
        int j;
        double rms;
    };

    class conf_clust {
    public:
        conf_clust(int conf, int clust) : conf(conf), clust(clust) {};
        int conf;
        int clust;
    };

protected:

    std::vector<Cluster::conf_clust> do_cluster(const std::list<rmsd_matrix> &rmsd_list, int conf_size) const;

    int do_sort_and_renumber_parallel(std::vector<conf_clust> &conf_clust_vector) const;

    void do_sort_clust_parallel(std::vector<conf_clust> &conf_clust_vector) const;

    void do_sort_conf_parallel(std::vector<conf_clust> &conf_clust_vector) const;

    int do_renumber_clust(std::vector<conf_clust> &conf_clust_vector) const;

    std::vector<Cluster::conf_clust> initialize_conf_clust_vector(int conf_size) const;

    std::list<Cluster::rmsd_matrix> do_calculate_rmsd_list_parallel();

    void setSetting(const Atom::AtomIndenter &atomIndenter, double cutoff);
};

std::unordered_map<int, std::vector<int>> do_find_frames_in_same_clust(const std::vector<Cluster::conf_clust> &clusts);

std::unordered_map<int, std::pair<int, double>> do_find_medium_in_clust(
        const std::vector<Cluster::conf_clust> &clusts, const std::list<Cluster::rmsd_matrix> &rmsd_list);

#endif //TINKER_CLUSTER_HPP
