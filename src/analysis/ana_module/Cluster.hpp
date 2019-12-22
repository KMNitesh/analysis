//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_CLUSTER_HPP
#define TINKER_CLUSTER_HPP


#include "utils/std.hpp"
#include "AbstractAnalysis.hpp"
#include "data_structure/atom.hpp"
#include "RMSDCal.hpp"

class Frame;

class Cluster : public AbstractAnalysis {
public:

    Cluster();

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void readInfo() override;

    [[nodiscard]] static std::string_view title() { return "Cluster Analysis(linkage)"; }

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


    template<typename ForwardRange>
    static std::vector<Cluster::conf_clust> do_cluster(const ForwardRange &rmsd_list, int conf_size, double cutoff);

    static int do_sort_and_renumber_parallel(std::vector<conf_clust> &conf_clust_vector);

protected:

    AmberMask ids;
    std::vector<std::shared_ptr<Atom>> group;


    double cutoff = 0.0;
    int steps = 0; // current frame number
    std::map<std::pair<int, int>, double> rmsd_map;

    std::vector<std::vector<std::tuple<double, double, double>>> coords;

    double rmsvalue(int index1, int index2);

    static void do_sort_clust_parallel(std::vector<conf_clust> &conf_clust_vector);

    static void do_sort_conf_parallel(std::vector<conf_clust> &conf_clust_vector);

    static int do_renumber_clust(std::vector<conf_clust> &conf_clust_vector);

    static std::vector<Cluster::conf_clust> initialize_conf_clust_vector(int conf_size);

    std::list<Cluster::rmsd_matrix> do_calculate_rmsd_list_parallel();

    void setSetting(const Atom::AmberMask &atomIndenter, double cutoff);
};

std::unordered_map<int, std::vector<int>> do_find_frames_in_same_clust(const std::vector<Cluster::conf_clust> &clusts);

std::unordered_map<int, std::pair<int, double>> do_find_medium_in_clust(
        const std::vector<Cluster::conf_clust> &clusts, const std::list<Cluster::rmsd_matrix> &rmsd_list);

#endif //TINKER_CLUSTER_HPP
