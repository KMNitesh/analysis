//
// Created by xiamr on 9/20/19.
//

#ifndef TINKER_COORDINATIONSTRUCTURECLASSIFICATION_HPP
#define TINKER_COORDINATIONSTRUCTURECLASSIFICATION_HPP

#include "std.hpp"
#include "AbstractAnalysis.hpp"
#include "atom.hpp"
#include "Cluster.hpp"

class Frame;

class CoordinationStructureClassification : public AbstractAnalysis {
public:

    CoordinationStructureClassification();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    [[nodiscard]] static std::string_view title() { return "Coordination Structure Classification based on RMSD"; }

    [[nodiscard]] static double calculateRmsdOfTwoStructs(std::vector<std::tuple<double, double, double>> &c1,
                                                          std::vector<std::tuple<double, double, double>> &c2);

protected:

    static std::pair<int, int> find_min_distance_pair(std::vector<std::tuple<double, double, double>> &c);

    static void fill_coord(double x[], double y[], double z[], int index1, int index2,
                           std::vector<std::tuple<double, double, double>> &c);

    static void permutation(double x1[], double y1[], double z1[],
                            double x2[], double y2[], double z2[], int start, int end/*not included*/);

    [[nodiscard]] std::map<int, std::list<Cluster::rmsd_matrix>> do_calculate_rmsd_list_parallel();

    AmberMask metal_mask;
    AmberMask Ow_atom_mask;

    std::shared_ptr<Atom> metal;
    std::vector<std::shared_ptr<Atom>> Ow_atoms;
    double cutoff2;
    double rmsd_cutoff;
    std::map<int, std::deque<std::pair<int, std::vector<std::tuple<double, double, double>>>>> systems;
    int nframe = 0;

    std::unordered_map<int, int> frame_cn_mapping;

    bool output_rms_matrix = false;

};


#endif //TINKER_COORDINATIONSTRUCTURECLASSIFICATION_HPP
