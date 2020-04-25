//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_RMSDCAL_HPP
#define TINKER_RMSDCAL_HPP

#include <boost/graph/breadth_first_search.hpp>
#include <map>
#include <unordered_set>

#include "AbstractAnalysis.hpp"
#include "data_structure/atom.hpp"
#include "dsl/AmberMask.hpp"
#include "utils/PBCUtils.hpp"
#include "utils/common.hpp"
#include "utils/xtc_writer.hpp"

class RMSDCal : public AbstractAnalysis {
public:
    RMSDCal();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    [[nodiscard]] static std::string_view title() { return "RMSD Calculator"; }

    static double rmsfit(double x1[], double y1[], double z1[], double x2[], double y2[], double z2[], int n_rms_calc);

    static double rms_max(double x1[], double y1[], double z1[], double x2[], double y2[], double z2[], int n_rms_calc);

    static void jacobi(int n, double a[4][4], double d[], double v[4][4]);

    static void quatfit(int n1, double x1[], double y1[], double z1[], int n_rms_calc, double x2[], double y2[],
                        double z2[], int nfit);

    static void center(int n1, double x1[], double y1[], double z1[], double mid[], int nfit);

private:
    std::deque<double> rmsds;
    bool first_frame = true;

    std::unique_ptr<double[]> x1, y1, z1;
    std::unique_ptr<double[]> x2, y2, z2;

    double rmsvalue(std::shared_ptr<Frame> &frame);

    void save_frame_coord(double x[], double y[], double z[], const std::shared_ptr<Frame> &frame);

    AmberMask mask_for_superpose;
    AmberMask mask_for_rmscalc;

    PBCUtils::MolPair mols;

    std::set<std::shared_ptr<Atom>> atoms_for_superpose;
    std::set<std::shared_ptr<Atom>> atoms_for_rmscalc;

    std::unique_ptr<XTCWriter> writer;
    void update(const std::shared_ptr<Frame> &frame);

    void save_superposed_frame(double *x, double *y, double *z, const std::shared_ptr<Frame> &frame);

    void saveJson(std::ostream &os) const;
};

#endif // TINKER_RMSDCAL_HPP
