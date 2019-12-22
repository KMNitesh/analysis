//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_RMSFCAL_HPP
#define TINKER_RMSFCAL_HPP

#include <memory>
#include <string>
#include <map>
#include <unordered_set>

#include "utils/common.hpp"
#include "AbstractAnalysis.hpp"
#include "data_structure/atom.hpp"
#include "RMSDCal.hpp"

class Frame;

class RMSFCal : public AbstractAnalysis {
public:

    RMSFCal();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    [[nodiscard]] static std::string_view title() { return "RMSF Calculator"; }

    ~RMSFCal() override;

protected:

    std::unique_ptr<std::ostream> pdb_ostream;

    AmberMask mask_for_superpose;
    AmberMask mask_for_rmsfcalc;

    std::vector<std::shared_ptr<Atom>> atoms_for_superpose;
    std::vector<std::shared_ptr<Atom>> atoms_for_superpose_and_rmsfcalc;
    std::vector<std::shared_ptr<Atom>> atoms_for_rmsfcalc;

    AmberMask mask_for_first_frame_output;
    std::vector<std::shared_ptr<Atom>> atoms_for_first_frame_output;

    int steps = 0; // current frame number

    static double rmsfit(double x1[], double y1[], double z1[],
                         double x2[], double y2[], double z2[], int nfit) {
        return RMSDCal::rmsfit(x1, y1, z1, x2, y2, z2, nfit);
    }

    static void jacobi(int n, double a[4][4], double d[], double v[4][4]) {
        RMSDCal::jacobi(n, a, d, v);
    }

    static void quatfit(int n1, double x1[], double y1[], double z1[],
                        int n2, double x2[], double y2[], double z2[], int nfit) {
        RMSDCal::quatfit(n1, x1, y1, z1, n2, x2, y2, z2, nfit);
    }

    static void center(int n_for_center, double x[], double y[], double z[], double mid[], int nfit);

    double *x_avg = nullptr, *y_avg = nullptr, *z_avg = nullptr;
    double *x1 = nullptr, *y1 = nullptr, *z1 = nullptr;
    double *x2 = nullptr, *y2 = nullptr, *z2 = nullptr;

    [[nodiscard]] double rmsvalue(int index) const;

    std::vector<std::vector<std::tuple<double, double, double>>> coords;

    void append_pdb(const std::vector<std::tuple<double, double, double>> &f_coord);

    std::tuple<double, double, double> save_coord(double *x, double *y, double *z, const std::shared_ptr<Frame> &frame);

    std::vector<std::tuple<double, double, double>>
    append_coord(double x[], double y[], double z[], int nfit1, int nfit2, int n);

    void calculate_average_structure();

    void allocate_array_memory();

    void saveJson(std::ostream &os) const;
};

#endif //TINKER_RMSFCAL_HPP
