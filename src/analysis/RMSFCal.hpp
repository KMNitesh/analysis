//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_RMSFCAL_HPP
#define TINKER_RMSFCAL_HPP

#include <memory>
#include <string>
#include <map>
#include <unordered_set>

#include "common.hpp"
#include "BasicAnalysis.hpp"
#include "atom.hpp"
#include "RMSDCal.hpp"

class Frame;

class RMSFCal : public BasicAnalysis {


    Atom::AtomIndenter ids;

    std::unordered_set<std::shared_ptr<Atom>> group;


    void find_matched_atoms(std::shared_ptr<Frame> &frame);

    int steps = 0; // current frame number
    bool first_frame = true;

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

    static void center(int n, double x[], double y[], double z[], double mid[], int nfit) {
        mid[0] = mid[1] = mid[2] = 0.0;
        double norm = 0.0;
        for (int i = 0; i < nfit; ++i) {
            mid[0] += x[i];
            mid[1] += y[i];
            mid[2] += z[i];
            norm += 1.0;
        }
        mid[0] /= norm;
        mid[1] /= norm;
        mid[2] /= norm;

        for (int i = 0; i < n; ++i) {
            x[i] -= mid[0];
            y[i] -= mid[1];
            z[i] -= mid[2];
        }
    }

    double x_avg[ATOM_MAX], y_avg[ATOM_MAX], z_avg[ATOM_MAX];
    double x1[ATOM_MAX], y1[ATOM_MAX], z1[ATOM_MAX];
    double x2[ATOM_MAX], y2[ATOM_MAX], z2[ATOM_MAX];

    double rmsvalue(int index);

    std::map<int, std::map<int, double>> x, y, z;

public:
    RMSFCal() {
        enable_outfile = true;
    }

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    static const std::string title() {
        return "RMSF Calculator";
    }
};

#endif //TINKER_RMSFCAL_HPP
