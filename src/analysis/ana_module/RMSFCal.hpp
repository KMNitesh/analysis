//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_RMSFCAL_HPP
#define TINKER_RMSFCAL_HPP

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <map>
#include <memory>
#include <string>
#include <unordered_set>

#include "AbstractAnalysis.hpp"
#include "ana_module/RMSDCal.hpp"
#include "data_structure/atom.hpp"
#include "dsl/AmberMask.hpp"
#include "utils/PBCUtils.hpp"
#include "utils/common.hpp"

class Frame;

class RMSFCal : public AbstractAnalysis {
public:
    RMSFCal();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    [[nodiscard]] static std::string_view title() { return "RMSF Calculator"; }

protected:
    std::unique_ptr<std::ostream> pdb_ostream;

    AmberMask mask_for_superpose;
    AmberMask mask_for_rmsfcalc;

    std::vector<std::shared_ptr<Atom>> atoms_for_superpose;
    std::vector<std::shared_ptr<Atom>> atoms_for_superpose_and_rmsfcalc;
    std::vector<std::shared_ptr<Atom>> atoms_for_rmsfcalc;

    PBCUtils::MolPair mols;

    AmberMask mask_for_first_frame_output;
    std::vector<std::shared_ptr<Atom>> atoms_for_first_frame_output;

    bool output_residue_average = false;

    struct Residue {
        Residue(std::size_t num, std::string name, std::shared_ptr<Molecule> mol)
            : num(num), name(std::move(name)), mol(std::move(mol)) {}
        std::size_t num;
        std::string name;
        std::shared_ptr<Molecule> mol;
        std::set<std::shared_ptr<Atom>> atoms;
    };
    std::vector<Residue> residues;

    int steps = 0;  // current frame number

    static double rmsfit(double x1[], double y1[], double z1[], double x2[], double y2[], double z2[], int nfit) {
        return RMSDCal::rmsfit(x1, y1, z1, x2, y2, z2, nfit);
    }

    static void jacobi(int n, double a[4][4], double d[], double v[4][4]) { RMSDCal::jacobi(n, a, d, v); }

    static void quatfit(int n1, double x1[], double y1[], double z1[], int n2, double x2[], double y2[], double z2[],
                        int nfit) {
        RMSDCal::quatfit(n1, x1, y1, z1, n2, x2, y2, z2, nfit);
    }

    static void center(int n_for_center, double x[], double y[], double z[], double mid[], int nfit);

    std::vector<double> x1, y1, z1;
    std::vector<double> x2, y2, z2;

    [[nodiscard]] double rmsvalue(int index) const;

    std::vector<std::array<
        boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::variance>>,
        3>>
        acc;

    void append_pdb();

    void save_coord(double *x, double *y, double *z, const std::shared_ptr<Frame> &frame);

    void update_accumulators(double x[], double y[], double z[], int nfit1, int nfit2, int n);

    void allocate_array_memory();

    void saveJson(std::ostream &os) const;

private:
    void add_atom_to_residue_vector(std::shared_ptr<Atom> &atom);
};

#endif  // TINKER_RMSFCAL_HPP
