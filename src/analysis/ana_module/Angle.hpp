
#ifndef TINKER_ANGLE_HPP
#define TINKER_ANGLE_HPP

#include <Eigen/Eigen>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

#include "AbstractAnalysis.hpp"
#include "data_structure/atom.hpp"
#include "dsl/AmberMask.hpp"
#include "utils/PBCUtils.hpp"
#include "utils/std.hpp"

class Frame;
class IntertiaVector;

class Angle : public AbstractAnalysis {
public:
    Angle();

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void readInfo() override;

    [[nodiscard]] static std::string_view title() { return "Angle between two groups (mass-weighted)"; }

    enum class AxisType { MIN, MAX };

    void setParameters(const IntertiaVector &v1, const IntertiaVector &v2, const std::string &out);

protected:
    Eigen::Matrix3d calculate_inertia(std::shared_ptr<Frame> &frame,
                                      const std::vector<std::shared_ptr<Atom>> &atom_group, PBCUtils::MolPair &mols);

    [[nodiscard]] std::tuple<double, double, double> calculate_axis(AxisType type,
                                                                    const Eigen::Matrix3d &inertia) const;

    AmberMask mask1, mask2;

    AxisType type1, type2;

    std::vector<std::shared_ptr<Atom>> atom_group1, atom_group2;

    PBCUtils::MolPair mol1, mol2;

    std::deque<double> deque;
};

#endif // TINKER_ANGLE_HPP
