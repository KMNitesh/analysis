#ifndef TINKER_PROTEINDIHEDRAL_HPP
#define TINKER_PROTEINDIHEDRAL_HPP

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

#include "AbstractAnalysis.hpp"
#include "data_structure/atom.hpp"
#include "utils/PBCUtils.hpp"
#include "utils/std.hpp"



class Frame;

class ProteinDihedral : public AbstractAnalysis {
public:
    ProteinDihedral();

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void readInfo() override;

    [[nodiscard]] static std::string_view title() { return "Protein Dihedral Angle Analysis"; }

protected:

    AmberMask mask, init_mask;
    std::vector<std::shared_ptr<Atom>> atom_sequence;

    std::vector<boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::variance>>> dihedrals;

};

#endif // TINKER_PROTEINDIHEDRAL_HPP