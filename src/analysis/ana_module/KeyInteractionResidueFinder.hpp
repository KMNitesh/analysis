
#ifndef KEYINTERACTIONRESIDUEFINDER_HPP
#define KEYINTERACTIONRESIDUEFINDER_HPP

#include "ana_module/AbstractAnalysis.hpp"
#include "data_structure/atom.hpp"
#include "dsl/AmberMask.hpp"

class Frame;

class KeyInteractionResidueFinder : public AbstractAnalysis {
public:
    KeyInteractionResidueFinder();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    [[nodiscard]] static std::string_view title() { return "Key Residue Finder for Protein-DNA interaction"; }

private:
    AmberMask sidechain_N_O_mask, dna_P_mask;
    double cutoff;

    std::vector<std::shared_ptr<Atom>> sidechain_N_O_atoms, dna_P_atoms;

    std::deque<std::map<std::string, double>> residue_distance_evolution;
};

#endif // KEYINTERACTIONRESIDUEFINDER_HPP