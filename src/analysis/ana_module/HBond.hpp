//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_HBOND_HPP
#define TINKER_HBOND_HPP

#include <map>
#include <memory>
#include <string>
#include <unordered_set>
#include <utility>

#include "AbstractAnalysis.hpp"
#include "data_structure/atom.hpp"
#include "dsl/AmberMask.hpp"

class Frame;

enum class Symbol { Hydrogen, Carbon, Nitrogen, Oxygen, Phosphorus, Sulfur, Sodium, X, Unknown };

Symbol which(const std::shared_ptr<Atom> &atom);

enum class Selector { Acceptor, Donor, Both };

enum class HBondType { VMDVerion, GMXVersion };

class HBond : public AbstractAnalysis {
public:
    HBond();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    [[nodiscard]] static std::string_view title() { return "Hydrogen Bond"; }

    void setParameters(const AmberMask &donor, const AmberMask &acceptor, double distance, double angle,
                       std::string criteria, const std::string &outfilename);

private:
    void Selector_Both(const std::shared_ptr<Frame> &frame);

    void Selector_Donor_Acceptor(const std::shared_ptr<Frame> &frame);

    bool check_hbond(const std::shared_ptr<Atom> &donor, const std::shared_ptr<Atom> &hydrogen,
                     const std::shared_ptr<Atom> &acceptor, const std::shared_ptr<Frame> &frame);

    void record_hbond(const std::shared_ptr<Atom> &donor, const std::shared_ptr<Atom> &acceptor,
                      const std::shared_ptr<Frame> &frame);

    AmberMask mask1, mask2;

    std::vector<std::array<std::shared_ptr<Atom>, 2>> donor_hydrogens;
    std::vector<std::shared_ptr<Atom>> acceptors;

    HBondType hbond_type = HBondType::VMDVerion;

    Selector mode;
    double donor_acceptor_dist_cutoff;
    double angle_cutoff;

    std::deque<std::map<std::array<std::shared_ptr<Atom>, 2>, double>> hbdist;

    std::deque<std::map<std::string, double>> residue_map;

    inline static const std::set<Symbol> donor_acceptor_symbols{Symbol::Oxygen, Symbol::Nitrogen, Symbol::Sulfur};

    void output_distance_statistics(std::ostream &os);
};

#endif  // TINKER_HBOND_HPP
