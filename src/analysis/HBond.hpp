//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_HBOND_HPP
#define TINKER_HBOND_HPP

#include <memory>
#include <unordered_set>
#include <string>
#include <map>

#include "common.hpp"
#include "AbstractAnalysis.hpp"
#include "atom.hpp"

class Frame;


enum class Symbol {
    Hydrogen,
    Carbon,
    Nitrogen,
    Oxygen,
    Phosphorus,
    Sulfur,
    Sodium,
    X,
    Unknown
};


Symbol which(const std::shared_ptr<Atom> &atom);

enum class Selector {
    Acceptor,
    Donor,
    Both
};

enum class HBondType {
    VMDVerion,
    GMXVersion
};

class HBond : public AbstractAnalysis {

public:
    HBond() {
        enable_outfile = true;
        enable_forcefield = true;
    }

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    static const std::string title() {
        return "Hydrogen Bond";
    }

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

private:
    std::map<int, int> hbonds;

    HBondType hbond_type = HBondType::VMDVerion;


    void Selector_Both(std::shared_ptr<Frame> &frame);

    void Selector_Donor(std::shared_ptr<Frame> &frame);

    void Selector_Acceptor(std::shared_ptr<Frame> &frame);

    Atom::AmberMask ids1;
    Atom::AmberMask ids2;

    std::unordered_set<std::shared_ptr<Atom>> group1;
    std::unordered_set<std::shared_ptr<Atom>> group2;

    Selector mode;
    double donor_acceptor_dist_cutoff;
    double angle_cutoff;
    int steps = 0;

};


#endif //TINKER_HBOND_HPP
