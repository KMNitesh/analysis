//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_SHELLDENSITY_HPP
#define TINKER_SHELLDENSITY_HPP

#include <list>
#include <map>
#include <memory>
#include <string>
#include <unordered_set>
#include <utility>

#include "AbstractAnalysis.hpp"
#include "data_structure/atom.hpp"
#include "dsl/AmberMask.hpp"

class Frame;

class ShellDensity : public AbstractAnalysis {
public:
    ShellDensity();

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    [[nodiscard]] static std::string_view title() { return "Shell Density function"; }

protected:
    AmberMask mask1;
    AmberMask mask2;

    std::unordered_set<std::shared_ptr<Atom>> group1;
    std::unordered_set<std::shared_ptr<Atom>> group2;

    double distance_width;

    int distance_bins;

    std::map<int, std::size_t> hist;

    std::size_t nframe = 0;

    void saveJson(std::ostream &os) const;
};

#endif  // TINKER_SHELLDENSITY_HPP
