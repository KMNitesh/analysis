//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_GREENKUBO_HPP
#define TINKER_GREENKUBO_HPP

#include <map>
#include <memory>
#include <string>
#include <unordered_set>

#include "AbstractAnalysis.hpp"
#include "data_structure/atom.hpp"
#include "dsl/AmberMask.hpp"

class Frame;

// Use Green-Kubo equation to calculate self-diffuse coefficients
class GreenKubo : public AbstractAnalysis {
public:
    GreenKubo();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    [[nodiscard]] static std::string_view title() { return "Green-Kubo"; }

private:
    AmberMask ids;
    std::unordered_set<std::shared_ptr<Atom>> group;

    double timestep;
    std::size_t steps = 0;

    std::map<int, double> vecx_map;
    std::map<int, double> vecy_map;
    std::map<int, double> vecz_map;
};

#endif // TINKER_GREENKUBO_HPP
