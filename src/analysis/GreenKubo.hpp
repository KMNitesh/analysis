//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_GREENKUBO_HPP
#define TINKER_GREENKUBO_HPP

#include <memory>
#include <unordered_set>
#include <string>
#include <map>

#include "common.hpp"
#include "AbstractAnalysis.hpp"
#include "atom.hpp"

class Frame;


// Use Green-Kubo equation to calculate self-diffuse coefficients
class GreenKubo : public AbstractAnalysis {

public:
    GreenKubo() {
        enable_read_velocity = true;
        enable_tbb = true;
        enable_outfile = true;
    }

    void process(std::shared_ptr<Frame> &) override;

    void print(std::ostream &os) override;

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void readInfo() override;

    static const std::string title() {
        return "Green-Kubo";
    }

private:
    Atom::AmberMask ids;
    std::unordered_set<std::shared_ptr<Atom>> group;

    double timestep;
    std::size_t steps = 0;

    std::map<int, double> vecx_map;
    std::map<int, double> vecy_map;
    std::map<int, double> vecz_map;
};

#endif //TINKER_GREENKUBO_HPP
