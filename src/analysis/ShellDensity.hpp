//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_SHELLDENSITY_HPP
#define TINKER_SHELLDENSITY_HPP

#include <memory>
#include <unordered_set>
#include <string>
#include <map>
#include <list>
#include <utility>

#include "common.hpp"
#include "BasicAnalysis.hpp"
#include "atom.hpp"

class Frame;


class ShellDensity : public BasicAnalysis {
public:
    ShellDensity() { enable_outfile = true; }

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    static const std::string title() {
        return "Shell Density function";
    }

protected:
    Atom::AtomIndenter ids1;
    Atom::AtomIndenter ids2;

    std::unordered_set<std::shared_ptr<Atom>> group1;
    std::unordered_set<std::shared_ptr<Atom>> group2;

    double distance_width;

    int distance_bins;

    std::map<int, std::size_t> hist;

    std::size_t nframe = 0;

};


#endif //TINKER_SHELLDENSITY_HPP
