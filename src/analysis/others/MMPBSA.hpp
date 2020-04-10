
#ifndef TINKER_MMPBSA_HPP
#define TINKER_MMPBSA_HPP

#include "utils/std.hpp"

class MMPBSA {
public:
    [[nodiscard]] static std::string_view title() { return "Amber MMPBSA.py Decomposition data post-analysis"; }

    static void process(const std::string &topology_filename);

    struct Component {
        Component() = default;
        Component(std::initializer_list<std::string> list) {
            assert(list.size() == 3);
            auto it = std::begin(list);
            avg = std::stod(*it++);
            stddev = std::stod(*it++);
            stderr = std::stod(*it);
        }
        double avg = 0.0, stddev = 0.0, stderr = 0.0;

        const Component &operator+=(const Component &rhs) {
            avg += rhs.avg;
            stddev += rhs.stddev;
            stderr += rhs.stderr;
            return *this;
        }
    };

    struct ResidueComponent {
        Component internal, vdW, electrostatic, polar_solv, nonpolar_solv, total;

        const ResidueComponent &operator+=(const ResidueComponent &rhs) {
            internal += rhs.internal;
            vdW += rhs.vdW;
            electrostatic += rhs.electrostatic;
            polar_solv += rhs.polar_solv;
            nonpolar_solv += rhs.nonpolar_solv;
            total += rhs.total;
            return *this;
        }
    };

    struct Residue {
        enum TYPE { R, L } location;
        std::string name;
        uint no;
    };

    friend std::istream &operator>>(std::istream &, std::vector<std::pair<Residue, ResidueComponent>> &);

    friend std::ostream &operator<<(std::ostream &, const std::vector<std::pair<Residue, ResidueComponent>> &);
};

#endif // TINKER_MMPBSA_HPP