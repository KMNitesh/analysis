//
// Created by xiamr on 8/26/19.
//

#ifndef TINKER_CLUSTERVOLUME_HPP
#define TINKER_CLUSTERVOLUME_HPP

#include "std.hpp"
#include <boost/multi_array.hpp>
#include "BasicAnalysis.hpp"
#include "atom.hpp"
#include "HBond.hpp"

class Frame;

class ClusterVolume : public BasicAnalysis {
public:
    ClusterVolume();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    static std::string title() { return "vdW Volume"; }

protected:

    enum class ATOM_Category : uint8_t {
        EMPTY, DEST, OTHER
    };

    static double getVdwRadii(const std::shared_ptr<Atom> &atom);

    void fill_atom(double grid_x_step, double grid_y_step, double grid_z_step,
                   const std::shared_ptr<Atom> &atom, ATOM_Category category);

    std::vector<std::tuple<int, int, int>>
    generate_neighbor_grids(double grid_x_step, double grid_y_step, double grid_z_step) const;

    AmberMask atom_mask;

    std::unordered_set<std::shared_ptr<Atom>> atom_group, other_atoms;

    int grid_x;
    int grid_y;
    int grid_z;

    std::deque<std::tuple<double, double, double, double>> volumes;


    inline static std::unordered_map<Symbol, double> vdWRadiis{
            {Symbol::Hydrogen,   1.10},
            {Symbol::Carbon,     1.70},
            {Symbol::Nitrogen,   1.55},
            {Symbol::Oxygen,     1.52},
            {Symbol::Phosphorus, 1.80},
            {Symbol::Sulfur,     1.80}
    };

    size_t countFilledGridPoints() const;

    bool fill_space(double grid_x_step, double grid_y_step, double grid_z_step);

    boost::multi_array<ATOM_Category, 3> grid;
};


#endif //TINKER_CLUSTERVOLUME_HPP
